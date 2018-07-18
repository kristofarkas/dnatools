import os
import tempfile
import subprocess

import numpy as np

import parmed as pmd

from openeye.oechem import *
from openeye.oedepict import *
from openeye.oedocking import *
from openeye.oeiupac import OECreateIUPACName

# from yank.utils import TLeap


class Intercalator:

    atom_expression = OEExprOpts_Aromaticity | OEExprOpts_AtomicNumber | \
                      OEExprOpts_RingMember | OEExprOpts_StringType
    bond_expression = OEExprOpts_Aromaticity | OEExprOpts_BondOrder | OEExprOpts_RingMember

    def __init__(self, mol=None):
        self._mol: OEMCMolBase = mol or OEMol()
        self._frcmod: pmd.amber.AmberParameterSet = None
        self._mcss: OEMCSSearch = None

    @property
    def mol(self):
        return self._mol

    @mol.setter
    def mol(self, new_mol):
        self._mol = new_mol

    @property
    def mcss(self):
        return self._mcss

    @classmethod
    def from_smiles(cls, smiles, assign_aromaticity=True, add_hydrogens=True):
        i = Intercalator()
        m = i._mol
        OESmilesToMol(m, smiles)

        i.sanitize(assign_aromaticity, add_hydrogens)

        return i

    def sanitize(self, assign_aromaticity=True, add_hydrogens=True):
        if assign_aromaticity:
            OEAssignAromaticFlags(self.mol, OEAroModelOpenEye)

        if add_hydrogens:
            OEAddExplicitHydrogens(self.mol)

        self.mol.SetTitle(OECreateIUPACName(self.mol))

        # TODO: only assign TRIPOS names and types if we want to write this molecule as MOL2 and
        # there are no names and types already assign. Also we would probably want to make these
        # changes only on a copy of the molecule, not the molecule itself!

        OETriposAtomNames(self.mol)
        OETriposAtomTypeNames(self.mol)
        OETriposBondTypeNames(self.mol)

    def prepare_for_mcs(self, approximate=False):
        self._mcss = OEMCSSearch(self._mol, self.atom_expression, self.bond_expression,
                                 OEMCSType_Default if not approximate else OEMCSType_Approximate)
        self._mcss.SetMCSFunc(OEMCSMaxAtomsCompleteCycles())
        self._mcss.SetMinAtoms(6)

    def total_charge(self):
        charge = 0.
        for atom in self._mol.GetAtoms():
            charge += atom.GetPartialCharge()
        return charge

    def assign_amber_atomtypes(self):
        with tempfile.TemporaryDirectory() as tempdir:
            name = os.path.join(tempdir, 'mol.mol2')
            out_name = os.path.join(tempdir, 'gaff.mol2')

            self.write_mol2(name)
            subprocess.run(['antechamber',
                            '-i', name, '-fi', 'mol2',
                            '-o', out_name, '-fo', 'mol2',
                            '-at', 'gaff2', '-du', 'n', '-an', 'n', '-j', '1', '-pf', 'y', '-dr', 'n'])
            gaff = Intercalator.read_mol2(out_name)
            for a, b in zip(self._mol.GetAtoms(), gaff.GetAtoms()):
                a.SetType(b.GetType())

    def generate_conformers(self, max_confs=800, rms_threshold=1.0):

        self._mol.SetDimension(3)

        from openeye.oeomega import OEOmega

        om = OEOmega()
        om.SetMaxConfs(max_confs)
        om.SetIncludeInput(True)
        om.SetStrictStereo(False)
        om.SetCanonOrder(False)
        om.SetSampleHydrogens(True)
        om.SetEnergyWindow(15.0)
        om.SetRMSThreshold(rms_threshold)

        success = om(self._mol)

        if not success:
            print('Failed to generate multiple conformers from molecule.')

    def assign_am1bcc_charges(self):
        from openeye.oequacpac import OEAM1BCCELF10Charges, OEAssignCharges

        am1bcc = OEAM1BCCELF10Charges()
        am1bcc.SetReturnSelectedConfs(True)
        OEAssignCharges(self._mol, am1bcc)

    @property
    def frcmod(self):

        if self._frcmod is not None:
            return self._frcmod

        with tempfile.TemporaryDirectory() as tempdir:
            name = os.path.join(tempdir, 'mol.mol2')
            out_name = os.path.join(tempdir, 'mol.frcmod')

            self.write_mol2(name)
            subprocess.run(['parmchk2', '-i', name, '-f', 'mol2', '-o', out_name, '-s', '2'])

            self._frcmod = pmd.load_file(out_name)

        return self._frcmod

    def write_mol2(self, name):
        ofs = oemolostream(name)
        ofs.SetFlavor(OEFormat_MOL2, OEOFlavor_MOL2_Forcefield)
        OEWriteMol2File(ofs, self._mol)

    @staticmethod
    def read_mol2(name):
        mol = OEMol()
        ifs = oemolistream(name)
        ifs.SetFlavor(OEFormat_MOL2, OEIFlavor_MOL2_Forcefield)
        OEReadMolecule(ifs, mol)
        return mol


class HybridIntercalator(Intercalator):

    def __init__(self, pattern: Intercalator, target: Intercalator, threshold=0.01):
        super(HybridIntercalator, self).__init__()
        self.target = target
        self.threshold = threshold

        self.pattern = pattern
        self.pattern.prepare_for_mcs()

        matches = self.pattern.mcss.Match(self.target.mol, True)
        self._match = matches.next()
        self._filter_match()
        self._align_molecules()
        self._create_hybrid()

    def _align_molecules(self):
        overlay = True
        rotate, translate = OEDoubleArray(9), OEDoubleArray(3)
        OERMSD(self.pattern.mcss.GetPattern(), self.target.mol, self._match, overlay, rotate, translate)
        OERotate(self.target.mol, rotate)
        OETranslate(self.target.mol, translate)

    def _filter_match(self):
        match = OEMatch()

        for m in self._match.GetAtoms():

            ratio = m.pattern.GetPartialCharge() / m.target.GetPartialCharge()
            max_change = max(abs(0.5 - 0.5 * ratio), abs(0.5 - 0.5 / ratio))

            if max_change < self.threshold:
                match.AddPair(m.pattern, m.target)

        for m in self._match.GetBonds():
            b0 = m.pattern.GetBgn() in match.GetPatternAtoms()
            b1 = m.pattern.GetEnd() in match.GetPatternAtoms()
            b2 = m.target.GetBgn() in match.GetTargetAtoms()
            b3 = m.target.GetEnd() in match.GetTargetAtoms()
            if all([b0, b1, b2, b3]):
                match.AddPair(m.pattern, m.target)

        self._match = match

    def _create_hybrid(self):

        # Using GraphMol as this will have one conformation only.

        hybrid_atom_map = OEAtomArray(self.pattern.mol.GetMaxAtomIdx())
        hybrid = OEGraphMol()
        OECopyMol(hybrid, self.pattern.mol, hybrid_atom_map)

        # Add atoms from target

        atom_map = {}

        for atom in self.target.mol.GetAtoms():
            if atom not in self._match.GetTargetAtoms():
                atom_map[atom] = hybrid.NewAtom(atom)

        # Add bonds from target

        for bond in self.target.mol.GetBonds():
            if bond not in self._match.GetTargetBonds():
                bgn = bond.GetBgn()
                end = bond.GetEnd()

                if (bgn not in self._match.GetTargetAtoms()) and (end not in self._match.GetTargetAtoms()):
                    b = hybrid.NewBond(atom_map[bgn], atom_map[end], bond.GetOrder())
                    b.SetType(bond.GetType())
                elif bgn not in self._match.GetTargetAtoms():
                    pat = next(pair.pattern for pair in self._match.GetAtoms() if pair.target == end)
                    h_end = hybrid.GetAtom(OEHasAtomIdx(pat.GetIdx()))
                    b = hybrid.NewBond(h_end, atom_map[bgn], bond.GetOrder())
                    b.SetType(bond.GetType())
                elif end not in self._match.GetTargetAtoms():
                    pat = next(pair.pattern for pair in self._match.GetAtoms() if pair.target == bgn)
                    h_bgn = hybrid.GetAtom(OEHasAtomIdx(pat.GetIdx()))
                    b = hybrid.NewBond(h_bgn, atom_map[end], bond.GetOrder())
                    b.SetType(bond.GetType())

        # Average charges and rescale appearing/disappearing part by charge adjustment factor.

        charge_adjustment = 0

        for p in self._match.GetAtoms():
            new_charge = np.mean([p.pattern.GetPartialCharge(), p.target.GetPartialCharge()])
            charge_adjustment += (p.pattern.GetPartialCharge() - p.target.GetPartialCharge()) / 2
            hybrid_atom_map[p.pattern.GetIdx()].SetPartialCharge(new_charge)

        for hybrid_atom in hybrid.GetAtoms():
            if hybrid_atom in atom_map.values():
                hybrid_atom.SetPartialCharge(hybrid_atom.GetPartialCharge() - charge_adjustment / len(atom_map))
            if (hybrid_atom in hybrid_atom_map) and not (
                    hybrid_atom in [hybrid_atom_map[a.GetIdx()] for a in self._match.GetPatternAtoms()]):
                hybrid_atom.SetPartialCharge(hybrid_atom.GetPartialCharge() + charge_adjustment / (
                            self.pattern.mol.NumAtoms() - len([_ for _ in self._match.GetPatternAtoms()])))

        # Check charges

        target_charge = 0
        pattern_charge = 0

        for a in hybrid.GetAtoms():
            if (a in atom_map.values()) or (a in [hybrid_atom_map[a.GetIdx()] for a in self._match.GetPatternAtoms()]):
                target_charge += a.GetPartialCharge()
            if a in hybrid_atom_map:
                pattern_charge += a.GetPartialCharge()

        self._mol = hybrid

    def draw_hybrid(self, name):

        p, t = self.pattern.mcss.GetPattern(), self.target.mol

        OEPrepareDepiction(p, False, False)
        OEPrepareDepiction(t, False, False)

        image = OEImage(1000, 500)

        rows, cols = 1, 2
        grid = OEImageGrid(image, rows, cols)

        opts = OE2DMolDisplayOptions(grid.GetCellWidth(), grid.GetCellHeight(), OEScale_AutoScale)
        opts.SetTitleLocation(OETitleLocation_Hidden)

        refscale = OEGetMoleculeScale(p, opts)
        fitscale = OEGetMoleculeScale(t, opts)
        opts.SetScale(min(refscale, fitscale))
        opts.SetHydrogenStyle(OEHydrogenStyle_ExplicitAll)

        refdisp = OE2DMolDisplay(p, opts)
        fitdisp = OE2DMolDisplay(t, opts)

        refabset = oechem.OEAtomBondSet(self._match.GetPatternAtoms(), self._match.GetPatternBonds())
        OEAddHighlighting(refdisp, oechem.OEBlueTint, OEHighlightStyle_BallAndStick, refabset)

        fitabset = oechem.OEAtomBondSet(self._match.GetTargetAtoms(), self._match.GetTargetBonds())
        OEAddHighlighting(fitdisp, oechem.OEBlueTint, OEHighlightStyle_BallAndStick, fitabset)

        refcell = grid.GetCell(1, 1)
        OERenderMolecule(refcell, refdisp)

        fitcell = grid.GetCell(1, 2)
        OERenderMolecule(fitcell, fitdisp)

        OEWriteImage(name, image)

    @property
    def charge_correction(self):
        diff = 0

        for m in self._match.GetAtoms():
            diff += (m.pattern.GetPartialCharge() - m.target.GetPartialCharge()) / 2

        return diff


class DNA:
    def __init__(self, dna: OEMolBase, box: OEBox):
        self._receptor: OEMolBase = OEGraphMol()
        OEMakeReceptor(self._receptor, dna, box)

        self._dock: OEDock = OEDock(OEDockMethod_Chemgauss4, OESearchResolution_High)
        self._dock.Initialize(self._receptor)

    @classmethod
    def with_original_intercalator(cls, dna: OEMolBase, original_intercalator: Intercalator, padding: float=2.0):

        box = OEBox(original_intercalator.mol)

        max_side = max(cls.get_box_dims(box))

        extend = [(max_side-c)/2 + padding for c in cls.get_box_dims(box)]

        OEBoxExtend(box, *extend)

        return cls(dna, box)

    @classmethod
    def with_intercalators(cls, dna: OEMolBase, original_intercalator: Intercalator, *intercalators: Intercalator,
                           padding: float=2.0):
        """Convenience initializer for multiple intercalators and when there was an original intercalator
        inside the crystal structure.

        :param dna: The DNA structure only.
        :param original_intercalator: The original intercalator that was in the crystal structure.
        :param intercalators: The list of intercalator you will want to dock
        :param padding: padding to add to the box created around the original intercalator
        :return: a DNA instance correctly initialized for docking
        """
        # This is the box of the original intercalator. We want to dock somewhere around here.
        box = OEBox(original_intercalator.mol)

        max_side = max(d for xyz in [cls.get_box_dims(OEBox(mol.mol)) for mol in intercalators + [original_intercalator]] for d in xyz)

        extend = [(max_side-c)/2 + padding for c in cls.get_box_dims(box)]

        OEBoxExtend(box, *extend)

        return cls(dna, box)

    @staticmethod
    def get_box_dims(box):
        return OEBoxXDim(box), OEBoxYDim(box), OEBoxZDim(box)

    def dock_intercalator(self, intercalator: Intercalator):
        """
        Dock molecule into DNA base-pair sequence.
        :param intercalator: molecule to be docked
        :return: docked molecule, rmsd
            Note, rmsd might only be relevant if the original molecule was, let's say, at a crystal structure
            docked site. Otherwise the rmsd from some arbitrary position is irrelevant.
        """
        docked_mol = OEGraphMol()
        self._dock.DockMultiConformerMolecule(docked_mol, intercalator.mol)
        sd_tag = OEDockMethodGetName(OEDockMethod_Chemgauss4)
        OESetSDScore(docked_mol, self._dock, sd_tag)
        self._dock.AnnotatePose(docked_mol)

        return Intercalator(docked_mol), OERMSD(intercalator.mol, docked_mol)


def load_oe_molecules(file_path, molecule_idx=None):
    """Read one or more molecules from a file.
    Requires OpenEye Toolkit. Several formats are supported (including
    mol2, sdf and pdb).
    Parameters
    ----------
    file_path : str
        Complete path to the file on disk.
    molecule_idx : None or int, optional, default: None
        Index of the molecule on the file. If None, all of them are
        returned.
    Returns
    -------
    molecule : openeye.oechem.OEMol or list of openeye.oechem.OEMol
        The molecules stored in the file. If molecule_idx is specified
        only one molecule is returned, otherwise a list (even if the
        file contain only 1 molecule).
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    from openeye import oechem
    extension = os.path.splitext(file_path)[1][1:]  # Remove dot.

    # Open input file stream
    ifs = oechem.oemolistream()
    if extension == 'mol2':
        mol2_flavor = (oechem.OEIFlavor_Generic_Default |
                       oechem.OEIFlavor_MOL2_Default |
                       oechem.OEIFlavor_MOL2_Forcefield)
        ifs.SetFlavor(oechem.OEFormat_MOL2, mol2_flavor)
    if not ifs.open(file_path):
        oechem.OEThrow.Fatal('Unable to open {}'.format(file_path))

    # Read all molecules.
    molecules = []
    for mol in ifs.GetOEMols():
        molecules.append(oechem.OEMol(mol))

    # Select conformation of interest
    if molecule_idx is not None:
        return molecules[molecule_idx]

    return molecules


def load_oe_molecule(file_path):
    return load_oe_molecules(file_path, molecule_idx=0)
