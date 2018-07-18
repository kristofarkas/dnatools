import os
import tempfile
import subprocess

import numpy as np

import parmed as pmd

from openeye.oechem import *
from openeye.oedepict import *
from openeye.oeiupac import OECreateIUPACName


class Intercalator:

    atom_expression = OEExprOpts_Aromaticity | OEExprOpts_AtomicNumber | \
                      OEExprOpts_RingMember | OEExprOpts_StringType
    bond_expression = OEExprOpts_Aromaticity | OEExprOpts_BondOrder | OEExprOpts_RingMember

    def __init__(self):
        super(Intercalator, self).__init__()
        self._mol = OEMol()
        self._frcmod = None
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

        if assign_aromaticity:
            OEAssignAromaticFlags(m, OEAroModelOpenEye)

        if add_hydrogens:
            OEAddExplicitHydrogens(m)

        m.SetTitle(OECreateIUPACName(m))

        # TODO: only assign TRIPOS names and types if we want to write this molecule as MOL2 and
        # there are no names and types already assign. Also we would probably want to make these
        # changes only on a copy of the molecule, not the molecule itself!

        OETriposAtomNames(m)
        OETriposAtomTypeNames(m)
        OETriposBondTypeNames(m)

        return i

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
    def force_field_modification(self):

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


class HybridTopology:

    def __init__(self, pattern: Intercalator, target: Intercalator, threshold=0.01):
        super(HybridTopology, self).__init__()
        self.target = target
        self.threshold = threshold

        self.pattern = pattern
        self.pattern.prepare_for_mcs()

        self.hybrid = Intercalator()

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

        self.hybrid.mol = hybrid

    def draw_hybrid(self, name):

        p, t =  self.pattern.mcss.GetPattern(), self.target.mol

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
