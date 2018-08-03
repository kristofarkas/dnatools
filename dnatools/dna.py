from openeye.oechem import *
from openeye.oedocking import *

from .intercalator import Intercalator


class DNA:
    def __init__(self, dna: OEMolBase, box: OEBox):
        self._dna: OEMolBase = dna
        self._receptor: OEMolBase = OEGraphMol()
        OEMakeReceptor(self._receptor, dna, box)

        self._dock: OEDock = OEDock(OEDockMethod_Chemgauss4, OESearchResolution_High)
        self._dock.Initialize(self._receptor)

    @property
    def dna(self):
        return self._dna

    @classmethod
    def with_original_intercalator(cls, dna: OEMolBase, original_intercalator: Intercalator, padding: float=2.0):

        box = OEBox(original_intercalator.mol)

        max_side = max(cls.get_box_dims(box))

        extend = [(max_side-c)/2 + padding for c in cls.get_box_dims(box)]

        OEBoxExtend(box, *extend)

        return cls(dna, box)

    @classmethod
    def with_intercalators(cls, dna: OEMolBase, original_intercalator: Intercalator, intercalators: [Intercalator],
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

        max_side = max(d for xyz in [cls.get_box_dims(OEBox(mol.mol)) for mol in (intercalators + [original_intercalator])] for d in xyz)

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

