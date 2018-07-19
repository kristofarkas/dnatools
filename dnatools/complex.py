import os

import parmed as pmd
from yank.utils import TLeap, write_oe_molecule

from .topology import DNA, Intercalator, HybridIntercalator


class Complex:
    def __init__(self, sim_dir='.'):
        self.dna: DNA = None
        self.intercalator: Intercalator = None
        self.sim_dir = sim_dir

    def prepare_for_simulation(self):

        dna = os.path.join(self.sim_dir, "dna.pdb")
        intercalator = os.path.join(self.sim_dir, "intercalator.mol2")
        frcmod = os.path.join(self.sim_dir, "intercalator.frcmod")

        pdb = os.path.join(self.sim_dir, "complex.pdb")
        prmtop = os.path.join(self.sim_dir, "complex.prmtop")

        write_oe_molecule(self.dna.dna, dna)
        self.intercalator.write_mol2(intercalator)
        self.intercalator.frcmod.write(frcmod)

        t = TLeap()

        t.add_commands("set default FlexibleWater on")
        t.load_parameters("leaprc.DNA.bsc1")
        t.load_parameters("leaprc.gaff2")
        t.load_parameters("leaprc.water.tip4pew")

        # Load intercalator
        t.load_unit("intercalator", intercalator)
        t.load_parameters(frcmod)

        # Load dna
        t.load_unit("dna", dna)

        # Combine, solvate, ions, save.
        t.combine("complex", "intercalator", "dna")
        t.solvate("complex", "TIP4PEWBOX", 10.0)
        t.add_ions("complex", "Na+")

        t.save_unit("complex", prmtop)
        t.save_unit("complex", pdb)

        t.run()

        if isinstance(self.intercalator, HybridIntercalator):
            dn: pmd.structure.Structure = pmd.load_file(pdb)
            for atom, alchemical_tag in zip(dn.residues[0].atoms, self.intercalator.alchemical_tag):
                atom.bfactor = alchemical_tag

            tags = os.path.join(self.sim_dir, "tags.pdb")

            dn.write_pdb(tags, renumber=False, increase_tercount=False)

        return pdb, prmtop
