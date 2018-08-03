import os
import subprocess

import parmed as pmd
from yank.utils import TLeap, write_oe_molecule

from .dna import DNA
from .intercalator import Intercalator, HybridIntercalator


class Complex:
    def __init__(self, sim_dir='.'):
        self.dna: DNA = None
        self.intercalator: Intercalator = None
        self.sim_dir = sim_dir

    _water_model = {
        'tip3': ("leaprc.water.tip3p", "TIP3PBOX"),
        'tip4': ("leaprc.water.tip4pew", "TIP4PEWBOX")
    }

    def prepare_for_simulation(self, water_model: str, clearance: float):
        """
        Prepare the complex for simulation.

        Create all the input files for NAMD, run TLeap solvate the complex and so on.
        :param water_model: One of ['tip3', 'tip4']
        :param clearance: distance of water box from end of complex, ing angstrom
        """

        dna = os.path.join(self.sim_dir, "dna.pdb")
        intercalator = os.path.join(self.sim_dir, "intercalator.mol2")
        frcmod = os.path.join(self.sim_dir, "intercalator.frcmod")

        intercalator_pdb = os.path.join(self.sim_dir, "intercalator.pdb")
        intercalator_prmtop = os.path.join(self.sim_dir, "intercalator.prmtop")

        pdb = os.path.join(self.sim_dir, "complex.pdb")
        prmtop = os.path.join(self.sim_dir, "complex.prmtop")

        write_oe_molecule(self.dna.dna, dna)
        self.intercalator.write_mol2(intercalator)
        self.intercalator.frcmod.write(frcmod)

        t = TLeap()

        if water_model == 'tip4':
            # TLeap won't write water parameters unless this is on.
            # NAMD needs those parameters for TIP4P support.
            # Nonetheless NAMD has bug (in some version(s)) that fail with tip4p water.
            t.add_commands("set default FlexibleWater on")
        t.load_parameters("leaprc.DNA.bsc1")
        t.load_parameters("leaprc.gaff2")
        t.load_parameters(self._water_model[water_model][0])

        # Load intercalator
        t.load_unit("intercalator", intercalator)
        t.load_parameters(frcmod)

        # Load dna
        t.load_unit("dna", dna)

        # Combine, solvate, ions, save.
        t.combine("complex", "intercalator", "dna")

        t.solvate("intercalator", self._water_model[water_model][1], clearance)
        t.add_ions("intercalator", "Na+")

        t.save_unit("intercalator", intercalator_prmtop)
        t.save_unit("intercalator", intercalator_pdb)

        t.solvate("complex", self._water_model[water_model][1], clearance)
        t.add_ions("complex", "Na+")

        t.save_unit("complex", prmtop)
        t.save_unit("complex", pdb)

        t.run()

        if isinstance(self.intercalator, HybridIntercalator):

            pdbs = {"intercalator_tags.pdb": intercalator_pdb,
                    "tags.pdb": pdb}

            for name, pdb in pdbs.items():
                dn: pmd.structure.Structure = pmd.load_file(pdb)
                for atom, alchemical_tag in zip(dn.residues[0].atoms, self.intercalator.alchemical_tag):
                    atom.bfactor = alchemical_tag

                tags = os.path.join(self.sim_dir, name)

                dn.write_pdb(tags, renumber=False, increase_tercount=False)
                subprocess.run(["sed", "-i.bak",  "s/HETATM/ATOM  /g", tags])

        return pdb, prmtop
