class Topology:
    def __init__(self):
        pass


class HybridTopology(Topology):
    def __init__(self, mol, other):
        super(HybridTopology, self).__init__()
        self.mol = mol
        self.other = other

    def align_molecules(self):
        pass