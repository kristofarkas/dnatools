import os


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