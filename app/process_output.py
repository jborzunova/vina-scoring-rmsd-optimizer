from rdkit import Chem
from rdkit.Chem import AllChem
import os
from parameters import *
import pandas as pd
from logger import log 
from spyrmsd.molecule import Molecule
from spyrmsd import io, rmsd


def extract_first_pose(filename):
    # Extract the first molecule from the SDF file
    # Read the file and record the index marking the end of the ligand's first pose
    end_index = 0
    with open(filename) as h:
        lines = h.readlines()
        for ind, line in enumerate(lines):
            if end_index == 0 and line.find('$$$$') != -1:
                end_index = ind
    # write new file
    new_filename = filename[:-4] + '_1.sdf'
    with open(new_filename, 'w') as h:
        h.writelines(lines[:end_index+1])
    return new_filename


def pdbqt2sdf(lig_pdbqt):
    # Read pdbqt 
    lig_pdb = lig_pdbqt.replace('_prep_out.pdbqt', '.pdb')
    lig_sdf = lig_pdbqt.replace('.pdbqt', '.sdf')
    os.system(f'/home/jborzunova/.conda/envs/unidock/bin/python3 /home/jborzunova/.local/bin/mk_export.py {lig_pdbqt} -o {lig_sdf}')

    # extract the first pose
    first_pose_sdffile = extract_first_pose(lig_sdf)

    # get smiles string for proper bonds definition
    data = pd.read_csv(DF)
    smiles = data.query('file==@lig_pdb')['smiles'].iloc[0] 
   
    # assign bond orders
    suppl = Chem.SDMolSupplier(first_pose_sdffile)
    molecule = next((m for m in suppl if m is not None), None)  # We assume the SDF file has a single molecule (it was extracted earlier)
    reference_molecule = Chem.MolFromSmiles(smiles)
    molecule = AllChem.AssignBondOrdersFromTemplate(reference_molecule, molecule)
    lig_sdf_corr = first_pose_sdffile.replace('.sdf', '_correct_bonds.sdf')
    w = Chem.SDWriter(lig_sdf_corr)
    w.write(molecule)
    w.close()
    #log(f"Conversion successful! The sdf file is saved as {lig_sdf_corr}")
    return lig_sdf_corr


def get_sdffile_info(sdf_file):
    # Load the molecule from the SDF file
    suppl = Chem.SDMolSupplier(sdf_file)
    # Iterate through the molecules in the SDF file (there can be multiple, but there should be one)
    for rdkit_mol in suppl:
        if rdkit_mol is not None:  # Make sure the molecule is valid
            Chem.SanitizeMol(rdkit_mol)
            mol = Molecule.from_rdkit(rdkit_mol)
            mol.strip()
            coords = mol.coordinates
            anum = mol.atomicnums
            adj = mol.adjacency_matrix
        else:
            log("Warning: Failed to load molecule.")
    return coords, anum, adj


def calc_rmsd(lig_file1, lig_file2):
    # this function calculates rmsd between docking pose amd a native one
    coords_ref, anum_ref, adj_ref = get_sdffile_info(lig_file2)
    coords_mol, anum_mol, adj_mol = get_sdffile_info(lig_file1)
    rmsd_val = rmsd.symmrmsd(coords_ref, coords_mol,
                         anum_ref, anum_mol,
                         adj_ref, adj_mol)
    return rmsd_val


def rmsd_from_output(lig_pdbqt):
    lig_sdf = pdbqt2sdf(lig_pdbqt)
    ligand_native_file = lig_pdbqt.replace('_out.pdbqt', '.sdf')
    #log(f'I calculate RMSD between {lig_sdf} and {PATH_TO_NATIVE_LIGANDS+ligand_native_file}')
    RMSD = calc_rmsd(lig_sdf, PATH_TO_NATIVE_LIGANDS+ligand_native_file)
    log(f'RMSD is {RMSD}')
    return RMSD


if __name__ == '__main__':
    os.chdir('./uni2/workdir')
    # Input PDBQT file and desired output SDF file name
    rmsd_from_output('ligand_DLW_6FKJ_B_B_501_prep_out.pdbqt')

