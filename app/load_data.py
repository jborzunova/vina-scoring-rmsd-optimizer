import os
import pandas as pd
from parameters import *


def define_ref_ligands(receptor, data):
    txt = receptor.replace('protein', 'ligand').replace('_esswat_prep.pdbqt', '.pdb') 
    N = int(data.query('file == @txt')['n_cluster'].item()) # the cluster number (of ligand)
    ligands_pdb = list(data.query('n_cluster == @N')['file'])
    ligands = [file.replace('.pdb', '_prep.pdbqt') for file in ligands_pdb]
    return N, ligands  # N is for folder labels

def define_fold_ligands(data):
    N = 200
    ligands_pdb = list(data.query(f'fold_{FOLD_NUMBER} == "train"')['file'])
    ligands = [file.replace('.pdb', '_prep.pdbqt') for file in ligands_pdb]
    return N, ligands  # N is for further sorting

def define_all_ligands(data):
    N = 200
    ligands_pdb = list(data.query(f'cl_center == True')['file'])
    ligands = [file.replace('.pdb', '_prep.pdbqt') for file in ligands_pdb]
    return N, ligands  # N is for further sorting


def load_pairs(df_path=DF):
    print(f"Reading: {df_path}")
    data = pd.read_csv(df_path)
    LIG_PROTEIN_PAIRS = []
    """Loads ligand-protein pairs from the given path."""
    if MODE == 0:  # ref docking of train
        data_proteins = data.query(f'cl_center')
        receptors = ["protein" + f[6:-4] + "_esswat_prep.pdbqt" for f in data_proteins["file"].to_list()]
        unsorted_pairs = []
        for receptor in receptors:
            n, ligands = define_ref_ligands(receptor, data)
            unsorted_pairs.append((n, ligands, receptor))
        # Sort by cluster number
        LIG_PROTEIN_PAIRS = sorted(unsorted_pairs, key=lambda x: x[0])
        # Drop the integer and keep only the ligand list and protein file
    elif MODE == 1:  # one file
        receptor = 'protein_WNY_5S4P_B_B_505_esswat_prep.pdbqt'
        n, ligands = define_ref_ligands(receptor, data)
        LIG_PROTEIN_PAIRS = [(n, ligands, receptor)]  
    elif MODE == 2:  # 2 files
        receptors = ['protein_YJ7_7LZ7_B_B_501_esswat_prep.pdbqt', 'protein_WNY_5S4P_B_B_505_esswat_prep.pdbqt']
        for receptor in receptors:
            n, ligands = define_ref_ligands(receptor, data)
            LIG_PROTEIN_PAIRS.append((n, ligands, receptor))
    elif MODE == 3:
        receptor = 'protein_89C_5XKH_B_B_504_esswat_prep.pdb'
        n, ligands = define_fold_ligands(data)
        LIG_PROTEIN_PAIRS = [(n, ligands, receptor)]
    elif MODE == 4:
        receptor = 'protein_89C_5XKH_B_B_504_esswat_prep.pdb'
        n, ligands = define_all_ligands(data)
        LIG_PROTEIN_PAIRS = [(n, ligands, receptor)]
    return LIG_PROTEIN_PAIRS


def create_ligands_file(ligands, output_filename=LIGANDS_TXT):
    """
    Creates a ligand list file for Uni-Dock's --ligand_index input.
    Returns:
        str: Path to the generated ligand index file.
    """
    with open(output_filename, 'w') as f:
        for ligand in ligands:
            f.write(f"{PATH_TO_LIGANDS}{ligand}\n")
    print(f"Ligand index file written to: {output_filename}")


def create_receptorname_file(receptor, output_filename=RECEPTOR_TXT):
    with open(output_filename, 'w') as f:
        f.write(f"{PATH_TO_PROTEINS}{receptor}")
    print(f"Receptor name written to file: {output_filename}")


def get_receptor_name(output_filename=RECEPTOR_TXT):
    with open(output_filename, 'r') as f:
        lines = f.readlines()
    return lines[0]


if __name__ == "__main__":
    batches = load_pairs()
    print(batches)

