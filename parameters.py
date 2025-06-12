import subprocess
import numpy as np 
import pandas as pd

JOBNAME = 'uni2'
PARTITION = 'gpu'
OPTIMIZ = 1  # 0 for dual_annealing and 1 for differential evolution algorithm 
NO_LOCAL_SEARCH = False  # this parameter is relevant for dual_annealing only
TIME = '0-00:01:00'
MODE = 4  # 0 for customizing of all center receptors of dataset with reference docking; 1 for one receptor; 2 for two receptors; 3 for 5XKH and train_{FOLD_NUMBER}
FOLD_NUMBER = 0
X0 = np.array([-0.035579, -0.005156, 0.840245, -0.035069, -0.587439])  # default Vina parameters
NNODES = 1
NTASKS = 1 # total for all nodes
CPUS_PER_TASK = 1
ROOT_FOLDER = '/home/jborzunova/colchicine_site/customization/'
DF = '/home/jborzunova/colchicine_site/data.csv'
WORK_FOLDER = ROOT_FOLDER + JOBNAME
BOUNDS = [(-1,0), (-1,0), (0,1), (-1,0), (-1,0)]
MINIMIZER_KWARGS = {
    "method": "BFGS",
    "bounds": BOUNDS,
    "options": {"gtol": 1e-3}
}  # this parameter is relevant for dual_annealing only
MAXITER = int(2e9)
TARGET_VALUE = 1.0
LIGANDS_TXT = 'ligands.txt'
RECEPTOR_TXT = 'receptor_name.txt'
STDOUT = 'stdout'
STDERR = 'stderr'
LOG_FILE = 'log'
DOCK_FILE = 'dock.score'  # file name
PATH_TO_PROTEINS = '/home/jborzunova/colchicine_site/protein_prep/'
PATH_TO_LIGANDS = '/home/jborzunova/colchicine_site/ligand_prep_3Dgen/'
PATH_TO_NATIVE_LIGANDS = '/home/jborzunova/colchicine_site/ligand_prep_3Dnative/'
PATH_TO_BOX = '/home/jborzunova/colchicine_site/'
LEAKAGE_PARAMETER = 0.1
WORKDIR = 'workdir'
