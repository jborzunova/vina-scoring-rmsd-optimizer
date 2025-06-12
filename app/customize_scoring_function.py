import os
import pandas as pd
import numpy as np
import itertools
import pickle
from scipy.optimize import minimize
import scipy
import time
import json
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from parameters import *
from load_data import *
from process_output import rmsd_from_output
from concurrent.futures import ProcessPoolExecutor, as_completed
from logger import *


def leaky_relu(f, alpha=LEAKAGE_PARAMETER):
    return f-2.5+0.25 if f > 2.5 else alpha * f


def unidock_custom_run(x):
    cmd_line = f'unidock --receptor {RECEPTOR} --ligand_index ../{LIGANDS_TXT} --dir . --search_mode fast --seed 2007 \
                  --center_x 15 --center_y 68 --center_z 39 \
                  --size_x 14 --size_y 22 --size_z 30 \
                  --weight_gauss1 {x[0]} \
                  --weight_gauss2 {x[1]} \
                  --weight_repulsion {x[2]} \
                  --weight_hydrophobic {x[3]} \
                  --weight_hydrogen {x[4]} > /dev/null 2>&1 \n' 
    log(os.getcwd())
    log(cmd_line)
    os.system(cmd_line) # in case you want it silent:  > /dev/null 2>&1    


def callback(x, f, context):  # Dual_annealing algorithm
    raw_mean_RMSD = HISTORY.get_raw_rmsd(x)
    log(f"Best Parameters so far: {x}, Raw mean RMSD value: {raw_mean_RMSD}, Objective Value: {f}, success = {context}")
    # Append the data to VALUES and save
    VALUES.append(x, raw_mean_RMSD, f)
    VALUES.save('loc_min_history.pkl')
    HISTORY.save('all_history.pkl')
    if f <= leaky_relu(TARGET_VALUE):
        log(f"Stopping optimization: reached target value.")
        return True
    else:
        return False


def callback_diffevol(intermediate_result):  # Differential Evolution algorithm. Callback is different for different algorithms
    f = intermediate_result.fun
    x = intermediate_result.x
    success = intermediate_result.success
    return callback(x, f, success)

####################################################################################################
def docking_functional(x):
    # x = [xi] - эти коэффициенты будут оптимизироваться
    # Это веса для энергетических слагаемых в скоринг-функции Vina
    # Цель: минимизировать общий RMSD для всех лигандов колхицинового сайта
    log(f'x is {x}')
    unidock_custom_run(x)  # it runs on GPU in parallel for all ligands and proteins

    # ---- Calculate RMSDs based on generated files ----
    all_rmsds = []
    with ProcessPoolExecutor() as executor:
        output_files = [file for file in os.listdir('.') if file.endswith('_out.pdbqt')]
        # Submit tasks for each file
        futures = {executor.submit(rmsd_from_output, output): output for output in output_files}
        # Gather results as they are completed
        for future in as_completed(futures):
            try:
                RMSD = future.result()
                all_rmsds.append(RMSD)
            except Exception as e:
                log(f'Error calculating RMSD: {e}')

    # ---- Gathering results ----
    raw_mean_rmsd = np.mean(all_rmsds)
    result = np.mean([leaky_relu(rmsd) for rmsd in all_rmsds])
    HISTORY.append(x, raw_mean_rmsd, result)
    log(f'mean_RMSD = {raw_mean_rmsd}')
    log(f'mean_leaky_relu_RMSD = {result}')
    return result
####################################################################################################
if __name__ == "__main__":

    # ---- Time ----
    start_time = time.time()
    
    # ---- Data Loading ----
    RECEPTOR = get_receptor_name()
    # ligands are already in the file 
    log(f'receptor name is {RECEPTOR}')

   # ---- CHDIR ----
    os.chdir(WORKDIR)
    log(f'I am in {os.getcwd()}')
 
    # ---- Processing ----
    # Calling minimize using the global optimization method Dual Annealing (with early stopping)
    HISTORY = History()  # stores the history of all objective function evaluations. Needed to output the raw value np.mean(RMSDs) instead of np.mean([leaky_relu(rmsd) for rmsd in RMSDs]).
    # Additionally, this will make it easier to resume calculations later.
    VALUES = History()  # stores the values of local minima
    if OPTIMIZ == 0:
        result = scipy.optimize.dual_annealing(docking_functional, BOUNDS, x0=X0, seed=2007, maxiter=MAXITER, no_local_search=NO_LOCAL_SEARCH, minimizer_kwargs=MINIMIZER_KWARGS, callback=callback)
    elif OPTIMIZ == 1:
        result = scipy.optimize.differential_evolution(docking_functional, BOUNDS, x0=X0, seed=2007, maxiter=MAXITER, tol=1e-3, disp=True, callback=callback_diffevol)
    else:
        log('Nothing to do. Assign OPTIMIZ to 0 or 1')

    log('')
    log("--- Work Time: %s seconds ---" % round((time.time() - start_time), 1))

