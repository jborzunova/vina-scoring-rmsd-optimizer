# conda activate unidock
import os
from parameters import *
import inspect
from app.load_data import * 


# SLURM script template
slurm_script = f"""#!/bin/bash
#SBATCH -t {TIME} # timelimit here 3 days
#SBATCH -p {PARTITION}        # partition (see sinfo command)
#SBATCH --gres=gpu:1
#SBATCH -J {JOBNAME}     # name
#SBATCH -o {STDOUT}    # stdout filename %j will be replaced with job id
#SBATCH -e {STDERR}    # stderr filename
#SBATCH -N {NNODES}
#SBATCH --ntasks={NTASKS}   # in total across all requested nodes
#SBATCH --cpus-per-task={CPUS_PER_TASK}

# Activate conda environment
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate unidock_env

# Run the Python script with MPI
python3 ../app/customize_scoring_function.py > {LOG_FILE}
"""
batches = load_pairs()
# Distribute by index parity
queue1 = [batch for batch in batches if batch[0] % 2 == 0]
queue2 = [batch for batch in batches if batch[0] % 2 == 1]
print('queue1:', queue1)
print('queue2:', queue2)

# Function to submit a chain of jobs
def submit_chain(queue, queue_id):
    prev_jobid = None
    for batch in queue:
        i, ligands, receptor = batch
        job_index = f"{i}"
        
        # Create folders and files
        work_path = WORK_FOLDER + f'_{job_index}'
        if not os.path.exists(work_path):
            os.makedirs(work_path)
        os.chdir(work_path)
        if not os.path.exists(WORKDIR):
            os.makedirs(WORKDIR)

        create_ligands_file(ligands)
        create_receptorname_file(receptor)
        os.system('cp ../parameters.py .')

        # Write SLURM script
        with open("slurm_job.sh", "w") as f:
            f.write(slurm_script)
        
        # Submit with or without dependency
        if prev_jobid is None:
            result = subprocess.run(["sbatch", "slurm_job.sh"], stdout=subprocess.PIPE, text=True)
        else:
            result = subprocess.run(["sbatch", f"--dependency=afterany:{prev_jobid}", "slurm_job.sh"], stdout=subprocess.PIPE, text=True)

        # Extract job ID
        prev_jobid = result.stdout.strip().split()[-1]
        print(f"Submitted job {prev_jobid} for batch {job_index}")
        
        os.chdir(ROOT_FOLDER)

# Submit both chains
submit_chain(queue1, "queue1")
submit_chain(queue2, "queue2")

