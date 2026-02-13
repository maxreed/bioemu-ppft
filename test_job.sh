#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --gres=gpu:1              # Number of GPUs (per node)
#SBATCH --mem=500M                # memory (per node)
#SBATCH --time=0-01:00            # time (DD-HH:MM)
#SBATCH --output=test_job.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

python my_scripts/test_bioemu_with_kras_wt.py
