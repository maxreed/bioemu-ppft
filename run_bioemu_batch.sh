#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:h100:1       # Number of GPUs (per node)
#SBATCH --mem=16G               # memory (per node)
#SBATCH --time=2-00:00          # time (DD-HH:MM)
#SBATCH --output=run_bioemu_batch.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

python my_scripts/run_bioemu_batch.py --m3a_dir data/msa_kras_mutants --out_dir outputs

