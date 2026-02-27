#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:h100:1       # Number of GPUs (per node)
#SBATCH --mem=16G               # memory (per node)
#SBATCH --time=0-01:00          # time (DD-HH:MM)
#SBATCH --output=test_job_newModel.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

python my_scripts/test_bioemu_with_kras_wt.py
python my_scripts/new_loss.py --xtc test_RAS_WT_finetuned_v5/samples.xtc --pdb test_RAS_WT_finetuned_v5/topology.pdb --name WT_finetuned_v5
python my_scripts/new_loss.py --xtc test_RAS_WT_finetuned_v4/samples.xtc --pdb test_RAS_WT_finetuned_v4/topology.pdb --name WT_finetuned_v4
python my_scripts/new_loss.py --xtc test_RAS_wt/samples.xtc --pdb test_RAS_wt/topology.pdb --name WT_bioemu1p1
