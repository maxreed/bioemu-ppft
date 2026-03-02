#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:h100:1       # Number of GPUs (per node)
#SBATCH --mem=16G               # memory (per node)
#SBATCH --time=0-10:00          # time (DD-HH:MM)
#SBATCH --output=run_bioemu_batch_v6_400_states.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

python my_scripts/run_bioemu_batch.py --m3a_dir data/msa_kras_mutants --out_dir outputs/v6_allVariants --num_states 400 \
	--ckpt_path /home/pmmreed/links/scratch/2026-03-01_bioemu_refinement_KRAS_CEST_attempt6/bioemu_finetune_output_v6/checkpoints/bioemu_cest_finetuned_v6_epoch199_properFormat.ckpt \
	--model_config_path /home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/checkpoints/bioemu-v1.1/config.yaml

