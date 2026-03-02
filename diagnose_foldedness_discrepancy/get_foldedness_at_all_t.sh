#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:h100:1       # Number of GPUs (per node)
#SBATCH --mem=16G           # memory (per node)
#SBATCH --time=0-02:00      # time (DD-HH:MM)
#SBATCH --output=get_foldedness_at_all_t.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

python diagnose_foldedness_trajectory.py \
    --msa_file /home/pmmreed/links/projects/def-awoolley/pmmreed/bioemu_project/data/msa_kras_mutants/WT.a3m \
    --ckpt_path /home/pmmreed/links/scratch/2026-03-01_bioemu_refinement_KRAS_CEST_attempt6/bioemu_finetune_output_v6/checkpoints/bioemu_cest_finetuned_v6_epoch199_properFormat.ckpt \
    --model_config_path /home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/checkpoints/bioemu-v1.1/config.yaml \
    --cache_so3 /home/pmmreed/sampling_so3_cache \
    --cache_embeds_dir /home/pmmreed/.bioemu_embeds_cache \
    --output_dir foldedness_trajectory_diagnostic \
    --num_samples 30 \
    --N_full 35 \
    --N_probe 1 \
    --max_t 0.99 \
    --eps_t 0.001

