#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:h100:1       # Number of GPUs (per node)
#SBATCH --mem=16G               # memory (per node)
#SBATCH --time=0-12:00          # time (DD-HH:MM)
#SBATCH --output=refine_bioemu_CEST_v16.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source /home/pmmreed/links/projects/def-awoolley/pmmreed/bioemu_project/bioemu_env/bin/activate

# wandb might try to write to the main server - which is impossible on a compute mode.
# not sure if this is an actual problem, but i did this pre-emptively.
export WANDB_MODE=offline

# this should help reduce GPU memory requirements
export PYTORCH_ALLOC_CONF=expandable_segments:True

python finetune_foldedness_v4b_hybrid_loss.py \
    /home/pmmreed/.bioemu_embeds_cache \
    /home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/data/kras_mutants_data.csv \
    /home/pmmreed/links/scratch/2026-03-27_bioemu_refinement_KRAS_CEST_attempt16/bioemu_finetune_output_v16 \
    /home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/data/kras_mutant_sequences.csv \
    --ckpt_path /home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/checkpoints/bioemu-v1.1/checkpoint.ckpt \
    --model_config_path /home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/checkpoints/bioemu-v1.1/config.yaml \
    --cache_so3 /home/pmmreed/sampling_so3_cache \
    --structures_dir /home/pmmreed/links/projects/def-awoolley/pmmreed/bioemu_project/outputs/bioemu1p1_KRAS \
    --batch_size 2 \
    --gpus 1 \
    --seed 42 \
    --n_replications 8 \
    --max_epochs 100 \
    --trainable_components model_nn.x1d_proj model_nn.x2d_proj model_nn.rp_proj model_nn.st_module.encoder.layers model_nn.st_module.diff_head \
    --mid_t 0.5 \
    --N_rollout 17 \
    --record_grad_steps 15 16 17 \
    --precision bf16-mixed \
    --dsm_weight 0.2


