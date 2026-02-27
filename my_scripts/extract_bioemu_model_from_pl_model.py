import torch
from pathlib import Path

lightning_ckpt_path = '/home/pmmreed/links/scratch/2026-02-26_bioemu_refinement_KRAS_CEST_attempt4/bioemu_finetune_output_v4/checkpoints/bioemu-foldedness-epoch=10-val/kras_CEST_refined_v4_epoch10.ckpt'
output_path = '/home/pmmreed/links/scratch/2026-02-26_bioemu_refinement_KRAS_CEST_attempt4/bioemu_finetune_output_v4/checkpoints/kras_CEST_refined_v4_epoch10_bioemu.ckpt'

# Load the Lightning checkpoint
lightning_ckpt = torch.load(lightning_ckpt_path, map_location="cpu", weights_only=True)

# Extract the nested state dict and strip the "score_model." prefix
lightning_state_dict = lightning_ckpt["state_dict"]
prefix = "score_model."
bioemu_state_dict = {
    k[len(prefix):]: v 
    for k, v in lightning_state_dict.items() 
    if k.startswith(prefix)
}

# Save in the flat format bioemu's load_model expects
torch.save(bioemu_state_dict, output_path)
print(f"Saved {len(bioemu_state_dict)} parameter tensors to {output_path}")
