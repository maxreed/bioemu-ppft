#!/usr/bin/env python
"""
Test BioEmu inference on GPU. (You need to request a GPU node before running this).
This script exists because Rorqual compute nodes can't access the internet.
This script assumes you have already generated single and pair embeddings for your sequence (you can use get_colabfold_embeds if you have BioEmu
installed), then copied those into $HOME/.bioemu_embeds_cache.
If you have done that, BioEmu will look in that directory for the query sequence embeddings (it will ALWAYS do this by default). Naming is done via
hashlib, so it'll find it so long as naming is left unchanged from what get_colabfold_embeds uses.
Because this avoids having to query the MSA server (which requires internet) this means that it allows BioEmu to run on a compute node.
Now, you can also just feed in an MSA (made by colabfold_batch or whatever other program you like). That ALSO allows BioEmu to run on a compute
node, and it's probably the better thing to do. That having been said, because embeddings are only generated once, I really don't think it matters
what approach you use.
Look in test_job.sh to see the example of how to do this with an input MSA instead of pre-computed embeddings.
"""
import os
from bioemu.sample import main as sample

# Set paths
checkpoint = 'checkpoints/bioemu-v1.1/checkpoint.ckpt'
config = 'checkpoints/bioemu-v1.1/config.yaml'
output_dir = 'outputs/bioemu_test_PYP_M100E'

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

print("Testing BioEmu inference...")
print(f"Checkpoint: {checkpoint}")
print(f"Config: {config}")
print(f"Output: {output_dir}")

# Generate structures
sample(
    sequence='EHVAFGSEDIENTLAKMDDGQLDGLAFGAIQLDGDGNILQYNAAEGDITGRDPKQVIGKNFFKDVAPCTDSPEFYGKFKEGVASGNLNTMFEYTFDYQETPTKVKVHMKKALSGDSYWVFVKRV',      # This is PYP M100E, because why not?
    num_samples=25,
    output_dir=output_dir,
    ckpt_path=checkpoint,
    model_config_path=config,
    # use_msa_server=False,     # Add this if needed (update: you don't need to)
)

print(f"\nDone! Check {output_dir} for PDB files")
