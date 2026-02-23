#!/usr/bin/env python
"""
Test BioEmu inference on GPU
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
    sequence='EHVAFGSEDIENTLAKMDDGQLDGLAFGAIQLDGDGNILQYNAAEGDITGRDPKQVIGKNFFKDVAPCTDSPEFYGKFKEGVASGNLNTMFEYTFDYQETPTKVKVHMKKALSGDSYWVFVKRV',      # Chignolin (10 residues)
    num_samples=25,
    output_dir=output_dir,
    ckpt_path=checkpoint,
    model_config_path=config,
    # use_msa_server=False,     # Add this if needed
)

print(f"\nDone! Check {output_dir} for PDB files")
