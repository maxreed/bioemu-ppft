from bioemu.sample import main as sample

ckpt_path = '/home/pmmreed/links/scratch/2026-02-26_bioemu_refinement_KRAS_CEST_attempt4/bioemu_finetune_output_v4/checkpoints/bioemu-foldedness-epoch=10-val/kras_CEST_refined_v4_epoch10.ckpt'

model_config_path = '/home/pmmreed/links/scratch/2026-02-25_bioemu_refinement_KRAS_CEST_attempt1/checkpoints/bioemu-v1.1/config.yaml'

sample(sequence='data/msa_kras_mutants/WT.a3m', num_samples=1000, output_dir='test_RAS_WT_finetuned_v4', ckpt_path=ckpt_path, model_config_path=model_config_path)

