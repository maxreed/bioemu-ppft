from bioemu.sample import main as sample
import glob, os, argparse

def main():
    p = argparse.ArgumentParser(
        description="Run bioemu on every .a3m file in a given directory."
    )
    p.add_argument('--m3a_dir',   required=True, help="Input directory containing MSA files.")
    p.add_argument('--out_dir',   required=True, help="Output directory for bioemu results.")
    p.add_argument('--num_states',   required=False, default=1000, help="Number of attempted states to generate for each protein sequence.")
    p.add_argument('--ckpt_path',   required=False, default=None, help="Checkpoint, i.e. the model version in use.")
    p.add_argument('--model_config_path',   required=False, default=None, help="Config file for the current model.")
    args = p.parse_args()

    num_states = int(args.num_states)
    msa_files = glob.glob(os.path.join(args.m3a_dir, "*" + ".a3m"))

    for msa_file in msa_files:
        print("Running BioEmu for " + os.path.basename(msa_file)[:-4] + "...") # -4 in the msa_file index removes the .a3m ending
        sample(sequence=msa_file, num_samples=num_states, output_dir= args.out_dir + "/KRAS_" + os.path.basename(msa_file)[:-4], ckpt_path=args.ckpt_path, model_config_path=args.model_config_path)

if __name__ == "__main__":
    main()

