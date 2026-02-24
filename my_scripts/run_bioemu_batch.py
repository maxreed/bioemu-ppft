from bioemu.sample import main as sample
import glob, os, argparse

def main():
    p = argparse.ArgumentParser(
        description="Run bioemu on every .a3m file in a given directory."
    )
    p.add_argument('--m3a_dir',   required=True, help="Input directory containing MSA files.")
    p.add_argument('--out_dir',   required=True, help="Output directory for bioemu results.")
    p.add_argument('--num_states',   required=False, default=1000, help="Number of attempted states to generate for each protein sequence.")
    args = p.parse_args()

    num_states = int(args.num_states)
    msa_files = glob.glob(os.path.join(args.m3a_dir, "*" + ".a3m"))
    print(msa_files)

    #sample(sequence='data/msa_kras_mutants/Y40F.a3m', num_samples=1000, output_dir='test_RAS_Y40F')

if __name__ == "__main__":
    main()

