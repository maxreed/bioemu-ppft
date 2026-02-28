#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1       # we actually don't need a GPU for this one
#SBATCH --mem=16G               # memory (per node)
#SBATCH --time=0-01:00          # time (DD-HH:MM)
#SBATCH --output=test_job_newModel.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

MUTANTS=("Y71F" "Y4F" "V160I" "V103L" "R149K" "M72V" "M72F" "L79M" "L6M" "L53V" "L159V" "L113I" "K104A" "F90Y" "A155C" "A130V")

for s in "${MUTANTS[@]}"; do
  python my_scripts/new_loss.py --xtc outputs/bioemu1p1_KRAS/KRAS_$s/samples.xtc --pdb outputs/bioemu1p1_KRAS/KRAS_$s/topology.pdb --name KRAS_$s
done
