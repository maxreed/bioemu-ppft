#!/bin/bash
#SBATCH --account=def-awoolley
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G                # memory (per node)
#SBATCH --time=0-01:00          # time (DD-HH:MM)
#SBATCH --output=assess_state_occupancy_v21.log

module load StdEnv/2023 python/3.11 cuda/12.2 gcc/12.3
source bioemu_env/bin/activate

STATE1_CIF=/home/pmmreed/links/projects/def-awoolley/pmmreed/bioemu_project/potential_references/6mbq.cif
STATE2_CIF=/home/pmmreed/links/projects/def-awoolley/pmmreed/bioemu_project/potential_references/6xi7.cif
ENSEMBLE_DIR=outputs/v21_allVariants

MUTANTS=("WT" "Y4F" "L6M" "V7I" "V8I" "V9I" "Q22N" "L23M" "Y40F" "V44A" "C51A" "L52M" "L52V" "L53M" "L53V" "I55V" "T58S" "Q61R" "Y64A" "R68L" "Y71F" "M72F" "M72L" "M72V" "R73H" "R73S" "E76A" "E76Q" "F78Y" "L79M" "L79V" "C80I" "F82Y" "I84M" "I84V" "S89T" "F90W" "F90Y" "D92N" "I93M" "I93V" "H95K" "Y96F" "Y96L" "Y96W" "R97H" "E98D" "E98Q" "I100V" "K101R" "V103A" "V103L" "K104A" "K104N" "K104T" "V109I" "V109M" "P110A" "M111F" "V112I" "L113I" "V114I" "V125I" "A130V" "L133M" "L133V" "I139V" "F141Y" "I142F" "R149H" "R149K" "A155C" "F156Y" "T158L" "L159V" "V160A" "V160I" "E162L" "I163M" "I163V" "H166Q" "K167E")

# Print CSV header once
echo "name,n_frames,binary_occupancy,soft_occupancy"

for s in "${MUTANTS[@]}"; do
    python my_scripts/assess_state_occupancy.py \
        --xtc  ${ENSEMBLE_DIR}/KRAS_${s}/samples.xtc \
        --pdb  ${ENSEMBLE_DIR}/KRAS_${s}/topology.pdb \
        --name KRAS_${s} \
        --state1_cif ${STATE1_CIF} \
        --state2_cif ${STATE2_CIF}
done
