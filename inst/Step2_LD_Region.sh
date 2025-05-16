#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --job-name=LD_Block
#SBATCH --mem=32000
#SBATCH --partition=main,mulan
#SBATCH --array=1-1368%400
#SBATCH --output=/net/fantasia/home/borang/MALGC/out/LD_block_%a.out
#SBATCH --error=/net/fantasia/home/borang/MALGC/out/LD_block_%a.err

bash

# Set variables (adjust these paths accordingly)
INPUT_DIR="/net/fantasia/home/borang/MALGC/pipeline_example"

ANCESTRY_1="EUR"
ANCESTRY_2="EAS"
TRAIT="LDL"
LD_BLOCK_FILE="/net/fantasia/home/borang/MALGC/ld_blocks/grch37.eur.eas.loci.bed"
R_SCRIPT="/net/fantasia/home/borang/MALGC/MALGC_software/Data_Process/Step2_LD_Region.R"


let k=0
for block in {1..1368}; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript ${R_SCRIPT} \
  --input_dir ${INPUT_DIR} \
  --ancestry_1 ${ANCESTRY_1} \
  --ancestry_2 ${ANCESTRY_2} \
  --trait ${TRAIT} \
  --ld_block_file ${LD_BLOCK_FILE} \
  --block_index ${SLURM_ARRAY_TASK_ID}
fi
done
