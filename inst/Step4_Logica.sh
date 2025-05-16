#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --job-name=LDL
#SBATCH --mem=32000
#SBATCH --partition=main,mulan
#SBATCH --array=1-1368%400
#SBATCH --output=/net/fantasia/home/borang/MALGC/out/LD_block_%a.out
#SBATCH --error=/net/fantasia/home/borang/MALGC/out/LD_block_%a.err

bash

let k=0
for block in {1..1368}; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
	Rscript /net/fantasia/home/borang/MALGC/MALGC_software/Data_Process/Step4_Logica.R ${block}
fi
done
