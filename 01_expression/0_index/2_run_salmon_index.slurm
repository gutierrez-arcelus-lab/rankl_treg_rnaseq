#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J Salmon-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/Salmon-index.%j.out
#SBATCH -e logs/Salmon-index.%j.err

source /programs/biogrids.shrc
export SALMON_X=1.10.1

cd "${SLURM_SUBMIT_DIR}"

GENTROME="./data/reference/gentrome.fa.gz"     
DECOYS="./data/reference/decoys.txt"
OUT="./salmon_vM39_pri_k21"

rm -rf "${TEMP_WORK}/temp/salmon" 

# Run indexing
salmon index \
    -t "$GENTROME" \
    -i "$OUT" \
    --decoys "$DECOYS" \
    -k 21 \
    -p "$SLURM_CPUS_PER_TASK" \
    --tmpdir "${TEMP_WORK}/temp/salmon"
