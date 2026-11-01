#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=24GB
#SBATCH -t 06:00:00
#SBATCH -p bch-compute
#SBATCH -J Salmon-quant
#SBATCH --array=1-38          
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/salmon-quant_%a.out
#SBATCH -e logs/salmon-quant_%a.err

source /programs/biogrids.shrc
export SALMON2_X=2.8.0

cd "${SLURM_SUBMIT_DIR}"

META="./data/metadata_salmon_merged.tsv"
INDEX="../../01_prepare_ref/salmon_vM39_pri_k21"
OUTBASE="./output/salmon"

TARGET_ROW=$((SLURM_ARRAY_TASK_ID + 1))
SAMPLE=$(awk -F'\t' -v r="$TARGET_ROW" 'NR==r{print $1}' "$META")
READ1=$(awk -F'\t' -v r="$TARGET_ROW" 'NR==r{print $2}' "$META")
READ2=$(awk -F'\t' -v r="$TARGET_ROW" 'NR==r{print $3}' "$META")

salmon quant \
    -i "$INDEX" \
    -l A \
    -1 ${READ1//,/ } \
    -2 ${READ2//,/ } \
    --gcBias --seqBias \
    --numBootstraps 100 \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "${OUTBASE}/${SAMPLE}"

