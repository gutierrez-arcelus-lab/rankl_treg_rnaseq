#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=24GB
#SBATCH -t 06:00:00
#SBATCH -p bch-compute
#SBATCH -J Salmon-quant
#SBATCH --array=1-65
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/salmon-quant_%a.out
#SBATCH -e logs/salmon-quant_%a.err

source /programs/biogrids.shrc
export SALMON2_X=2.8.0

cd "${SLURM_SUBMIT_DIR}"

META="./data/metadata_salmon_runs.tsv"
INDEX="../../01_prepare_ref/salmon_vM39_pri_k21"
TRIM="../../00_fastq_qc/output/fastq_trimmed"
OUTBASE="./output/salmon"

TARGET_ROW=$((SLURM_ARRAY_TASK_ID + 1))
RUN=$(awk -F'\t' -v r="$TARGET_ROW" 'NR==r{print $1}' "$META")

salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "${TRIM}/${RUN}_R1.fastq.gz" \
    -2 "${TRIM}/${RUN}_R2.fastq.gz" \
    --gcBias --seqBias \
    --numBootstraps 100 \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "${OUTBASE}/${RUN}"
