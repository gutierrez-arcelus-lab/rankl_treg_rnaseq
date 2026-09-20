#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=32gb
#SBATCH -t 06:00:00
#SBATCH -p bch-compute
#SBATCH -J Salmon-quant
#SBATCH --array=1-96           
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/Salmon-quant_%a.out
#SBATCH -e logs/Salmon-quant_%a.err

source /programs/biogrids.shrc
export SALMON_X=1.10.1

cd "${SLURM_SUBMIT_DIR}"

META="../../00_fastq_qc/data/metadata_fastp.tsv"
INDEX="../0_index/salmon_vM39_pri_k21"
TRIM="../../00_fastq_qc/output/fastq_trimmed" 
OUTBASE="./output/salmon"

TARGET_ROW=$((SLURM_ARRAY_TASK_ID + 1))
WELL=$(awk -F'\t' -v r="$TARGET_ROW" 'NR==r{print $1}' "$META")
SAMPLE=$(awk -F'\t' -v r="$TARGET_ROW" 'NR==r{print $2}' "$META")
RUN="${WELL}_${SAMPLE}"

salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "${TRIM}/${RUN}_R1.fastq.gz" \
    -2 "${TRIM}/${RUN}_R2.fastq.gz" \
    --numBootstraps 100 \
    --gcBias \
    --seqBias \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "${OUTBASE}/${RUN}"
