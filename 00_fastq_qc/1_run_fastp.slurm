#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=4gb
#SBATCH --time=00:30:00
#SBATCH -p bch-compute
#SBATCH --array=1-96
#SBATCH --job-name=fastp
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/fastp_%a.out
#SBATCH -e logs/fastp_%a.err

source /programs/biogrids.shrc
export FASTP_X=0.23.4

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

META="./data/metadata_fastp.tsv"

TARGET_ROW=$((SLURM_ARRAY_TASK_ID + 1))

WELL=$(awk -F'\t' -v row="$TARGET_ROW" 'NR==row{print $1}' "$META")
SAMPLE=$(awk -F'\t' -v row="$TARGET_ROW" 'NR==row{print $2}' "$META")
FQ1=$(awk -F'\t' -v row="$TARGET_ROW" 'NR==row{print $4}' "$META")
FQ2=$(awk -F'\t' -v row="$TARGET_ROW" 'NR==row{print $5}' "$META")

OUTDIR="./output/fastp"              # QC reports (json/html)
TRIMDIR="./output/fastq_trimmed"     # trimmed reads for Salmon

# fastp takes ONE file per read, so concatenate the 4 lanes first.
FQ1_TMP="${TEMP_WORK}/fastq/${WELL}_${SAMPLE}_R1.fastq.gz"
FQ2_TMP="${TEMP_WORK}/fastq/${WELL}_${SAMPLE}_R2.fastq.gz"

cat ${FQ1//,/ } > "$FQ1_TMP"
cat ${FQ2//,/ } > "$FQ2_TMP"

fastp \
    -i "$FQ1_TMP" -I "$FQ2_TMP" \
    -o "${TRIMDIR}/${WELL}_${SAMPLE}_R1.fastq.gz" \
    -O "${TRIMDIR}/${WELL}_${SAMPLE}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --length_required 21 \
    --thread "$SLURM_CPUS_PER_TASK" \
    --json "${OUTDIR}/${WELL}_${SAMPLE}.fastp.json" \
    --html "${OUTDIR}/${WELL}_${SAMPLE}.fastp.html"
