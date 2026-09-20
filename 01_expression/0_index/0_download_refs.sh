#!/usr/bin/bash

REF=./data/reference
mkdir -p "$REF" logs
cd "$REF"
 
BASE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39
GENOME_GZ=GRCm39.primary_assembly.genome.fa.gz          # PRI genome (decoy)
TX_ALL_GZ=gencode.vM39.transcripts.fa.gz                # ALL-level transcript fasta
GTF_PRI_GZ=gencode.vM39.primary_assembly.annotation.gtf.gz   # PRI GTF (for id list)
 
# Downloads (skip if present)
[[ -f $GENOME_GZ  ]] || wget -q ${BASE}/${GENOME_GZ}
[[ -f $TX_ALL_GZ  ]] || wget -q ${BASE}/${TX_ALL_GZ}
[[ -f $GTF_PRI_GZ ]] || wget -q ${BASE}/${GTF_PRI_GZ}
