library(tidyverse)
library(Biostrings)

ref <- "./data/reference"
gtf <- file.path(ref, "gencode.vM39.primary_assembly.annotation.gtf.gz")

gencode_gene <-
    gtf |>
    rtracklayer::import(feature.type = "gene") |>
    as_tibble()

gencode_tx <-
    gtf |>
    rtracklayer::import(feature.type = "transcript") |>
    as_tibble()

# ALL-level transcript fasta; strip pipe header to the versioned ENSMUST id
fasta <- readDNAStringSet(file.path(ref, "gencode.vM39.transcripts.fa.gz"), "fasta")
names(fasta) <- str_extract(names(fasta), "^([^|]+)", group = 1)

# Safety check: every PRI transcript must be present in the ALL fasta
all(gencode_tx$transcript_id %in% names(fasta))

# Subset ALL -> PRI
# This shows that the original transcript fasta is actually already 'PRI'
fasta_pri <- fasta[gencode_tx$transcript_id]

writeXStringSet(fasta_pri,
                file.path(ref, "gencode.vM39.primary_assembly.transcripts.fa.gz"),
                compress = TRUE)

# Gentrome = PRI transcripts + primary-assembly genome (transcripts first)
fasta_genome <- readDNAStringSet(file.path(ref, "GRCm39.primary_assembly.genome.fa.gz"))
names(fasta_genome) <- str_extract(names(fasta_genome), "\\S+")

gentrome <- c(fasta_pri, fasta_genome)
writeXStringSet(gentrome, file.path(ref, "gentrome.fa.gz"), compress = TRUE)

# Decoys = genome sequence names
write_lines(names(fasta_genome), file.path(ref, "decoys.txt"))

# tx2gene for tximport 
gencode_tx |>
    select(transcript_id, gene_id) |>
    write_tsv(file.path(ref, "tx2gene.tsv"))

# Gene metadata
gencode_gene |>
    select(gene_id, gene_name, gene_type) |>
    write_tsv(file.path(ref, "gencode.vM39.gene_metadata.tsv"))

