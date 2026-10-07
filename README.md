
# RANKL in regulatory T cells: bulk RNA-seq analysis

This repository contains the code for the bulk RNA-seq analyses in:

> *Manuscript title* (in preparation). Authors.

## Data

Libraries were prepared with SMART-seq2 (low input) by the Broad
Institute and sequenced as 38-bp paired-end reads. Each biological
sample was split into two independently prepared libraries (technical
replicates).

| Experiment | Cells | Comparison | Design |
|----|----|----|----|
| 1 | Treg and Teff | Treg vs Teff | Paired by mouse, `~ mouse + sex + genotype` |
| 2 | iTregs | Foxp3<sup>WT</sup> vs Foxp3<sup>ΔRANKL</sup> | Unpaired, `~ sex + genotype` |
| 3 | iTregs | IL-1β vs IL-6 | Paired by mouse, `~ mouse + condition` |
| 4 | iTregs | RANKL<sup>+</sup> vs RANKL<sup>−</sup> | Paired by mouse, `~ mouse + condition` |
| 5 | Thymic Tregs | Foxp3<sup>WT</sup> vs Foxp3<sup>ΔRANKL</sup> | Unpaired, `~ sex + genotype` |

Raw data: *GEO accession to be added.*

## Pipeline

1.  **Read QC and trimming:** fastp 0.23.4 (adapter detection for
    paired-end reads, minimum read length 21).
2.  **Reference:** GENCODE vM39 primary-assembly transcripts with the
    GRCm39 genome as decoys (`build_gentrome.R`).
3.  **Quantification:** Salmon 1.10.1, decoy-aware index with k = 21,
    `--gcBias --seqBias`, one run per library.
4.  **Import and QC:** tximport
    (`countsFromAbundance = "lengthScaledTPM"`). Libraries with \< 1
    million fragments are removed, and the two technical libraries of
    each sample are summed. Sample QC uses PCA and correlation of
    principal components with technical covariates.
5.  **Differential expression:** edgeR quasi-likelihood pipeline
    (primary), with DESeq2 as a replication. Genes are filtered with
    `filterByExpr()` using the experimental groups, normalized with TMM,
    and tested within each experiment.
6.  **Gene set enrichment:** fgsea with MSigDB Hallmark gene sets for
    mouse (msigdbr), ranked by the signed quasi-likelihood statistic.
