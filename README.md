
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

1.  **Read QC and trimming** (`00_fastq_qc`): fastp 0.23.4 (adapter
    detection for paired-end reads, minimum read length 21), summarized
    with MultiQC. Libraries with \< 1 million read pairs after filtering
    are excluded.
2.  **Reference** (`01_prepare_ref`): GENCODE vM39 primary-assembly
    transcripts with the GRCm39 genome as decoys; Salmon 2.8.0
    decoy-aware index with k = 21.
3.  **Quantification and differential expression** (`02_analysis`), in
    two branches that differ in how technical replicates are handled:
    -   `merge_by_sample` (main analysis): the passing libraries of each
        biological sample are quantified together in one Salmon run
        (`--gcBias --seqBias --numBootstraps 100`). Import with tximport
        (estimated counts) and edgeR `DGEListFromTximport(divide = TRUE)`;
        edgeR quasi-likelihood tests.
    -   `merge_by_run`: each library is quantified separately (same
        settings and import); limma-voom (`voomLmFit()`) with the
        biological sample as a blocking factor.

    In both, genes are filtered with `filterByExpr()` using the
    experimental groups, normalized with TMM, and tested within each
    experiment. Sample QC uses PCA and correlation of principal
    components with technical covariates.
4.  **Figures** (`03_figures`): built from the `merge_by_sample`
    results. Gene set enrichment uses fgsea with MSigDB Hallmark gene
    sets for mouse (msigdbr).

## Repository structure

```
00_fastq_qc/        fastp, MultiQC
01_prepare_ref/     reference download and preparation, Salmon index
02_analysis/
  merge_by_sample/  0_make_samplesheet.R, 1_run_salmon.slurm, 2_process.R,
                    3_pca.R, 4_run_edger.R
  merge_by_run/     same steps, one Salmon run per library
03_figures/         figure scripts (merge_by_sample results)
```

Scripts are run from their own folder, in numbered order.

