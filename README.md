# WGCNA_metaNET
This repository will cover the analytical workflow use to conduct co-expression analysis of RNAseq data. Weighted geene co-expression network analysis (https://doi.org/10.1186/1471-2105-9-559) was used as a framework for this analysis. All analyses were performed on the CAMH SCC cluster given the exceptionally high computational requirements of permuted co-expression analyses. Additionally, the default linear algebra backend in R is too slow to perform this analysis in a reasonable time - thus R installed with openBLAS was used to speed-up matrix multiplication steps. 

These analyses were conducted on the Sibille Lab "PITT Tetrad" dataset. Laser-dissected populations of L2/3 Pyramidal, L5/6 Pyramidal, SST-neurons, PV-neurons, and VIP-neurons were collected (130 cells/sample) and sequenced on the NovaSeq platform. Cells were dissected from the subgenual anterior cingulate cortex, from a post-mortem cohort containing 19 matched tetrads of major depressive disorder (MDD), bipolar disorder (BPD), schizophrenia (SCZ), and control subjects. Here, analyses focused on (1) changes in co-expression patterns within each cell-type across disorder contrasts (i.e. psych groups v.s. controls) and (2) changes in transcriptome-wide co-expression between cell-types via eigengene networks.

1. Data import organization, QC, and parameter identification.
   * `QC and input generation.R` and `Differential WGCNA paramater testing.R`
   * Following the WGCNA guidelines, data is quality controlled for noisy/0-count genes and outlier samples.
   * A range of WGCNA parameters are explored and dendrograms/# of genes clustering in modules are analyzed.

2. Generate observed modules.
   * `Observed modules.R`
   * Generates WGCNA modules and associated output files (module eigengenes, dendrograms, module colour vectors, etc).

3. Pre-generate permutation analysis blocks.
   * `Permutation list generation.R`
   * This uses a set.seed(12345) in R to generate randomized, but standardized, permutations of subject IDs for permutation testing.
   * Pre-generation is necessary since co-expression analysis will be parallelized as 200-permutation chunks.

4. Perform permutation analysis.
   * Performed by `Density thresholded permutation test.R`
   * As an added level of confidence, in each permutation the differences between psych and control groups are compared at 20 eigengene correlation p-value thresholds.
   * Differences which are significantly different across 5 or more consecutive thresholds will be considered significant for the cell-cell connection.
   * Parallelized using `Differential_WGCNA_Permutations_MDD_parallelization.sh`, etc.

5. Combine permutation analysis and generate final outputs.
   * `Merge permutation output.R`- Combine all 200-chunk outputs into single file.
   * `Density threshold final results.R` - Generates empirical p-values at each threshold for all cell-types and provide detailed output plots visualizing permutated and observed results.
   * `Within cell-type coexpression results.R` - Generates facetted histograms showing module size distribution in each cell-type and experimental condition. Additionally, generates module-trait correlations and exports plots of module membership v.s. gene significance for age-correlated modules.
