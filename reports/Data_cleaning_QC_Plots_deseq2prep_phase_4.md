
# Phase 4 — Data Cleaning & DESeq2 Preparation

**Project Title:** Transcriptomics-Based Drug Target Discovery in Diabetic Cardiomyopathy (T2D iPSC‑CM)

**Purpose of this note:** reproducible, human-readable record of everything performed in Phase 4: data validation, de-duplication, group-aware filtering, DESeq2 object creation, normalization/VST, and QC plots. Save this file in `reports/` and commit to your repo for methods/provenance.

**Date:** 2025-09-14  
---

## Short executive summary

* **Inputs:** aligned raw counts (`data/processed/expression_matrix_aligned.csv`) and metadata (`data/processed/design_matrix.csv`).
* **Main actions:** load + sanity-check, collapse duplicate genes (if any), drop all-zero genes, apply group-aware low-count filter (`count >= 10` in `>= 3` replicates within a group), build a `DESeqDataSet` with `design = ~ Treatment * Oxygen`, compute size factors and VST, run QC plots (library sizes, PCA, sample distance heatmap, VST boxplots), and export analysis-ready artifacts.
* **Numbers (your run):** start = **62,266** genes → after removing all-zero = **33,258** → after group-aware filtering = **18,090** genes.
* **Outputs (saved):** `data/processed/counts_synced.csv`, `data/processed/design_matrix_synced.csv`, `data/processed/counts_filtered.csv`, `data/processed/dds_phase4_raw.rds`, `data/processed/vsd_phase4.rds`, `data/processed/vst_matrix.csv`, and QC figures in `results/qc/`.

---

## Files & directories (project layout)

* `data/processed/expression_matrix_aligned.csv` — aligned raw counts (genes × samples; first column `gene`).
* `data/processed/design_matrix.csv` — sample metadata (Sample, Treatment, Oxygen, Replicate, Group).
* `data/processed/counts_synced.csv` — counts after alignment and reordering to metadata.
* `data/processed/design_matrix_synced.csv` — metadata cleaned & saved.
* `data/processed/counts_filtered.csv` — group-aware filtered counts (analysis-ready).
* `data/processed/dds_phase4_raw.rds` — saved DESeqDataSet object (pre-DESeq fit).
* `data/processed/vsd_phase4.rds` — VST-transformed object for QC & visualization.
* `data/processed/vst_matrix.csv` — raw numeric VST matrix (genes × samples) for sharing.
* `results/qc/` — PNGs: `library_sizes.png`, `pca_group.png`, `sample_distance_heatmap.png`, `vst_boxplot.png`.
* `reports/Phase4_DataCleaning_DESeq2_Prep.md` — this file.

---

## Environment & packages

**R** (the kernel used in Jupyter/VSCode). Important packages used:

* `readr`, `dplyr`, `tibble`, `stringr` — data IO and manipulation.
* `DESeq2` — DESeqDataSet, normalization, `DESeq()` pipeline.
* `vsn` — provides methods used by `vst()` for stable transforms (Bioconductor).
* `ggplot2` — plotting.
* `pheatmap` — heatmaps.
* `reshape2` (only used in optional ggplot boxplot transformation).

**Important:** install Bioconductor packages (DESeq2, vsn) *in the R installation your Jupyter kernel uses* (not necessarily in Conda unless the R kernel points to Conda's R). Example install snippet:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","vsn"))
```

---

# Step-by-step record, explanation, and reasoning

Each numbered step below matches the code ran. For each: what we did, why we did it, expected types/shape of inputs & outputs, typical checks, and biological significance for T2D cardiomyopathy drug discovery.

---

## STEP 1 — Load, align, and validate counts & metadata

**Code overview (conceptual):** read metadata, coerce factors and reference levels, read counts CSV, ensure `gene` column exists and move it to rownames, check that sample column names match metadata (set equality), reorder counts columns to metadata order.

**Why:** DESeq2 expects count matrix columns to correspond to sample rows in `colData`. Mistatches or misordered samples will produce incorrect modeling.

**Input:** `expression_matrix_aligned.csv` (format: `gene,GSMxxxx,...`), `design_matrix.csv` with `Sample` matching the GSM columns in counts.

**Checks performed:**

* `setequal(colnames(counts_mat), meta$Sample)` — same sample names present.
* `identical(colnames(counts_mat), meta$Sample)` — same order (reorder if needed).
* `all(counts_mat >= 0)`, no NAs, integer-like values.
* `!anyDuplicated(rownames(counts_mat))` — gene IDs unique (if not, handle below).

**Outputs:** `counts_synced.csv`, `design_matrix_synced.csv` saved to `data/processed/` for provenance.

**Biological significance:** correct sample matching is foundational — any downstream DE gene list depends on sample labels being accurate. For drug discovery, mislabeling could produce false targets.

**Quick notes:** include a small `print(head(...))` to confirm visually and `getwd()` to ensure relative file paths are correct.

---

## STEP 2 — De-duplicate genes, drop all-zero genes, group-aware filtering

### 2.1 De-duplicate (if present)

**Code behavior:** detect duplicated gene names (`duplicated(rownames(counts_mat))`). If length > 0, collapse duplicates by summing counts across duplicate rows with `group_by(gene) %>% summarise(across(..., sum))`.

**Why:** duplicates can stem from transcript-level vs gene-level aggregation, or multiple annotations mapping to the same symbol. Summing is conservative and preserves total gene counts.

**Check:** compare row counts before/after collapse.

### 2.2 Drop all-zero genes

**Behavior:** compute `rowSums(counts_mat)`, remove rows with sum == 0.

**Why:** genes with zero counts across all samples contain no information and only increase multiple testing burden and computation.

### 2.3 Filter lowly expressed genes (group-aware; chosen option)

**Filter description:** keep genes that have `count >= min_count` (10) in `>= min_samples` (3) replicates **within at least one group** (CN, CH, IRN, IRH). This is the union across groups.

**Why group-aware filter:** retains genes that are reproducibly expressed in at least one condition — this preserves condition-specific genes (e.g., genes only induced under IRH), which are critical for target discovery. A global filter (across all samples) would lose such condition-specific signals.

**Your numbers:** start `62,266` → after dropping all-zero `33,258` → after group-aware filtering `18,090`.

**Outputs:** `counts_filtered.csv` (group-aware filtered counts).

**Biological significance:** filtering reduces noise and increases statistical power. In drug discovery, condition-specific genes (e.g., hypoxia-induced mitochondrial stress genes) might be potential targets — group-aware filtering helps retain them.

---

## STEP 3 — Re-load processed files & sanity re-check

**What:** reload `counts_filtered.csv` and `design_matrix_synced.csv` to ensure eventual artifacts were saved correctly. Re-check integer-like counts, no negatives, and columns match metadata.

**Why:** reproducibility: Phase 5 should be able to restart from these saved inputs without rerunning all pre-processing.

**Checks:** `all_integer` function, `any(is.na())`, `any(counts_mat < 0)`, `identical(colnames(counts_mat), meta$Sample)` and reorder if needed.

**Note:** always set `rownames(meta) <- meta$Sample` before creating `DESeqDataSetFromMatrix()` if matching by rownames.

---

## STEP 4 — Create `DESeqDataSet` (dds) and set factor references

**Core code & concept:**

```r
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = meta,
                              design = ~ Treatment * Oxygen)
```

**Line-by-line explanation:**

* `countData`: an integer matrix (genes × samples). We ensured with `round(as.matrix(...))` that values are integer-like.
* `colData`: sample metadata, where rows are samples and rownames should match column names of `countData`.
* `design = ~ Treatment * Oxygen`: the statistical model. The `*` expands to main effects (`Treatment`, `Oxygen`) plus their interaction (`Treatment:Oxygen`). This allows testing (a) effect of insulin resistance, (b) effect of hypoxia, and (c) whether the insulin resistance effect changes under hypoxia (interaction).

**Why set factor references:**

* The baseline (reference) level is the first factor level and controls interpretation of log2 fold changes. We set:

```r
dds$Treatment <- relevel(dds$Treatment, ref = "Control")
dds$Oxygen    <- relevel(dds$Oxygen, ref = "Normoxia")
```


So fold changes default to `Insulin Resistant vs Control` and `Hypoxia vs Normoxia`.

**Pre-filter:** `dds <- dds[rowSums(counts(dds)) > 0, ]` drops any remaining all-zero genes.

**Save:** `saveRDS(dds, file = "../data/processed/dds_phase4_raw.rds")`.

**Checks:**
- `dim(dds)` (genes × samples), `table(dds$Treatment, dds$Oxygen)` (sample counts per group), `colnames(counts(dds)) == rownames(colData(dds))`.

**Biological significance:** building the right model at this stage sets up correct contrasts in Phase 5 (for example: IR effect across oxygen, oxygen effect across treatment, and the interaction for shared or conditional effects). For a drug discovery study you’ll later extract genes whose differential expression is specifically induced in IRH (insulin resistant cardiomyocytes under hypoxia) — these can be candidate targets.

---

## STEP 5 — Size factors & Variance Stabilizing Transform (VST)

**Code:**
```r
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
````

**What `estimateSizeFactors` does:** computes a sample-level scaling factor (median-of-ratios) to correct for library/compositional differences. Inspect with `sizeFactors(dds)`.

**What `vst()` does:** transforms raw counts into continuous values where variance does not depend on mean — this is ideal for PCA, clustering and plotting. Use `blind=FALSE` to preserve biological differences while estimating parameters.

**Outputs saved:** `vsd_phase4.rds` and `vst_matrix.csv` (the matrix returned by `assay(vsd)`).

**Checks:** inspect `dim(vsd)`, `head(assay(vsd))`, and a mean–sd plot if desired.

**Why important for drug discovery:** VST enables reliable identification of clustering patterns and inspection of differential patterns across conditions. It does *not* replace statistical testing on raw counts — but it’s essential for visualization and downstream exploratory analyses (co-expression modules, heatmaps of top DE genes, etc.).

---

## STEP 6 — QC & diagnostics (plots + interpretation)

**Plots created:**

1. **Library sizes barplot** (`library_sizes.png`) — check for samples with low sequencing depth or extreme imbalance.
2. **PCA (PC1 vs PC2)** (`pca_group.png`) — check how samples cluster by Group/Treatment/Oxygen.
3. **Sample distance heatmap** (`sample_distance_heatmap.png`) — hierarchical clustering of samples based on VST distances; checks replicates and outlier detection.
4. **VST Boxplots** (`vst_boxplot.png`) — verifies similarity of distribution after VST; medians and IQRs should be comparable across samples.

**How to interpret (quick guide):**

* Library sizes: typical bulk RNA-seq depth ~20–30M. One sample at ~30M while others ~20M is ok — normalization will account for it.
* PCA: if PC1 or PC2 separates treatment/oxygen groups clearly, your biological conditions produce strong transcriptional signals. Example from our run: PC1 (72%) separated CH (Control-Hypoxia) far from others; PC2 (16%) separated IRH vs CH — indicates hypoxia is strong driver and insulin resistance adds additional effects.
* Heatmap: look for dark-blue blocks along off-diagonal — tight replicates. Orange/yellow between blocks indicates group separation. No single-light-row means no clear outliers.
* Boxplots: aligned medians and similar spread across samples after VST indicate successful normalization.

**Biological consequence:** These QC checks confirm that your data is reproducible and the observed structure is biological (hypoxia and insulin resistance), meaning DE results in Phase 5 are more likely to be valid leads for target discovery.

---

## Provenance & reproducibility notes

* Saved artifacts: `dds_phase4_raw.rds`, `vsd_phase4.rds`, `counts_filtered.csv`, `counts_synced.csv`, `design_matrix_synced.csv`, `vst_matrix.csv`, QC PNGs.
* Save `sessionInfo()` to `data/processed/sessionInfo.txt` to document package versions (critical for reproducibility).
* Store large `.rds` files with Git LFS if committing to GitHub.

**Git commit example:**

```bash
git add reports/Phase4_DataCleaning_DESeq2_Prep.md
git add data/processed/design_matrix_synced.csv data/processed/counts_synced.csv data/processed/counts_filtered.csv
git commit -m "docs(phase4): cleaning, filtering, VST & QC outputs for GSE288708"
git push origin YOUR_BRANCH
```

---

## Biological & theoretical notes (T2D cardiomyopathy context — concise)

* **Biological context:** Type 2 diabetes (insulin resistance) causes metabolic remodeling in cardiomyocytes: mitochondrial dysfunction, altered substrate preference (fatty acid vs glucose), oxidative stress, calcium-handling defects, and increased susceptibility to hypoxia. These processes contribute to diabetic cardiomyopathy (DCM) and heart failure risk.

* **Hypoxia relevance:** in ischemic or poorly perfused myocardium, hypoxia triggers HIF signaling, metabolic reprogramming, and stress responses that intersect with insulin resistance—this joint effect is central to the model (IR ± hypoxia).

* **Why DE analysis matters for drug discovery:** DE genes that are reproducibly altered in IR or IR+Hypoxia could be upstream regulators, signaling nodes, or enzymes amenable to therapeutic modulation (drug targets). Following DE, prioritize genes by effect size, consistency, pathway membership (mitochondrial metabolism, HIF targets, ECM/remodeling, apoptosis, inflammation), and druggability (membrane localization, known inhibitors).

---

## Quick Questions & Answers (Phase 4 topics)

**Q: Why must DESeq2 get raw integer counts?**
A: DESeq2 models counts using the Negative Binomial; it expects discrete counts to estimate dispersion and mean–variance. TPM/FPKM or normalized values break model assumptions.

**Q: Why remove all-zero genes?**
A: They have zero variance — no information — only increase multiple testing burden and slow computations.

**Q: Why group-aware filtering?**
A: It keeps condition-specific genes (e.g., genes that are only highly expressed in IRH replicates) while removing scattered low counts that aren’t reproducible.

**Q: What is library size and why normalize?**
A: Library size = total reads per sample. Normalization (size factors) adjusts for these differences so counts are comparable across samples.

**Q: When should I re-run Phase 4?**
A: If you add/remove samples, correct metadata, or change mapping/filtering procedures — re-run Phase 4 and version artifacts.

---

## Common errors & fixes (quick)

* **`library(vsn)` fails:** install `vsn` inside the R used by your Jupyter kernel via BiocManager.
* **`ggsave()` directory missing:** ensure `dir.create("results/qc", recursive=TRUE)` executed before saving.
* **Mismatched sample names:** set `rownames(meta) <- meta$Sample` and reorder count columns to `meta$Sample` before `DESeqDataSetFromMatrix()`.
* **Non-integer counts:** verify you loaded raw counts, not TPM/FPKM. If tiny floating noise exists (e.g., `100.0000001`), rounding is acceptable; otherwise re-extract raw counts.

---

## Phase 5 preview — what you will do next

1. Fit the model: `dds <- DESeq(dds)`.
2. Check `resultsNames(dds)` to identify coefficient names for contrasts.
3. Extract contrasts:

   * `results(dds, contrast = c("Treatment","Insulin Resistant","Control"))`
   * `results(dds, contrast = c("Oxygen","Hypoxia","Normoxia"))`
   * Interaction or group-specific contrasts using `name=` with `resultsNames(dds)`.
4. Shrink LFCs (`lfcShrink`, prefer `apeglm` or `ashr` where available).
5. Create ranked gene lists and run pathway enrichment (GSEA/ORA). Tools: `clusterProfiler`, `fgsea`.
6. Prioritize druggable targets via cross-referencing DGIdb/DrugBank/ChEMBL, subcellular localization, and literature.

---

> Raw gene-level counts were loaded and aligned to sample metadata. We removed features with zero counts across all samples and applied a group-aware low-count filter (genes with ≥10 counts in ≥3 replicates within at least one experimental group), yielding 18,090 genes for analysis. A `DESeqDataSet` was constructed with design `~ Treatment * Oxygen`. Size factors were estimated and counts were variance-stabilized (VST) for visualization. Quality control included library size inspection, PCA on VST-transformed counts, sample-to-sample distance heatmap, and boxplots of VST distributions. All processed inputs and the `DESeqDataSet` were saved for reproducibility (files: `counts_filtered.csv`, `dds_phase4_raw.rds`, `vsd_phase4.rds`).

---

## Appendix — useful commands & checks

* Check dimensions: `dim(counts_filtered)` | `nrow(dds)` | `ncol(dds)`.
* Confirm factor levels: `levels(colData(dds)$Treatment)`.
* Show size factors: `sizeFactors(dds)`.
* Quick PCA data: `pcaData <- plotPCA(vsd, intgroup=c("Group","Treatment","Oxygen"), returnData=TRUE)`.
* List saved QC files: `list.files("results/qc")`.

---

## Final remarks

Phase 4 is complete. We now have a validated, QC-passed dataset that is ready for differential expression testing and downstream target prioritization. All saved artifacts are kept in `data/processed/` and `results/qc/`, record `sessionInfo()` and package versions, and use Git/LFS for large files.

