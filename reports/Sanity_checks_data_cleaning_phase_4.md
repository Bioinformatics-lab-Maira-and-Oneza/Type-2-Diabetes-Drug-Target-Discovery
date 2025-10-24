# Phase 4 — Data Cleaning & DESeq2 Preparation (T2D iPSC-CM Project)

**Short summary:**
This document records the code, logic, results and biological reasoning for Phase 4 (data cleaning + DESeq2 preparation) of the Transcriptomics-Based Drug Target Discovery project in iPSC-derived cardiomyocytes (T2D / insulin resistance ± hypoxia). It documents exactly what we ran, why we ran it, the outputs produced, and the choices made (thresholds, group-aware filtering). Save this in `reports/` for reproducibility and method reporting.

---

**Date:** 2025-09-14  
**Project:** Transcriptomics-Based Drug Target Discovery in Diabetic Cardiomyopathy (T2D)

---


## Files / paths used (project layout)

* `data/processed/expression_matrix_aligned.csv` — aligned raw counts (rows = gene IDs, cols = sample GSM IDs; first column `gene`).
* `data/processed/design_matrix.csv` — sample metadata (Sample, Treatment, Oxygen, Replicate, Group).
* `data/processed/counts_synced.csv` — saved synced counts (created during the pipeline).
* `data/processed/design_matrix_synced.csv` — saved synced metadata.
* `data/processed/counts_filtered.csv` — final group-aware filtered counts (analysis-ready).
* Notebook: `notebooks/03_Sanity_checks_data_cleaning_deseq2_prep.ipynb` (R kernel) — interactive record of steps.
* Reports: this file `reports/Sanity_checks_data_cleaning_phase_4.md`.

---

## High-level goals for Phase 4

1. Verify data type and ensure counts are integer read counts suitable for DESeq2.
2. Load and align the counts matrix and metadata; ensure every column (sample) in the counts matrix matches metadata `Sample` exactly and in the same order.
3. Sanity checks: non-negativity, integer-valued counts, no missing values, no duplicated gene IDs (or collapse them if present).
4. Remove uninformative features: drop genes with *all-zero* counts across samples.
5. Filter low-count genes using a biologically sensible rule (we used group-aware filtering: `count >= 10` in `>= 3` replicates within *at least one* group). Save the filtered count matrix for DE testing.
6. Produce provenance files so Phase 5 can start from a clean, validated input.

---

## What we did (record of your run)

**Start (raw)**: 62,266 gene rows in the original aligned matrix.

**After removing all-zero genes**: 33,258.

**After group-aware filtering (≥10 counts in ≥3 samples within at least one group)**: **18,090 genes** kept.

This is the dataset we exported to `data/processed/counts_filtered.csv` and will use to construct the DESeqDataSet (dds) in the next step.

---

## Step-by-step code flow (explained)

Below are the core chunks you ran and their intent. Each block is followed by a plain-language explanation of the input → operation → output.

### 1) Load metadata and set factor levels

```r
meta <- read_csv(meta_path, show_col_types = FALSE) |>
  filter(!is.na(Sample)) |>
  mutate(
    Sample    = as.character(Sample),
    Treatment = factor(Treatment, levels = c("Control","Insulin Resistant")),
    Oxygen    = factor(Oxygen,    levels = c("Normoxia","Hypoxia")),
    Replicate = factor(Replicate),
    Group     = factor(Group, levels = c("CN", "CH", "IRN", "IRH"))
  )
```

**Why:** ensures metadata types are explicit. `factor(..., levels=...)` sets reference ordering (Control & Normoxia as baselines). This will determine the interpretation of log2 fold changes later.

**Important:** this *does not* reorder your count matrix — it just sets the categorical labels and their order for modeling and plotting.

---

### 2) Load counts and align columns to metadata

```r
counts_df <- read_csv(counts_path, show_col_types = FALSE)
stopifnot("gene" %in% colnames(counts_df))
# Move gene column into rownames and produce counts_mat (data.frame/matrix)
counts_df <- counts_df |> relocate(gene)
gene_ids <- counts_df$gene
counts_mat <- counts_df |> select(-gene) |> as.data.frame()
# Ensure sample names match metadata
stopifnot(setequal(colnames(counts_mat), meta$Sample))
# Reorder count columns to match meta order exactly
counts_mat <- counts_mat[, meta$Sample, drop = FALSE]
rownames(counts_mat) <- gene_ids
```

**What it does:**

* Reads the counts CSV.
* Checks for `gene` column and uses it as row identifiers.
* Validates that all sample column names are present in metadata (`setequal`) and then *reorders* columns to the exact order of `meta$Sample` (so DESeq2's column order will align with `colData`).

**Why careful ordering matters:** DESeq2 matches samples by column order to `colData` (or by sample names if present), so columns must correspond 1:1 to metadata rows.

---

### 3) Sanity checks

```r
# Non-negative
stopifnot(all(counts_mat >= 0, na.rm = TRUE))
# Integers
is_integer_like <- function(x) all(!is.na(x)) && all(abs(x - round(x)) < 1e-8)
all_integer <- all(vapply(counts_mat, is_integer_like, logical(1)))
stopifnot(all_integer)
# No NAs
stopifnot(!any(is.na(counts_mat)))
# Unique gene IDs
stopifnot(!anyDuplicated(rownames(counts_mat)))
```

**Meaning:** DESeq2 expects non-negative integers (counts). We double-check integers (floating point representation can exist but values must be integer-like), and ensure no missing values.

**Pitfall:** if `any(!all_integer)` you must inspect why values are non-integer — maybe TPM/FPKM or normalized values were accidentally supplied instead of raw counts.

---

### 4) Duplicated gene IDs (defensive collapse)

```r
dup_genes <- rownames(counts_mat)[duplicated(rownames(counts_mat))]
if (length(dup_genes) > 0) {
  counts_mat <- as.data.frame(counts_mat) |>
    tibble::rownames_to_column(var = "gene") |>
    group_by(gene) |>
    summarise(across(-gene, ~ sum(.x, na.rm = TRUE))) |>
    tibble::column_to_rownames(var = "gene")
}
```

**Why:** some imports can double-report the same gene (multiple transcript rows, symbol mapping). Collapsing by summing preserves total gene counts (conservative) and removes duplicated rownames.

**Output:** either unchanged if zero duplicates, or a deduplicated `counts_mat` with fewer rows.

---

### 5) Drop all-zero genes

```r
row_sum <- rowSums(counts_mat)
all_zero <- row_sum == 0
counts_nonzero <- counts_mat[!all_zero, , drop = FALSE]
```

**Why:** genes with zero counts in all samples have no statistical information and only increase multiple-testing burden. Remove them.

**Your numbers:** Start 62,266 → after all-zero removal 33,258.

---

### 6) Low-count filtering (two options)

**Option A (global):** keep if `count >= 10` in `>= 3` samples across the entire dataset (anywhere among 20 samples)

```r
min_count <- 10
min_samples <- 3
keep_mask_A <- rowSums(counts_nonzero >= min_count) >= min_samples
counts_filtered_A <- counts_nonzero[keep_mask_A, , drop = FALSE]
```

**Option B (group-aware, chosen):** keep if `count >= 10` in `>= 3` replicates *within at least one group* (CN, CH, IRN, IRH)

```r
groups <- split(colnames(counts_nonzero), meta$Group)
per_group_keep <- lapply(groups, function(cols) {
  rowSums(counts_nonzero[, cols, drop = FALSE] >= min_count) >= min_samples
})
keep_mask_B <- Reduce(`|`, per_group_keep) # union across groups
counts_filtered <- counts_nonzero[keep_mask_B, , drop = FALSE]
```

**Why we chose Option B:** it preserves condition-specific genes (e.g. genes turned on only in IRH replicates) while discarding genes that have weak scattered expression across unrelated samples. For a drug-discovery study where condition-specific signals are key, group-aware filtering is recommended.

**Your numbers:** Option A kept 18,805 genes; Option B kept 18,090 genes. The difference (715 genes) are genes that had scattered counts but were not reproducible within any single condition.

---

## What we achieved and why it matters

* **Inputs validated:** counts are integer-like, non-negative, matrix aligned to metadata, no NAs.
* **Uninformative features removed:** all-zero genes removed (33k → 18k after all filters), improving statistical power.
* **Group-aware filtering retained condition-specific genes** (important to detect disease or hypoxia-specific genes for target discovery).
* **Provenance saved:** `counts_synced.csv`, `design_matrix_synced.csv`, `counts_filtered.csv` saved to `data/processed/` for reproducibility.

Biological consequence: the reduced gene list focuses hypothesis testing on genes that are actually expressed in our iPSC-derived cardiomyocytes and that show reproducible expression patterns in biologically meaningful groups. This reduces false-positive noise and improves the chance of discovering robust T2D/IR- and hypoxia-related targets.

---

## Biological & theoretical notes (T2D cardiomyopathy context)

* **Counts represent** the number of sequencing reads aligned to a gene feature (gene-level). They are discrete integer measures; negative or fractional values are invalid inputs for count-based methods.
* **Library size** (total reads per sample) varies (your data: \~22–30M). Because library sizes differ, normalization (DESeq2 size factor, or CPM scaling) is needed when comparing across samples. Prefiltering (counts ≥ 10) roughly corresponds to a small CPM threshold given these library sizes.
* **Why filter?** Low-count genes inflate multiple testing and produce unstable dispersion estimates. Removing them increases power for detecting real DE genes.
* **Why group-aware?** Condition-specific genes (like those only upregulated in IRH) are crucial in a drug-discovery pipeline — discarding them would remove potential drug targets.
* **T2D cardiomyopathy biology (concise):** insulin resistance and hyperglycemia can cause metabolic inflexibility, mitochondrial dysfunction, increased oxidative stress and susceptibility to hypoxia, altered calcium handling, ER stress and maladaptive remodeling (ECM changes). Genes in metabolism, mitochondrial function, hypoxia response (HIF targets), apoptosis, inflammation, and fibrosis pathways are good candidates to inspect in DE results.

---

## Practical caveats & recommended checks

* Confirm `identical(colnames(counts_filtered), meta$Sample[match(colnames(counts_filtered), meta$Sample)])` before building `dds`.
* If any samples have extremely low library size or are outliers on PCA, investigate mapping stats, RNA quality, and consider removing them or adding covariates.
* Keep a saved copy of the *pre-filter* gene list (for reproducibility & methods). We saved `counts_synced.csv` and `counts_filtered.csv`.
* If you discover a gene of biological importance in the `lost_genes` list, inspect and consider relaxing the threshold for that gene (not globally) or documenting why it was removed.

---

## Biological quiz & facts (study quick prep)

**Q1: Why must DESeq2 get raw integer counts?**
A1: DESeq2 models counts using the Negative Binomial distribution and estimates dispersion directly from integer counts. Using normalized (TPM/FPKM) or log-transformed values will violate model assumptions.

**Q2: Why remove all-zero genes before DE testing?**
A2: All-zero genes contain no information (zero variance) and only increase multiple testing burden. Removing them saves computation and improves power.

**Q3: What is a library size and why does it matter?**
A3: Library size = sum of reads per sample. Samples with different library sizes must be normalized (size factors or CPM) to compare expression levels across samples.

**Q4: What is CPM and when to use it?**
A4: Counts per million (CPM) scales counts by library size. CPM-based filters are library-size aware and are useful when library depth differs a lot.

**Q5: What is group-aware filtering and why did we choose it?**
A5: Group-aware filtering requires genes to be reproducibly expressed in a minimum number of replicates *within a group*. We chose it to keep condition-specific genes (e.g., IRH-specific), which are important for target discovery.

**Q6: What is the role of factor levels in modeling?**
A6: Factor levels set reference categories. The first level is the baseline for contrasts; setting `Control` and `Normoxia` as references ensures log2FC is relative to biologically meaningful baselines.

**Q7: Why do we collapse duplicated gene IDs (by summing)?**
A7: Duplicate rows for the same gene arise from multiple transcript records or mapping artifacts. Summing conservatively preserves total gene count while producing unique rownames.

**Q8: What biological processes are we likely to see change in T2D cardiomyopathy?**
A8: Expect changes in energy metabolism, fatty-acid oxidation vs glycolysis shifts, mitochondrial genes, oxidative stress response, hypoxia pathways (HIF1), calcium handling, ECM remodeling and inflammatory/fibrotic signaling.

---

## ✅Next steps (Phase 4 → Phase 5 transition)

1. Create a `DESeqDataSet` (`dds`) from `counts_filtered` and `meta`:
   `dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_filtered)), colData = meta, design = ~ Treatment * Oxygen)`
2. Compute VST (variance stabilizing transform): `vsd <- vst(dds, blind = FALSE)` and run QC plots (PCA, sample distances, boxplots).
3. If QC is clean, run `dds <- DESeq(dds)` and extract contrasts (Treatment, Oxygen, interaction, group contrasts) with `results()`.
4. Run independent filtering & multiple testing correction (DESeq2 does this automatically in `results()`).
5. Functional follow-up: pathway enrichment (GSEA/ORA), network & target prioritization (drugability, subcellular localization), literature cross-check.

---

