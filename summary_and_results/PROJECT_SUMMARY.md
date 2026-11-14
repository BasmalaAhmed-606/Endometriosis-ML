# Project Summary and Results

**Machine-learning prediction of *LYVE1* and *CCL19* expression in unseen endometriosis patients from single-cell RNA-seq**

| | |
|---|---|
| **Status** | Analysis complete; write-up and poster in progress |
| **Last updated** | 2026-11-14 |
| **Code and data pointers** | See [section 15](#15-repository-map-and-reproducibility) |

This document collects what was done, every main number, the reasons behind each design decision, and the limitations.
Numbers come from the outputs of the R pipeline and the Python notebook. Where a figure or table is stored in the repository, the
file name is given.

---

## Contents

1. [Research question and answer](#1-research-question-and-answer)
2. [Key numbers at a glance](#2-key-numbers-at-a-glance)
3. [Data](#3-data)
4. [Processing pipeline (R, Seurat)](#4-processing-pipeline-r-seurat)
5. [Target screening (feasibility)](#5-target-screening-feasibility)
6. [Machine-learning design](#6-machine-learning-design)
7. [Internal results](#7-internal-results)
8. [Robustness checks](#8-robustness-checks)
9. [ROC and precision-recall analysis](#9-roc-and-precision-recall-analysis)
10. [Genes that drive the predictions](#10-genes-that-drive-the-predictions)
11. [External validation (*LYVE1*)](#11-external-validation-lyve1)
12. [What the results do and do not support](#12-what-the-results-do-and-do-not-support)
13. [Limitations](#13-limitations)
14. [Comparison with the original study's methods and literature check](#14-comparison-with-the-original-studys-methods-and-literature-check)
15. [Repository map and reproducibility](#15-repository-map-and-reproducibility)
16. [Decision log](#16-decision-log)
17. [Not done and recommended next steps](#17-not-done-and-recommended-next-steps)
18. [Glossary](#18-glossary)
19. [References](#19-references)

---

## 1. Research question and answer

**Question.** Can machine learning predict the expression level of selected immune-regulatory genes (*LYVE1* in macrophages,
*CCL19* in perivascular/mural cells) in cells from patients that were **not used for training**, using the expression of other genes?

**Answer: partly.**

- In held-out patients, the models ranked cells by target-gene expression clearly better than simple baselines
  (mean per-patient Spearman about 0.60 for *LYVE1* and 0.57 for *CCL19*).
- The *LYVE1* result was supported in an independent cohort (mean Spearman 0.55 in five patients).
- Accuracy is moderate, varies strongly between patients, and the number of patients is small.
- *LYVE1* was predictable in several tissues; *CCL19* only where it is expressed (mainly peritoneal lesions). *CCL19* relied heavily
  on *CCL21* and was not validated externally.
- All predictions are **associations, not causal mechanisms**.

**Why this question.** Tan et al. (2022) described the cell populations of eutopic and ectopic endometriosis tissue, including a
perivascular population expressing *CCL19* and *LYVE1*-expressing macrophages, and focused on finding populations and gene
programmes. This project asks a different question: can a model predict the expression of a chosen gene in cells from patients it
has never seen? This gap was proposed by the team; it is not a limitation stated by the original authors, and novelty has not been
confirmed (see [section 14](#14-comparison-with-the-original-studys-methods-and-literature-check)).

---

## 2. Key numbers at a glance

| | *LYVE1* (macrophages) | *CCL19* (mural cells) |
|---|---|---|
| Cells used | 4,762 | 2,537 |
| Cells expressing the gene | 1,431 (30.1%) | 495 (19.5%) |
| Test patients | 10 internal, 5 external | 7 internal, none external |
| Mean per-patient Spearman, Elastic Net | **0.60** internal, **0.55** external | **0.57** internal |
| Baseline used | 0.48 (subtype mean); 0.14 (technical metrics only, external) | 0.11 (subtype mean) |
| ROC AUC (pooled test patients) | 0.83 to 0.84 internal; 0.86 external | 0.91 to 0.92 internal |
| Main caveat | Clusters partly circular; test patients selected for detectable *LYVE1* | Depends strongly on *CCL21*; only 7 patients |

> **How to read the baselines.** The baseline type differs between rows (subtype mean versus technical metrics only), so the gaps
> should not be compared across rows. Averages of per-patient values (Spearman) and pooled values (AUC) are different summaries and
> are labelled as such throughout this document.

---

## 3. Data

| Dataset | Content | Use |
|---|---|---|
| GSE179640 (Tan et al., 2022) | 39 files = 33 `.h5` + 2 samples as `mtx` (6 files); 35 samples | Training and internal evaluation |
| GSE213216 (Fonseca et al., 2023) | 45 single-cell tissue samples from 19 patients (about 15 GB); 16 samples from 14 patients screened | External validation (*LYVE1*) |
| GSE266265 (endometrial atlas) | One 38.9 GB `.rds` file | Considered, not used (too large for the available machine) |

**Samples used from GSE179640:** 31 samples from 14 individuals.

| Tissue | Abbreviation | Samples |
|---|---|---|
| Control endometrium | Ctrl | 3 |
| Eutopic endometrium | EuE | 9 |
| Ectopic peritoneal lesion | EcP | 8 |
| Adjacent peritoneal tissue | EcPA | 7 (includes one technical replicate of patient E09) |
| Ectopic ovarian lesion | EcO | 4 |

**Excluded samples.** Four EOR samples were excluded: EOR01 and EOR03 had zero detected *CCL19*/*LYVE1* cells; EOR02 and EOR04 did not
contain the genes by name and had unusual cell numbers (100,000 and 249,996). The reason was not investigated.

---

## 4. Processing pipeline (R, Seurat)

| Step | Result |
|---|---|
| Cells before QC | 140,433 |
| QC per sample (at least 500 genes, at least 1,000 UMIs, at most 100,000 UMIs, below 25% mitochondrial) | 111,060 |
| Doublets removed (scDblFinder) | 6,717 removed, **104,343 cells** |
| Genes | 27,722 shared by all samples, detected in at least 3 cells |
| Normalization | To the median UMI count (8,855), then log; 1,937 variable genes after removing mitochondrial, ribosomal, hemoglobin, cell-cycle and stress genes |
| Dimensions and batch correction | 30 principal components; Harmony on patient (clustering only); 23 clusters (resolution 0.3) |

**Practical decisions.** The 16 GB memory limit made the standard `merge()` fail (`std::bad_alloc`), so QC was done per sample and
count matrices were combined with `cbind()`. One sample (GSM6595261, patient E10) used a different gene annotation, hence the
common-gene set. Patient and tissue labels come from file names because the GEO metadata download failed.

### 4.1 Major cell types

Assigned manually from marker DotPlots. *CCL19* and *LYVE1* were **not** used as markers.

| Type | Cells |
|---|---|
| Stromal | 37,564 |
| Lymphocytes | 21,176 |
| Epithelial | 16,725 |
| Myeloid | 11,504 |
| Endothelial | 10,551 |
| Unassigned | 6,823 |

### 4.2 Macrophages (myeloid re-clustering)

- Myeloid cells were re-clustered into 14 clusters. Macrophages are clusters 1, 2, 5, 6 and 8.
- Cluster 2 (1,289 cells; FOLR2+, MRC1+; 55.7% *LYVE1*+) is called **Mac-LYVE1**; the other macrophages (3,473 cells) are
  **Other macrophage**.

### 4.3 Mural cells (stromal re-clustering)

- Stromal cells were re-clustered into 16 clusters. Mural clusters are 5 and 12.
- Cluster 5 (2,159 cells; MCAM, RGS5, GGT5, ABCC9, STEAP4 high; 22.1% *CCL19*+) is called **Prv-CCL19** in the code. It should be
  described as a **perivascular cluster**, because it also contains cells from control and ovarian tissue where *CCL19* is almost
  absent. Cluster 12 (378 cells) is **Other mural**.

> Cluster numbers come from one run and must be re-checked against the DotPlots if the pipeline is repeated.

---

## 5. Target screening (feasibility)

A patient counts as **usable** for a group if it has at least 50 cells and at least 10 expressing cells. These thresholds were taken
from the pipeline template; they are not a standard. Their purpose is to make a per-patient correlation meaningful.

| Gene | Group | Cells | Expressing | Usable patients |
|---|---|---|---|---|
| *LYVE1* | All macrophages | 4,762 | 1,431 (30.1%) | 10 |
| *LYVE1* | Mac-LYVE1 only | 1,289 | 718 (55.7%) | 4 |
| *CCL19* | All mural | 2,537 | 495 (19.5%) | 7 |
| *CCL19* | Prv-CCL19 only | 2,159 | 478 (22.1%) | 7 |
| *CCL19* | All stromal | 37,564 | 569 (1.5%) | 7 |

**Usable patients per tissue** (all macrophages for *LYVE1*, all mural cells for *CCL19*):

| | Ctrl | EuE | EcP | EcPA | EcO |
|---|---|---|---|---|---|
| *LYVE1* | 0 | 4 | 8 | 3 | 4 |
| *CCL19* | 1 | 0 | 3 | 5 | 0 |

*LYVE1* can be modelled in several tissues; *CCL19* only in a few patients, mainly peritoneal lesions.

**Decision.** The ML uses the **broad groups** (all macrophages, all mural cells) rather than Mac-LYVE1 or Prv-CCL19, because those
clusters were partly defined by the target genes themselves.

---

## 6. Machine-learning design

| Element | Choice | Reason |
|---|---|---|
| Task | One regression model per target gene | Each gene has its own cell group and its own question |
| Models | Elastic Net, Random Forest, XGBoost (fixed settings) | A linear, a bagged-tree and a boosted-tree model, to see whether the result depends on the model family |
| Evaluation | Leave-one-patient-out; every test patient is completely absent from training; only usable patients are tested, the others stay in training | Cells from the same patient resemble each other; a random split would overstate performance |
| Inputs | Log-normalized expression; the Harmony embedding is not used | Harmony removes patient differences, which would hide exactly what is being tested |
| Target gene | Removed from the inputs | Prevents the trivial answer |
| Feature selection | 2,000 most variable genes, chosen on training cells only | Selecting on all cells would leak test information |
| Main metric | Spearman correlation | Expression is zero-inflated and the interest is in ranking cells |
| Other metrics | R², RMSE, MAE, bias, AUC (expressing versus zero), Spearman among expressing cells | Describe absolute accuracy, calibration and the zero/non-zero split |
| Baselines | Subtype mean, cluster mean, technical metrics only (UMI count, genes detected, mitochondrial fraction) | Show whether the models learn more than cell type or cell quality |

The Python notebook is `scripts/03_ml_python/03_ml_prediction_and_external_validation.ipynb` (Google Colab). It reproduced the earlier
numbers exactly when re-run from the exported data.

---

## 7. Internal results

Mean over the test patients:

| Analysis | Model | Spearman | R² | RMSE | MAE | AUC | Spearman among expressing |
|---|---|---|---|---|---|---|---|
| *LYVE1*, all macrophages (10 patients) | Subtype mean (baseline) | 0.478 | -0.033 | 0.991 | 0.744 | 0.722 | 0.385 |
| | Elastic Net | 0.601 | 0.402 | 0.698 | 0.540 | 0.864 | 0.689 |
| | Random Forest | 0.578 | 0.408 | 0.720 | 0.535 | 0.844 | 0.699 |
| | XGBoost | 0.597 | 0.471 | 0.653 | 0.468 | 0.859 | 0.710 |
| *CCL19*, all mural (7 patients) | Subtype mean (baseline) | 0.109 | -0.274 | 1.331 | 0.969 | 0.542 | 0.162 |
| | Elastic Net | 0.569 | 0.342 | 0.876 | 0.593 | 0.901 | 0.475 |
| | Random Forest | 0.555 | 0.249 | 1.016 | 0.662 | 0.892 | 0.399 |
| | XGBoost | 0.553 | 0.338 | 0.914 | 0.556 | 0.890 | 0.439 |

The three models are close. With this few patients and no tuning, **no model is claimed to be better than another**.

### 7.1 Per-patient Spearman, *LYVE1* (macrophages)

| Patient | E01 | E02 | E03 | E04 | E05 | E06 | E07 | E09 | E10 | E11 |
|---|---|---|---|---|---|---|---|---|---|---|
| Subtype mean | 0.124 | 0.553 | 0.663 | 0.506 | 0.424 | 0.617 | 0.469 | 0.566 | 0.727 | 0.132 |
| Elastic Net | 0.596 | 0.472 | 0.381 | 0.623 | 0.618 | 0.798 | 0.579 | 0.763 | 0.862 | 0.318 |
| Random Forest | 0.571 | 0.388 | 0.397 | 0.588 | 0.605 | 0.772 | 0.568 | 0.698 | 0.861 | 0.329 |
| XGBoost | 0.605 | 0.455 | 0.405 | 0.620 | 0.615 | 0.786 | 0.527 | 0.746 | 0.872 | 0.342 |

- The models beat the subtype baseline in **8 of 10 patients**; in E02 and E03 the baseline was higher. Range: 0.32 (E11) to 0.86 (E10).
- E10 has only 54 macrophages.
- E11 is the weakest. All models over-predict its mean expression (bias +0.68 to +0.81), but the simple subtype baseline over-predicts it
  too (+0.58). This suggests a lower *LYVE1* level in E11 cells than in same-type cells of other patients, rather than a model failure
  (cause not tested).

### 7.2 Per-patient Spearman, *CCL19* (mural cells)

| Patient | C02 | E01 | E02 | E03 | E04 | E05 | E06 |
|---|---|---|---|---|---|---|---|
| Subtype mean | 0.034 | 0.035 | 0.160 | 0.339 | 0.180 | -0.107 | 0.120 |
| Elastic Net | 0.225 | 0.609 | 0.594 | 0.771 | 0.487 | 0.789 | 0.504 |
| Random Forest | 0.210 | 0.569 | 0.540 | 0.786 | 0.484 | 0.772 | 0.525 |
| XGBoost | 0.223 | 0.593 | 0.547 | 0.768 | 0.446 | 0.791 | 0.505 |

- All seven patients are above the baseline. Range: 0.23 (C02) to 0.79 (E05).
- C02 is the only control and has just 13 *CCL19*+ cells.
- Patient E03 alone contributes 258 of the 495 *CCL19*+ cells.

### 7.3 Calibration

The models rank cells well but estimate the absolute level less well: predictions shrink toward the mean, and the highest expression
values are under-estimated (`figures/fig2_calibration_deciles.png`).

---

## 8. Robustness checks

| Check | Result |
|---|---|
| *LYVE1* without the *LYVE1*-defined cluster (Other macrophages only; 9 patients) | Elastic Net 0.381, XGBoost 0.371, versus 0.037 for a cluster-mean baseline; inside every remaining cluster the correlation stayed positive (0.29 to 0.55) |
| Technical metrics only, *LYVE1* | Pooled Spearman 0.158 (with genes: 0.584) |
| Technical metrics only, *CCL19* inside Prv-CCL19 | Pooled Spearman 0.085 (with genes: 0.668) |
| *CCL19* without *CCL21* as input | Elastic Net 0.548 (R² 0.076), XGBoost 0.535 (R² 0.190): ranking changed little, absolute level got worse |
| Elastic Net tuned (nested; setting chosen inside the training patients) | *LYVE1* 0.602 (R² 0.405); *CCL19* 0.564 (R² 0.339): no gain over defaults |
| Expressing versus zero (AUC) | *LYVE1* 0.86; *CCL19* 0.89 to 0.90 |
| Ranking among expressing cells | *LYVE1* 0.69 to 0.71; *CCL19* 0.44 to 0.48 |

**What each check answers.**

- *Without the LYVE1-defined cluster* tests whether performance only comes from recognising a cluster that was partly defined by *LYVE1*.
- *Technical metrics only* tests whether cell quality alone explains the result.
- *Without CCL21* tests how much *CCL19* prediction depends on one gene.
- *Nested tuning* tests whether fixed settings were holding performance back.

**Interpretation notes.** The AUC is computed on all cells of the group, so part of it reflects membership of the expression-rich
subtype. The Spearman among expressing cells may partly reflect raw counts and cell depth; this was not tested.

---

## 9. ROC and precision-recall analysis

The models predict a value, so to draw these curves the task is turned into a classification: a cell is "expressing" if its count of
the target gene is above zero, and the predicted expression is the score (no threshold needed). The curves pool the test patients.
This analysis was requested by the supervisor.

**Files:** `figures/fig5_roc_curves.png`, `figures/fig6_precision_recall_curves.png`, `results/05_ml_results/roc_pr_pooled_summary.csv`,
`results/05_ml_results/roc_pr_per_patient.csv`.

| Analysis | Cells (patients) | Expressing | Model | ROC AUC | Average precision | Baseline AUC / AP |
|---|---|---|---|---|---|---|
| *LYVE1* internal | 4,416 (10) | 32% | Elastic Net | 0.833 | 0.715 | 0.529 / 0.397 (subtype mean) |
| | | | XGBoost | 0.839 | 0.722 | |
| *CCL19* internal | 1,802 (7) | 27% | Elastic Net | 0.920 | 0.840 | 0.363 / 0.230 (subtype mean) |
| | | | XGBoost | 0.910 | 0.806 | |
| *LYVE1* external | 1,663 (5) | 28% | Elastic Net | 0.860 | 0.745 | 0.494 / 0.283 (technical metrics) |
| | | | XGBoost | 0.860 | 0.746 | |

**Mean of the per-patient values (Elastic Net):** AUC 0.864 (*LYVE1* internal), 0.901 (*CCL19*), 0.849 (*LYVE1* external); the
per-patient baselines are 0.722, 0.542 and 0.582. Random Forest is close to the other two models (pooled AUC 0.813 and 0.900).

> **Pooled versus per-patient.** Pooled AUC (all test cells together) and the mean of per-patient AUC are different summaries and give
> slightly different numbers. Always state which one is quoted.

**Notes.**

- The average precision depends on the share of expressing cells (the no-skill level equals that share), so it should not be compared
  between genes.
- The pooled subtype-mean baseline for *CCL19* has an AUC below 0.5 because the training mean shifts when a patient is held out; the
  per-patient value (0.54) is the fair comparison (likely explanation, not tested).
- Part of the AUC reflects recognising the expression-rich cell type, and a zero count can be a dropout rather than a true absence.

---

## 10. Genes that drive the predictions

These are the genes with a non-zero Elastic Net coefficient in **every** leave-one-patient-out fold, always with the same sign.
The folds overlap heavily, so this is weak evidence of stability. A large coefficient means association with the target gene, not a cause.

| Target | Genes | Positive, largest | Negative, largest |
|---|---|---|---|
| *LYVE1* | 59 | *RNASE1, LILRB5, COLEC12, MARCO, AGFG1, EMB, CLEC10A, CYBB, FCGR2B, C1QC* | *APOE, GPNMB, ALCAM, LGALS1, IGSF6, HLA-DQA1* |
| *CCL19* | 54 | *CCL21* (coefficient 0.568, about 3.5 times the next gene *CXCL12* at 0.160), then *CFD, DEPP1, TNC, NBL1* | *TCF21, ZFP36, IGFBP2, MT2A, NRP1* |

- *FOLR2* was not selected in any fold for *LYVE1*; *MRC1* and *C1QA* were selected in 9 of 10 folds.
- For *CCL19*, mural markers other than *GGT5* were barely used.
- The original study (preprint) reports that Prv-CCL19 cells also express *CCL21*, consistent with the *CCL21* dependence.
- The biological interpretation of the *LYVE1* genes has not been checked against the literature.

---

## 11. External validation (*LYVE1*)

**Dataset:** GSE213216 (Fonseca et al., 2023). Selected samples were downloaded in three small batches.

**Design, fixed before looking at results**

| Rule | Definition |
|---|---|
| Test patient | At least 10 *LYVE1*+ macrophage-like cells |
| Macrophage-like cell | At least 2 of *C1QA, C1QB, C1QC, CD68, CD163, MS4A7* detected; *LYVE1, FOLR2* and *MRC1* not used |
| Processing | Same QC, doublet removal and normalization scale factor as the training data |
| Models | Trained on all GSE179640 macrophages, fixed settings, applied without any tuning on external data |
| Success criterion | Spearman positive and higher than the technical-only model in at least 4 of 5 patients |

**Screening result.** 16 samples from 14 patients were screened; 5 patients passed (GEO patients 7, 8, 9, 11 and 14, labelled X07 to X14
in the files). This gives 1,663 cells.

| Patient | Cells | *LYVE1*+ | Elastic Net | XGBoost | Technical only |
|---|---|---|---|---|---|
| X07 | 198 | 48 | 0.510 | 0.529 | 0.107 |
| X08 | 346 | 154 | 0.695 | 0.703 | 0.210 |
| X09 | 362 | 138 | 0.592 | 0.577 | 0.077 |
| X11 | 667 | 108 | 0.465 | 0.502 | 0.104 |
| X14 | 90 | 19 | 0.505 | 0.426 | 0.213 |
| **Mean** | | | **0.553** | **0.547** | **0.142** |

- Mean R²: Elastic Net 0.323, XGBoost 0.354 (a constant training-mean prediction gives -0.055). AUC about 0.85; pooled Spearman 0.58.
- **The pre-defined success criterion was met in 5 of 5 patients.**
- Performance was about 0.05 lower than in the internal evaluation, which is expected for independent data.

**Caveats.** Test patients were selected for detectable *LYVE1*; 4 of 5 are endometrioma samples; X14 has only 90 cells; hormonal
treatment and cycle phase are unknown; the macrophage definition differs from the training clusters. *CCL19* was not validated
externally. Figure: `figures/fig3_external_validation.png`.

---

## 12. What the results do and do not support

| The results support | The results do not support |
|---|---|
| In held-out patients, models rank cells by *LYVE1* and *CCL19* expression better than simple baselines | Any causal mechanism: all findings are associations |
| Performance is higher than technical metrics alone (0.16 and 0.09) | Validity of *CCL19* prediction outside the 7 internal patients (no external test) |
| *LYVE1* prediction held in an independent cohort (5 of 5 patients met the pre-defined criterion) | Accurate absolute expression levels (predictions shrink toward the mean) |
| *LYVE1* can be modelled in several tissues; *CCL19* mainly in peritoneal lesions | Use as a diagnostic or clinical tool |
| Performance remains after removing the *LYVE1*-defined cluster (0.38 versus 0.04) | Superiority of any one model over the others |
| | Novelty of the approach (literature check was not systematic) |
| | A single average as a summary: per-patient performance ranges from 0.23 to 0.86 |

---

## 13. Limitations

| # | Limitation | Effect |
|---|---|---|
| 1 | Few patients (10, 7 and 5); per-patient performance varies widely | Averages are uncertain; no confidence intervals were computed |
| 2 | Test patients were selected for detectable expression of the target gene (internally and externally) | Performance may be higher than in unselected patients |
| 3 | External validation covers *LYVE1* only, mostly endometrioma samples, and one patient has only 90 cells | *CCL19* and other tissues are untested externally |
| 4 | Annotation was done on all patients before the train/test split | Possible leakage through cluster definitions |
| 5 | Partial circularity: Mac-LYVE1 and Prv-CCL19 were partly shaped by the target genes | The analysis without the *LYVE1* cluster addresses this for *LYVE1* only |
| 6 | Doublet removal differs from the original study (scDblFinder only) | Residual doublets are possible |
| 7 | Four EOR samples were excluded without finding the reason | Possible selection effect |
| 8 | Random Forest and XGBoost were not tuned (XGBoost tuning was too slow on the available machine) | Model comparison is not definitive |
| 9 | The Elastic Net shrinks predictions toward the mean | The highest expression values are under-estimated |
| 10 | Baselines differ between analyses (subtype mean, cluster mean, technical metrics) | Gaps cannot be compared across analyses |
| 11 | Hormonal treatment and cycle phase are unknown for the external samples | Possible unmeasured confounding |
| 12 | Association, not causation | No mechanistic claim |

---

## 14. Comparison with the original study's methods and literature check

### 14.1 Methods compared with the original study

Only the preprint Methods were read; the published version could not be opened.

| Item | Original (preprint) | This project |
|---|---|---|
| QC thresholds | 500 genes, 1,000 UMIs, 100,000 UMIs, 25% mitochondrial, genes in at least 3 cells | Same |
| Normalization and variable genes | Median UMI, log, 2,000 genes, same exclusions | Same |
| Doublets | Scrublet plus cluster-level filtering | scDblFinder only |
| Harmony | Patient, tissue type and stage | Patient only (clustering only) |
| Clustering | Leiden (Scanpy) | Louvain (Seurat) |
| Annotation | Marker thresholds, statistical | Manual from DotPlots |

The figure of 108,497 cells for the study comes from a secondary table, not the paper itself; this project has 104,343 cells.

### 14.2 Literature check

A short search (two queries, not systematic) found machine-learning work on endometriosis diagnosis (bulk RNA-seq classifiers,
foundation models, a CNN on an atlas) and one study that used GSE179640 for deconvolution and a diagnosis classifier. It did not find
work that predicts the expression of a specific gene in unseen patients on these data. **This is not proof that none exists; a proper
PubMed search is needed before claiming novelty.**

---

## 15. Repository map and reproducibility

### 15.1 Where everything is

| What | Where |
|---|---|
| R pipeline | `scripts/01_r_pipeline/01_endometriosis_pipeline.R` |
| External data preparation | `scripts/02_external_prep/02_external_validation_prep.R` |
| ML notebook (Colab) | `scripts/03_ml_python/03_ml_prediction_and_external_validation.ipynb` |
| Target screening tables | `results/03_target_screening/` |
| ML tables and predictions | `results/05_ml_results/` (summary: `summary_all_analyses.csv`) |
| External validation outputs | `results/06_external_validation/` |
| Figures | `figures/` (final selection in `summary_and_results/figures/`) |
| Early exploratory scripts | `scripts/archive/` |

### 15.2 Figure index

| Figure | File |
|---|---|
| UMAP overview (clusters, tissue, patient) | `03_umap_overview.png` |
| UMAP by major cell type | `05_umap_celltype.png` |
| Major cell-type DotPlot | `04_major_celltype_dotplot.png` |
| Myeloid DotPlot | `09_myeloid_dotplot.png` |
| Stromal DotPlot | `07_stromal_dotplot_v2.png` |
| Per-patient Spearman | `fig1_spearman_per_patient.png` |
| Calibration by prediction decile | `fig2_calibration_deciles.png` |
| External validation | `fig3_external_validation.png` |
| *LYVE1* versus *CCL19* summary | `fig4_LYVE1_vs_CCL19_summary.png` |
| ROC curves | `fig5_roc_curves.png` |
| Precision-recall curves | `fig6_precision_recall_curves.png` |

### 15.3 Run order

1. R pipeline, section by section (restart R between sections on a 16 GB machine).
2. External data preparation (R).
3. ML notebook in Google Colab, run all cells.

### 15.4 Reproducibility notes

- Random seed 42 is used in the R pipeline and the Python notebook.
- The R scripts were assembled from code that was run step by step and **have not been run end to end as single files**.
- Cluster numbers used for annotation come from one run and must be re-checked against the DotPlots.
- The Python notebook was re-run from the exported data and reproduced the earlier numbers exactly.
- Large files (`.rds`, `.mtx`, raw data) are not tracked in Git.
- Patient and tissue labels come from file names (GEO metadata download failed).

---

## 16. Decision log

| Decision | Reason |
|---|---|
| Per-sample QC and `cbind()` instead of `merge()` | The standard merge ran out of memory on 16 GB |
| Common-gene set (27,722 genes) | One sample used a different gene annotation |
| Labels taken from file names | The GEO metadata download failed twice |
| ML on broad groups (all macrophages, all mural cells) | Mac-LYVE1 and Prv-CCL19 were partly defined by the target genes |
| Harmony used for clustering only | The ML must see patient differences, not a corrected embedding |
| *CCL19* and *LYVE1* not used as annotation markers | Avoid circularity |
| GSE266265 rejected as external dataset | 38.9 GB single file, too large for the available machine |
| GSE213216 used, in three small download batches | Fits the available machine |
| Only Elastic Net tuned | XGBoost tuning was too slow; Elastic Net tuning gave no gain |
| External validation rules fixed in advance | Avoid choosing criteria after seeing results |
| ROC and precision-recall curves added | Requested by the supervisor |

---

## 17. Not done and recommended next steps

### 17.1 Not done

- External validation of *CCL19*.
- Systematic literature review (including a check of the *LYVE1*-associated genes).
- Tuning of Random Forest and XGBoost.
- Investigation of the four excluded EOR samples.
- A larger external cohort covering more tissues.

### 17.2 Recommended next steps

| Step | Purpose |
|---|---|
| Bootstrap over patients for confidence intervals of the mean Spearman and AUC | Quantify uncertainty with so few patients |
| Use the same baselines in all analyses (subtype mean and technical metrics, internal and external) | Make gaps comparable |
| Partial Spearman among expressing cells, controlling for UMI count and genes detected | Test whether that correlation only reflects sequencing depth |
| Check that the "same sign in every fold" filter is applied in the code for the stable-gene lists | Make the 59 and 54 gene counts match the definition |
| Save the cluster-to-label tables and a session record (`sessionInfo()`) | Make annotation reproducible |
| Run the full workflow end to end on a clean machine | Confirm that numbers reproduce |

---

## 18. Glossary

| Term | Meaning |
|---|---|
| scRNA-seq | Single-cell RNA sequencing |
| UMI | Unique molecular identifier; roughly one count per original RNA molecule |
| QC | Quality control |
| Doublet | Two cells captured as one |
| Harmony | Method that corrects batch (here patient) differences in an embedding |
| HVG | Highly variable gene |
| LOPO | Leave-one-patient-out: each patient is held out in turn as the test set |
| Spearman correlation | Rank correlation between predicted and observed expression |
| R² | Share of variance explained; negative if worse than predicting the mean |
| RMSE / MAE | Root mean squared error / mean absolute error |
| AUC (ROC) | Probability that a random expressing cell gets a higher score than a random non-expressing cell |
| Average precision | Area under the precision-recall curve |
| Elastic Net | Linear regression with L1 and L2 penalties |
| Baseline | Simple reference prediction the models must beat |
| Circularity | Using information in the analysis that was itself derived from the target |
| Ctrl, EuE, EcP, EcPA, EcO | Control, eutopic endometrium, ectopic peritoneal lesion, adjacent peritoneal tissue, ectopic ovarian lesion |
| EOR | Sample group excluded from this analysis |

---

## 19. References

1. Tan Y, et al. Single-cell analysis of endometriosis reveals a coordinated transcriptional programme driving immunotolerance and
   angiogenesis across eutopic and ectopic tissues. *Nat Cell Biol.* 2022. (GSE179640)
2. Fonseca MAS, et al. Single-cell transcriptomic analysis of endometriosis. *Nat Genet.* 2023. (GSE213216)
