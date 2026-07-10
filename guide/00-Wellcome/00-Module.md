# Python for Single-Cell RNA-seq Bootcamp

*A practical, project-based bootcamp using Python — from programming fundamentals to advanced single-cell RNA-seq analysis.*

---

# Prerequisite Module — Python Fundamentals

## Goal
Build the Python foundation required for modern single-cell RNA-seq analysis.

## Topics

- Python basics
- NumPy
- Pandas
- Matplotlib
- AnnData
- Jupyter Notebook
- Conda environments
- Project organization

---

# Module 1 — Fundamentals of scRNA-seq Analysis (Scanpy)

## Goal
Understand the complete preprocessing workflow using **Scanpy**.

## Practical Topics

### Reading Datasets
- 10x Genomics
- `.h5`
- `.h5ad`
- MTX

### Understanding the AnnData Structure
- `X`
- `obs`
- `var`
- `obsm`
- `layers`
- `uns`
- `raw`

### Quality Control (QC)
- `total_counts`
- `n_genes_by_counts`
- `pct_counts_mt`
- Doublet detection

### Filtering

### QC Visualization
- Violin plots
- Scatter plots
- Histograms

### Normalization
- `normalize_total()`
- `log1p()`
- Pearson Residuals *(optional)*

### Feature Selection
- Highly Variable Genes (HVGs)

### Dimensionality Reduction
- Scaling
- PCA

### Graph Construction
- Neighborhood graph

### Visualization
- UMAP

### Clustering
- Leiden clustering

### Marker Gene Identification
- `rank_genes_groups()`

### Cell Type Annotation
- CellMarker
- PanglaoDB
- CellTypist
- scANVI *(optional)*

## Deliverable

- Fully annotated **AnnData** object

---

# Module 2 — Integration & Batch Effect Correction

## Goal
Integrate multiple datasets while preserving biological variation.

## Topics

### Why Batch Correction Matters

### Integration Methods
- Harmony
- BBKNN
- Scanorama
- MNN
- scVI

### Comparing Integration Methods

### Integration Diagnostics
- UMAP
- Batch mixing
- Biological conservation

### Evaluation Metrics
- kBET
- LISI
- Silhouette score

## Deliverable

- Integrated **AnnData** object

---

# Module 3 — Differential Expression & Biological Interpretation

## Goal
Identify biologically meaningful differential expression.

## Topics

### Differential Expression
- Subsetting cell populations
- Within-cell-type DE
- Across-cell-type DE
- Across-condition DE

### Pseudobulk Analysis
- DESeq2
- edgeR
- PyDESeq2

### Pathway Enrichment
- gseapy
- Gene Ontology (GO)
- KEGG
- Reactome
- GSEA

### Biological Interpretation

## Deliverable

- Differential expression report
- Pathway enrichment analysis

---

# Module 4 — Cell–Cell Communication

## Goal
Infer communication between cell populations.

## Tools

- LIANA
- CellPhoneDB
- CellChat
- NicheNet
- Squidpy *(optional)*

## Deliverable

- Cell–cell communication network

---

# Module 5 — Trajectory & Cell Fate Analysis

## Goal
Reconstruct dynamic biological processes and infer cellular trajectories.

## Tools

- PAGA
- DPT
- Palantir
- Monocle3
- scVelo *(RNA velocity)*

## Topics

- Graph abstraction
- Pseudotime
- Lineage inference
- Branch analysis
- RNA velocity
- Fate probability

## Deliverable

- Trajectory and lineage analysis

---

# Module 6 — Spatial Transcriptomics *(Optional)*

## Topics

- Squidpy
- Visium
- Neighborhood enrichment
- Spatial autocorrelation
- Ligand–receptor analysis in spatial context
- Spatial visualization

---

# Module 7 — Visualization & Publication-Ready Figures

## Topics

### Embedding Visualization
- UMAP customization

### Expression Visualization
- DotPlot
- Heatmap
- MatrixPlot
- ViolinPlot

### Statistical Visualization
- Volcano plots
- Pathway plots

### Figure Preparation
- Figure assembly
- Export for publication

---

# Module 8 — Reproducible Bioinformatics

## Topics

### Environment Management
- Conda
- Mamba
- `environment.yml`

### Version Control
- Git
- GitHub

### Workflow Management
- Snakemake *(optional)*

### Best Practices
- Logging
- Project structure
- Reproducible workflows

---

# Final Project

## Goal

Analyze a complete public **scRNA-seq** dataset from raw sequencing counts to biological interpretation.

## Deliverables

- Clean GitHub repository
- Reproducible analysis pipeline
- Fully annotated AnnData object
- Publication-quality figures
- Final report
- Complete project documentation

---

# Learning Outcome

By the end of this bootcamp, you will be able to:

- Build reproducible Python environments for bioinformatics.
- Process raw scRNA-seq datasets using Scanpy.
- Perform quality control and preprocessing.
- Integrate multiple datasets while correcting batch effects.
- Identify cell types and marker genes.
- Conduct differential expression and pathway enrichment analyses.
- Infer cell–cell communication.
- Perform trajectory and cell fate analyses.
- Analyze spatial transcriptomics datasets.
- Produce publication-quality visualizations.
- Develop fully reproducible bioinformatics workflows.
- Complete an end-to-end single-cell RNA-seq project suitable for research and publication.