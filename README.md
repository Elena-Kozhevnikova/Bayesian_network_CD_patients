# Bayesian_network_CD_patients

# :scientist: Studying Gene Regulatory Hierarchy in Crohn’s Disease with a Bayesian Network
> Elena Kozhevnikova  
> Institute of Molecular and Cellular Biology  
> Novosibirsk  
> elena.n.kozhevnikova@gmail.com  
> tg: @ezeste  

______________________________________________________________________________________________
## Summary
Inflammatory bowel diseases (IBD) are chronic inflammatory disorders of the gastrointestinal tract that include Crohn's disease and ulcerative colitis. IBD are chronic idiopathic disorders, known to be influenced by a variety of factors, including genetic predisposition and diet. Recent transcriptome analyses have uncovered thousands of genes associated with IBD revealing its complex genetic landscape. Several studies demonstrate strong connection between transcriptomic and metabolomic data, highlighting the interplay between gene expression and metabolism ([Massimino et al., 2021](https://www.nature.com/articles/s43588-021-00114-y)). Identifying key factors among the thousands of differentially expressed genes is essential for determining potential targets in drug development and advancing personalized healthcare. Here we applied weighted gene co-expression network analysis (WGCNA, [Langerfeld and Horvath, 2008](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)) together with Bayesian network (BN, [Puga et al., 2015](https://www.nature.com/articles/nmeth.3550)) prediction to identify key genetic and metabolic factors involved in regulation of the inflammatory state characteristic to IBD using open data ([Braun et al., 2024](https://www.nature.com/articles/s41467-024-48106-6)). We highlighted a list of genes contributing to actin polymerization and bundling, of which PREX1 was further confirmed using the transcriptome data obtained from _Muc2_ mice.

*Keywords* Crohn's disease, Bayesian Networks, WGCNA, gene regulation, metabolism

### The experimental approach described here was inspired by [Agrahari et al., 2018](https://www.nature.com/articles/s41598-018-24758-5#Sec14). 
### We used [*bnlearn*](https://www.bnlearn.com/) R package for Bayesian network learning and inference.

### Contents:
This repository contains the complete analytical pipeline, from raw data processing to Bayesian Network modeling and validation.

> **[00_Environment_Setup.R](code/00_setup_environment.R)** — Installation of required CRAN and Bioconductor packages (`edgeR`, `WGCNA`, `bnlearn`, etc.).
>
> **[Data_Processing_Bash](code/Data_Processing.md)** — Pre-processing of raw GEO datasets (Israel GSE199906 and China GSE233900) using `awk` and `parallel` to merge individual files into transcript-level TPM matrices.
>
> **[01_DGE_Analysis.R](code/01_DGE_analysis.R)** — Aggregation of transcripts to gene level, sample synchronization, and Differential Gene Expression (DGE) analysis using `edgeR`.
>
> **[02_WGCNA_Modules.R](code/02_WGCNA_modules.R)** — Weighted Gene Co-expression Network Analysis to identify modules associated with clinical traits like Calprotectin.
>
> **[03_BN_Computation.R](code/03_BN_computation.R)** — Construction of Bayesian Networks using bootstrap resampling for both gene-only and phenotype-augmented structures.
>
> **[04_Visualization.R](code/04_Visualization.R)** — High-resolution network rendering using `Rgraphviz` with custom color palettes and structural compression.
>
> **[05_GO_Enrichment.R](code/05_GO_enrichment.R)** — Biological interpretation of network hubs and WGCNA modules via Gene Ontology (GO) enrichment analysis.
>
> **[06_Prediction_Model.R](code/06_Prediction_Model.R)** — Evaluation of predictor influence, cross-cohort validation on the Chinese dataset, and classification performance (Confusion Matrix/ROC).

---

*Data files* used in this project include:  
1.```Israel_metadata.csv``` contain Israeli patients' metadata  
2.```china_metadata.csv``` contain Chinese patient metadata  

*Bigger data files* ```israel_transcript_tpm.csv``` containing Isreali patients' gene expression dataset and ```china_transcript_tpm.csv``` containing Chinese gene expression dataset are available by request. Alternatively, these data can be uploaded from its origin for [Israeli](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199906) and [Chinese](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233900) patients originally described by ([Braun et al., 2024](https://www.nature.com/articles/s41467-024-48106-6)) and processes using **[Data_Processing_Bash](code/Data_Processing.md)**.

