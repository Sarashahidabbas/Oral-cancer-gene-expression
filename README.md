# Paired Differential Gene Expression Analysis of Primary Oral Squamous Cell Carcinoma and Matched Lymph Node Metastases

## Abstract

**Background:** Cervical lymph node metastasis strongly predicts poor outcomes in oral squamous cell carcinoma (OSCC). Transcriptional comparisons between primary tumors and matched metastases help uncover mechanisms of spread and candidate markers while controlling for inter-patient variation.

**Methods:** Public bulk RNA-seq data from GEO series [GSE275870](https://nih.gov) (PubMed ID 39778089) provided 108 tumor samples from 72 patients, yielding 16,044 protein-coding genes after filtering, with paired differential expression evaluated via PyDESeq2 across 11 patients.

**Results:** Paired analysis identified 1,384 significant genes (862 upregulated and 522 downregulated in metastases, including top shifts in *C7*, *IL33*, *ADH1B*, *TIMD4*, and *IL24*).

**Conclusions:** Paired analysis isolates metastasis-associated transcripts while capturing local microenvironmental influences.

## Objective of study
This project compares gene expression in primary oral squamous cell carcinoma (OSCC) tumors and matched lymph node metastases. It uses public RNA-seq data and Python.

## Data

Source: GEO series [GSE275870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275870), file `GSE275870_counts.csv.gz`. The original study is on PubMed (PMID 39778089).

| Item | Count |
| --- | --- |
| Patients | 72 |
| Primary tumor samples | 87 |
| Lymph node metastasis samples | 21 |
| Patients with both sample types | 11 |
| Genes in file | 60,624 |
| Protein coding genes | 19,941 |

The dataset has no normal tissue and no survival data.

## Question

Which genes differ between metastases and primary tumors in the same patient?

## Methods

1. Keep protein coding genes with at least 10 counts in at least 20 samples. This leaves 16,044 genes.
2. Run a PCA on log2 counts per million for all 108 samples.
3. Select the 11 patients with both a primary tumor and a metastasis. Use the first sample of each type.
4. Test differential expression with PyDESeq2. The design is `~ patient + tissue`, so each patient acts as their own control.
5. Call genes significant at adjusted p below 0.05 and absolute log2 fold change above 1.

## Results

### Sample overview


![PCA of all samples](pca_all_samples.png)
Metastases (red) sit toward the left of PC1. Several primary tumors overlap with them.

### Differential expression, metastasis versus primary

![Volcano plot](volcano_metastasis_vs_primary.png) 
The paired test finds 1,384 significant genes. Of these, 862 are higher in metastases and 522 are lower.

| Gene | log2 fold change | Adjusted p |
| --- | --- | --- |
| C7 | 7.90 | 1.5e-20 |
| IL33 | 3.49 | 1.0e-13 |
| ADH1B | 6.07 | 5.1e-12 |
| TIMD4 | 4.67 | 2.1e-11 |
| IL24 | -3.82 | 3.4e-10 |
| ADGRG2 | 3.28 | 1.5e-08 |
| FDCSP | 6.56 | 1.5e-08 |
| PLA2G2D | 4.42 | 3.0e-08 |
| GALNT16 | 2.89 | 4.5e-08 |
| PDK4 | 1.76 | 5.8e-08 |

### Top 30 genes




 ![Heatmap of top 30 genes](heatmap_top30.png) 
Each column is one sample from the 11 paired patients. Values are z-scored per gene.

## Limitations

* Eleven pairs give limited power.
* Metastasis samples come from lymph nodes. Many top genes (C7, ADH1B, TIMD4, FDCSP) mark stromal and immune cells, so part of the signal reflects the surrounding node and not the tumor cells.
* The file has no clinical data, so survival analysis is not possible here.
* Only the first sample per patient and tissue type enters the paired test.

## Next steps

* Estimate immune and stromal content (xCell or ESTIMATE) and add it to the model.
* Run GO and KEGG enrichment with gseapy.
* Build a STRING protein network and rank hub genes.
* Test hub genes as metastasis markers with logistic regression and ROC curves, one sample per patient.
* Validate in TCGA-HNSC, including survival.

## How to run

```bash
pip install pandas numpy matplotlib scikit-learn pydeseq2
mkdir data
# place GSE275870_counts_csv.gz inside data/
python 01_qc_and_paired_de.py
```

##  Results

###  Global Structure and Paired Expression

PC1 accounted for 16.9% of variance and PC2 accounted for 9.7%. Paired testing revealed 1,384 significant transcripts (862 upregulated, 522 downregulated in metastases). Top candidate genes include upregulated markers such as *C7*, *IL33*, *ADH1B*, *TIMD4*, *ADGRG2*, *FDCSP*, *PLA2G2D*, *GALNT16*, and *PDK4*, alongside downregulated *IL24*.


## Repository layout

```
data/        input counts (not tracked)
results/     figures and the full differential expression table
01_qc_and_paired_de.py
README.md
```
##  Future Directions

Future extensions include microenvironment deconvolution (`xCell`, `ESTIMATE`), GO/KEGG pathway enrichment (`gseapy`), `STRING` network analysis, and validation in `TCGA-HNSC` cohorts.

---


##  Data Citation and License

* **Source Dataset:** NCBI GEO accession [GSE275870](https://nih.gov) (PubMed ID: `39778089`).
* **License:** Released under the **MIT License** (`LICENSE`).


