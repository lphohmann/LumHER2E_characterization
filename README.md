# LumHER2E_characterization
This repository contains the scripts that were used in the analyses of the Bioinformatics master thesis project "**Molecular and clinicopathological characterization of Luminal HER2-enriched breast cancer**".

## Project overview: Background and main findings
**Introduction:** Clinically, breast cancer is divided into subgroups based on the assessment of estrogen receptor (ER), progesterone receptor, and human epidermal growth factor receptor 2 (HER2) status. HER2-/ER+ tumors constitute the luminal subgroup. By using a PAM50-based gene expression profiling assay, luminal tumors can be stratified into the intrinsic molecular subtypes Luminal A, Luminal B and HER2-enriched (HER2E). Despite its association with poor patient outcome, the luminal HER2E subtype currently holds no implications for clinical treatment decisions other than as a marker for aggressive disease due to a lack of insight into its characterizing features. Therefore, the study at hand aimed to comprehensively characterize luminal HER2E breast cancer based on clinicopathological and molecular features to evaluate the response of patients to conventionally administered therapies and identify features that might be used to inform treatment decisions in the future. 

**Results:** The findings of the study at hand demonstrated that luminal HER2E breast cancer is a small but clinically important subgroup that is associated with a significantly faster disease recurrence than other luminal subtypes, regardless of current standard of care treatment. With highly proliferative characteristics akin to the Luminal B subtype, HER2E tumors are characterized by a low expression of ESR1, a high immune response, a high burden of copy number alterations and a high frequency of TP53 mutations.

# Project workflow and associated scripts

## Data cohorts
The analyses in this study were based on the primary breast cancer datasets of two large independent cohorts:
1. The Sweden Cancerome Analysis Network – Breast study (SCAN-B)
* provided clinical and RNA-seq data on 4413 clinically ER+/HER2- cases (*Luminal HER2E:* 79; *Luminal A:* 2856; *Luminal B:* 1249)
* data for SCAN-B was obtained from a submitted study by Staaf et al. (https://www.medrxiv.org/content/10.1101/2021.12.03.21267116v2)
* **Availability:** https://data.mendeley.com/datasets/yzxtxn4nmd/draft?a=efbf838d-a16c-4133-8315-eba181b75f3e

2. The Molecular Taxonomy of Breast Cancer International Consortium (METABRIC) 
* provided data on 1227 clinically HER2-/ER+ cases (*Luminal HER2E:* 58; *Luminal A:* 622; *Luminal B:* 350)
* **Availability:** data was downloaded from the CBioPortal website (https://www.cbioportal.org/) as pre-compiled data

## Clinicopathological variables 
### Description
The clinicopathological characterization of luminal HER2E breast cancer included analyses to determine differences of the clinical variables between luminal breast cancers belonging to the intrinsic molecular subtypes Luminal A, Luminal B, and HER2E. Fisher’s exact tests were employed to compare the differences in tumor grades, lymph node statuses, and HER2-low frequency. The clinicopathological variables age and tumor size were assessed using either a two-sided Student’s or Welch's t-test, depending on whether the assumption of equal variance was fulfilled or not. 
### Associated Scripts
* cp_variable_comparisons.R

## Survival analysis 
### Description
Using a Kaplan-Meier analysis with the log-rank test in combination with a univariate Cox’s proportional hazards model, the difference in recurrence-free intervals (RFI) after treatment between the intrinsic luminal subtypes was investigated. To evaluate the effect of the intrinsic molecular PAM50 subtype on RFI when taking additional variables into consideration, a multivariate Cox’s proportional hazards model including the variables PAM50 subtype, age, lymph node status, tumor size and tumor grade was constructed. 
### Associated Scripts
* surv_analysis.R

## Gene expression analyses
**Data preprocessing:** For all gene expression analyses the SCAN-B cohort data preprocessing consisted of log2-transforming the FPKM gene expression data with an offset of +1 and subsequently gene-wise scaling by applying a z-transformation. The METABRIC data was already preprocessed (including z-transformation) and therefore simply used as deposited.

## Metagene analysis 
### Description
The metagene analysis focused on investigating differences between the luminal subtypes in six transcriptional programs (termed basal, lipid, mitotic checkpoint, immune response, steroid response, and stroma) related to breast-cancer biology. Each of the transcriptional programs was assessed by analyzing the expression of an associated gene set (termed a metagene), which were previously defined by Fredlund et. al..
The basal metagene constituted genes such as basal cell keratins and therefore assessed basal transcriptomic characteristics. The lipid metagene was representative of adipocytic characteristics, the mitotic checkpoint metagene of cell cycle processes, the immune response metagene of immune response processes, and the stroma metagene of extracellular matrix-related processes. Furthermore, the steroid response (SR) metagene was constituted by a set of known ER status-related genes and therefore assessed processes related to the response to ER-targeted endocrine treatment. The metagene score for each sample was calculated as the median of the processed gene expression values. Differences in metagene scores between the intrinsic molecular subtypes were assessed using two-sided Welch's or Student’s t-tests.
### Associated Scripts
* metagene_analysis.R

## Differential gene expression analysis
### Description
The differential gene expression analysis was performed to identify differentially expressed genes between luminal breast cancers belonging to the HER2E and subtypes Luminal A and Luminal B. Two-sided t-tests in combination with false-discovery rate for multiple testing correction were employed to identify genes that differ significantly in their expression between the subtypes. Genes with a corrected p-value  0.05 were considered to be differentially expressed. The core set of differentially expressed genes was defined by genes that were both differentially expressed in the SCAN-B, as well as the METABRIC cohort. 
### Associated Scripts
* DE_analysis.R
* heatmap_gex.R

## Target gene expression analysis
### Description
Analysis of the expression of specifically selected genes (like ERBB2) was performed to address the question of sample misclassification, as well as to validate obtained results and potentially gaining further biological/molecular insight. Two-sided Student’s t-tests was performed to investigate differences between luminal intrinsic subtypes. To assess the number of clinically HER2-negative HER2E cases in SCAN-B that might potentially in fact be HER2-positive we defined misclassified samples as cases with ERBB2 mRNA expression values further than 2.7 standard deviations away from the mean of all HER2E cases.
### Associated Scripts
* singlegene_gex_analysis.R

## Copy number alteration analyses 
### Description
Copy number analyses were performed in the METABRIC cohort, as no data was available for SCAN-B. For 22544 genes deposited copy number data provided a copy number state of neutral (0), gain (1), amplification (2), loss (-1), or low-deletion/homozygous deletion (-2). The total frequencies of copy number alterations between subtypes were compared by summing up all alterations (gain, amplification, or loss) that occurred. In contrast, the comparison of spatial copy number alteration profiles differentiated between gain (copy number state >= 1) and loss (copy number state <= -1) alterations, which were assessed for each genomic position. Furthermore, identification of possible tumor drivers was investigated by assessing driver gain (copy number state > 1) and driver loss (copy number state < -1) for each genomic position. Fisher’s exact tests in combination with false-discovery rate for multiple testing correction were employed to assess significance.
### Associated Scripts
* CN_analyses.R

## Mutational enrichment analysis 
### Description
The enrichment of mutations in the subgroups was investigated to identify potential drivers specific to the luminal HER2E subgroup. For each gene the mutational frequency was calculated for the three subgroups and subsequently Fisher’s exact tests in combination with false-discovery rate for multiple testing correction were employed to assess significance (p  0.05). For SCAN-B expressed somatic variants similar to Brueffer et al. were used, and for METABRIC mutational calls from a targeted NGS DNA-based panel.
### Associated Scripts
* mut_analyses.R


