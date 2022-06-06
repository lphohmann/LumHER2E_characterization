# Script: DE analysis based on FPKM data based on scanb and metabric

# TODO next:

# empty environment
rm(list=ls())

# set working directory to where the data is
setwd("~/Desktop/MTP_project")

# packages
library(tidyverse)
library(matrixStats)
library(pheatmap)
library(Hmisc)
library(VennDiagram)
library(readxl)
library(ggfortify)
library(janitor)
library(biomaRt)

#######################################################################
# 1. set parameters
#######################################################################

# 1.1 input for which cohort the analysis should be run
cohort <- "SCANB" # Metabric or SCANB or TCGA

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

if (cohort=="SCANB") {
    # load annotation data
    load("Data/SCAN_B/fu.all_v1.1.Rdata")
    # load gex data
    load("Data/SCAN_B/lib_corrected/genematrix_noNeg.Rdata")
    # only include clinical ER+Her2- samples that are "LumA", "LumB", or "Her2"
    anno <- fu.all_v1.1 %>% filter(NCN.PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(!grepl('ERpHER2nLNn_None', evalGroupCM)) %>% filter(grepl('ERpHER2n', evalGroupCM)) %>% dplyr::rename(sampleID = rba, PAM50 = NCN.PAM50, fuV8 = flag_FUv8) %>% filter(fuV8==1) %>% dplyr::select(sampleID, PAM50, NHG)
    genematrix_noNeg <- as.data.frame(genematrix_noNeg)
    gex_data <- genematrix_noNeg[,colnames(genematrix_noNeg) %in% anno$sampleID]
    # remove the other files
    rm(fu.all_v1.1, genematrix_noNeg)
    
    # filter genes based on the stdev cutoff
    
    # 1.2 set stdev filtering cutoff value 
    stdev_cutoff <- 0.5
    
    # 1. log transformed FPKM data
    gex_data_log <- as.data.frame(log2(gex_data + 1))
    
    # filter based on stdev cutoff (do this before or after log transform?)
    gex_data_log <- gex_data_log %>% mutate(stdev=rowSds(as.matrix(.[colnames(gex_data_log)]))) %>% filter(stdev >= stdev_cutoff) %>% dplyr::select(-c(stdev))
    
    # mean center
    gex_data_log <- as.data.frame(t(apply(gex_data_log, 1, function(y) (y - mean(y)) / sd(y))))
    
    View(gex_data_log)
    
} else if (cohort=="Metabric") {
    # load annotation data
    load("Data/Metabric/Annotations/Merged_annotations.RData")
    anno <- as.data.frame(anno) %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(grepl('ERpHER2n', ClinGroup)) %>% dplyr::rename(sampleID=METABRIC_ID,NHG=Grade) %>% dplyr::select(sampleID,PAM50,NHG)
    # load gex data
    # load gex data
    gex_data <- as.data.frame(read.table("Data/Metabric/data_mRNA_median_all_sample_Zscores.txt", sep="\t")) %>% row_to_names(row_number = 1) # samples = 1906 - genes = 24368
    
    # remove rowns with na and duplicates
    gex_data <- na.omit(gex_data)
    gex_data <- gex_data[!duplicated(gex_data$Hugo_Symbol),]
    
    # also save the gene annotation data
    gene_anno <- gex_data[,1:2]
    gex_data <- gex_data %>% remove_rownames() %>% column_to_rownames(var="Hugo_Symbol") %>% dplyr::select(-Entrez_Gene_Id)
    gex_data <- gex_data[,colnames(gex_data) %in% anno$sampleID]
    
    # remove samples from anno that are not in the gex data
    anno <- anno %>% filter(sampleID %in% colnames(gex_data))

    # find the ensembl id corresponding to hugo symbol
    # apprach 1
    #BiocManager::install("AnnotationDbi")
    #library("AnnotationDbi")
    #BiocManager::install("org.Hs.eg.db")
    #library("org.Hs.eg.db")
    #gene_anno$ensembl_gene_id <- mapIds(org.Hs.eg.db,
    #                                    keys=gene_anno$Hugo_Symbol, 
    #                                    column="ENSEMBL",
    #                                    keytype="SYMBOL",
    #                                    multiVals="first")
    #sum(is.na(gene_anno$ensembl_gene_id)) # 8683 not matched genes
    
    # approach 2

    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org", ensemblRedirect = FALSE)
    
    #mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
    #                dataset = "hsapiens_gene_ensembl",
    #                    host = "www.ensembl.org")

    result <- getBM(filters = "hgnc_symbol",
                     attributes = c("hgnc_symbol","ensembl_gene_id"),
                     values = gene_anno$Hugo_Symbol, 
                     mart = mart)
    
    result <- result %>% dplyr::rename(Hugo_Symbol = hgnc_symbol)
    
    length(gene_anno$Hugo_Symbol)-length(result$ensembl_gene_id) #6578 not matched ids
    # merge
    gene_anno <- merge(gene_anno, result, by = "Hugo_Symbol")
    save(gene_anno, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/gene_anno.RData",sep = ""))
    # later i have to exclude non shared genes between the cohorts to define the core set
    # get correct format of files
    gex_data <- gex_data %>% mutate(across(where(is.character), as.numeric)) 
    
    # rename so i can unse the same script as for scanb
    gex_data_log <- gex_data
    #View(head(gex_data_log))

}


# 2. sample annotation
sample_info <- anno[c("sampleID","PAM50")]
#save(sample_info, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/DE_anno.RData", sep = ""))

#######################################################################
# 4. Perform t-test for each gene
#######################################################################

Her2_cols <- sample_info %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA_cols <- sample_info %>% filter(PAM50=="LumA") %>% pull(sampleID)
LumB_cols <- sample_info %>% filter(PAM50=="LumB") %>% pull(sampleID)

# initialize vector of stored p-values and expression differences (based on mean comparison)
H2vsLA_pvalue <- rep(0,nrow(gex_data_log))
H2vsLB_pvalue <- rep(0,nrow(gex_data_log))
H2vsLA_diff <- rep(0,nrow(gex_data_log))
H2vsLB_diff <- rep(0,nrow(gex_data_log))

# loop through the genes
# add progress bar
pb = txtProgressBar(min = 0, max = nrow(gex_data_log), initial = 0, style = 3) 
for (i in 1:nrow(gex_data_log)) {
    setTxtProgressBar(pb,i)
    # set vars
    hdata <- as.numeric(gex_data_log[i,Her2_cols])
    adata <- as.numeric(gex_data_log[i,LumA_cols])
    bdata <- as.numeric(gex_data_log[i,LumB_cols])
    # for Her2 vs LumA
    # equal variance check
    if (var.test(unlist(hdata),unlist(adata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = FALSE)
    } else {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = TRUE)
    }
    # save results
    H2vsLA_pvalue[i] <- H2vsLA_ttest_result$p.value
    H2vsLA_diff[i] <- H2vsLA_ttest_result$estimate[1]-H2vsLA_ttest_result$estimate[2]
    
    # for Her2 vs LumB
    # equal variance check
    if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
    } else {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
    }
    # save results
    H2vsLB_pvalue[i] <- H2vsLB_ttest_result$p.value
    H2vsLB_diff[i] <- H2vsLB_ttest_result$estimate[1]-H2vsLB_ttest_result$estimate[2]
    close(pb)
}

#######################################################################
# 5. Process the output
#######################################################################

# create the final output
results <- gex_data_log %>% add_column(H2vsLA_pvalue = H2vsLA_pvalue,H2vsLB_pvalue = H2vsLB_pvalue,H2vsLA_diff = H2vsLA_diff,H2vsLB_diff = H2vsLB_diff) %>% dplyr::select(H2vsLA_pvalue,H2vsLB_pvalue,H2vsLA_diff,H2vsLB_diff)

# adjust p value 
results$H2vsLA_padj <- p.adjust(results$H2vsLA_pvalue, "fdr")
results$H2vsLB_padj <- p.adjust(results$H2vsLB_pvalue, "fdr")

# up or down (refers to the situation in lumHer2 subtype; e.g. "up" indicates a gene whose expression is upregulated in lumHer2 compared to lumB or lumA)
results <- results %>% mutate(H2vsLA_de =
                                  case_when(H2vsLA_diff <= 0 ~ "down",
                                            H2vsLA_diff >= 0 ~ "up"))
results <- results %>% mutate(H2vsLB_de =
                                  case_when(H2vsLB_diff <= 0 ~ "down",
                                            H2vsLB_diff >= 0 ~ "up"))

# save the file
save(results, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/complete_DE_results.RData",sep = ""))
#load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))

# final result 
# Her2 vs LumA
H2vsLA_DEGs <- results %>% filter(H2vsLA_padj <= 0.05) %>% dplyr::select(H2vsLA_de, H2vsLA_padj) %>% rownames_to_column("Gene") #%>% pull(Gene)
H2vsLA_DEGs_up <- results %>% filter(H2vsLA_padj <= 0.05) %>% filter(H2vsLA_de == "up") %>% rownames_to_column("Gene") #%>% pull(Gene)
H2vsLA_DEGs_down <- results %>% filter(H2vsLA_padj <= 0.05) %>% filter(H2vsLA_de == "down") %>% rownames_to_column("Gene") #%>% pull(Gene)

# Her2 vs LumB
H2vsLB_DEGs <- results %>% filter(H2vsLB_padj <= 0.05) %>% dplyr::select(H2vsLB_de, H2vsLB_padj) %>% rownames_to_column("Gene") #%>% pull(Gene)
H2vsLB_DEGs_up <- results %>% filter(H2vsLB_padj <= 0.05) %>% filter(H2vsLB_de == "up") %>% rownames_to_column("Gene") #%>% pull(Gene)
H2vsLB_DEGs_down <- results %>% filter(H2vsLB_padj <= 0.05) %>% filter(H2vsLB_de == "down") %>% rownames_to_column("Gene") #%>% pull(Gene)

# compare overlap
venn.diagram(
    x = list(H2vsLB_DEGs$Gene, H2vsLA_DEGs$Gene),
    category.names = c("H2vsLB" , "H2vsLA"),
    filename = 'venn_diagramm.png',
    output=TRUE)

# #######################################################################
# # 6. Prepare export for the functional enrichment analysis
# #######################################################################
# 
# # i think i have to strip the endings of the ids because the webtools dont recognize these identifiers
# H2vsLA_DEGs_up$Gene <- gsub("\\..*","",H2vsLA_DEGs_up$Gene)
# H2vsLA_DEGs_down$Gene <- gsub("\\..*","",H2vsLA_DEGs_down$Gene)
# H2vsLB_DEGs_up$Gene <- gsub("\\..*","",H2vsLB_DEGs_up$Gene)
# H2vsLB_DEGs_down$Gene <- gsub("\\..*","",H2vsLB_DEGs_down$Gene)
# 
# # export for functional enrichment analysis
# #write.table(H2vsLB_DEGs_up$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLB_DEGs_up.txt",sep = ""), sep="\t", row.names = FALSE)
# #write.table(H2vsLB_DEGs_down$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLB_DEGs_down.txt",sep = ""), sep="\t", row.names = FALSE)
# #write.table(H2vsLA_DEGs_up$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLA_DEGs_up.txt",sep = ""), sep="\t", row.names = FALSE)
# #write.table(H2vsLA_DEGs_down$Gene, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/H2vsLA_DEGs_down.txt",sep = ""), sep="\t", row.names = FALSE)
# 
# #######################################################################
# # 7. PCA
# #######################################################################
# 
# pca_pam50anno <- sample_info %>% column_to_rownames(var="sampleID") %>% dplyr::rename(PAM50_subtype = PAM50)
# pca_data_analysis <- t(gex_data_log)
# pca_data_plotting <- merge(pca_data_analysis,pca_pam50anno,by=0) %>% column_to_rownames(var = "Row.names")
# 
# pca_res <- prcomp(pca_data_analysis) #, scale=T
# pdf(file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/pca_plot.pdf",sep =""))
# autoplot(pca_res, data=pca_data_plotting, colour="PAM50_subtype", main= "Principal component analysis of gene expression data") +
#     theme(legend.position=c(0.9,0.9))
# dev.off()
# 
# #######################################################################
# # 8. compare bonferroni vs. fdr DEGs
# #######################################################################
# 
# # load scan-b fdr DEGs
# output_dir <- "SCANB"
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))
# scanb_DEGs <- results
# scanb_DEGs <- scanb_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,
#                                     H2vsLB_padj_fdr = H2vsLB_padj)
# 
# # load tcga fdr DEGs
# output_dir <- "TCGA"
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))
# tcga_DEGs <- results
# tcga_DEGs <- tcga_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,
#                                   H2vsLB_padj_fdr = H2vsLB_padj)
# 
# # add bf adjustments
# # tcga
# tcga_DEGs$H2vsLA_padj_bf <- p.adjust(tcga_DEGs$H2vsLA_pvalue,"bonferroni")
# tcga_DEGs$H2vsLB_padj_bf <- p.adjust(tcga_DEGs$H2vsLB_pvalue, "bonferroni")
# # scanb
# scanb_DEGs$H2vsLA_padj_bf <- p.adjust(scanb_DEGs$H2vsLA_pvalue,"bonferroni")
# scanb_DEGs$H2vsLB_padj_bf <- p.adjust(scanb_DEGs$H2vsLB_pvalue, "bonferroni")
# 
# # DEGs
# H2vsLA_DEGs_tcga_fdr <- tcga_DEGs %>% dplyr::filter(H2vsLA_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_tcga_fdr <- tcga_DEGs %>% dplyr::filter(H2vsLB_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLA_DEGs_tcga_bf <- tcga_DEGs %>% dplyr::filter(H2vsLA_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_tcga_bf <- tcga_DEGs %>% dplyr::filter(H2vsLB_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# 
# H2vsLA_DEGs_scanb_fdr <- scanb_DEGs %>% dplyr::filter(H2vsLA_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_scanb_fdr <- scanb_DEGs %>% dplyr::filter(H2vsLB_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLA_DEGs_scanb_bf <- scanb_DEGs %>% dplyr::filter(H2vsLA_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# H2vsLB_DEGs_scanb_bf <- scanb_DEGs %>% dplyr::filter(H2vsLB_padj_bf <= 0.05) %>% rownames_to_column("Gene") %>% dplyr::pull(Gene)
# 
# # compare overlap
# # her2 vs lumA
# venn.diagram(
#     x = list(H2vsLA_DEGs_scanb_fdr, H2vsLA_DEGs_scanb_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_H2vsLA_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumA - SCANB)")
# 
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_fdr, H2vsLA_DEGs_tcga_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/tcga_H2vsLA_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumA - TCGA)")
# 
# # her2 vs lumB
# venn.diagram(
#     x = list(H2vsLB_DEGs_scanb_fdr, H2vsLB_DEGs_scanb_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_H2vsLB_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumB - SCANB)")
# 
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_fdr, H2vsLB_DEGs_tcga_bf),
#     category.names = c("FDR" , "BONFERRONI"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/tcga_H2vsLB_DEGs.png",
#     main = "Comparison FDR vs. BF DEGs (Her2/LumB - SCANB)")
# 
# # scanb vs tcga DEGs
# # FDR Her2-LumB
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_fdr, H2vsLB_DEGs_scanb_fdr),
#     category.names = c("TCGA DEGs" , "SCAN-B DEGs"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLB_DEGs_fdr.png",
#     main = "Shared Her2/LumB-DEGs between TCGA and SCAN-B (core set)")
# 
# # BF Her2-LumB
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_bf, H2vsLB_DEGs_scanb_bf),
#     category.names = c("TCGA" , "SCANB"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLB_DEGs_bf.png",
#     main = "Comparison Bonferroni-DEGs (Her2/LumB)")
# 
# # FDR Her2-LumA
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_fdr, H2vsLA_DEGs_scanb_fdr),
#     category.names = c("TCGA DEGs" , "SCAN-B DEGs"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLA_DEGs_fdr.png",
#     main = "Shared Her2/LumA-DEGs between TCGA and SCAN-B (core set)")
# 
# # BF Her2-LumA
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_bf, H2vsLA_DEGs_scanb_bf),
#     category.names = c("TCGA" , "SCANB"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/scanb_vs_tcga_H2vsLA_DEGs_bf.png",
#     main = "Comparison Bonferroni-DEGs (Her2/LumA)")
# 
# #######################################################################
# # 9. compare the DEGs identified using count data to the ones 
# # identified using FPKM data
# #######################################################################
# 
# # load tcga fdr DEGs
# output_dir <- "TCGA"
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/complete_DE_results.RData",sep = ""))
# tcga_fpkm_DEGs <- results
# tcga_fpkm_DEGs <- tcga_fpkm_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,
#                                   H2vsLB_padj_fdr = H2vsLB_padj)
# 
# # DEGs
# H2vsLA_DEGs_tcga_fpkm_fdr <- tcga_fpkm_DEGs %>% filter(H2vsLA_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% pull(Gene)
# H2vsLB_DEGs_tcga_fpkm_fdr <- tcga_fpkm_DEGs %>% filter(H2vsLB_padj_fdr <= 0.05) %>% rownames_to_column("Gene") %>% pull(Gene)
# 
# 
# # load count DEGs
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/gex_count_luma_diffExp.RData", sep =""))
# 
# load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/gex_count_lumb_diffExp.RData", sep =""))
# 
# H2vsLA_DEGs_tcga_counts <- res_table_luma %>% filter(padj <= 0.05) %>% dplyr::rename(Gene = genename_luma) %>% pull(Gene)
# H2vsLB_DEGs_tcga_counts <- res_table_lumb %>% filter(padj <= 0.05) %>% dplyr::rename(Gene = genename_lumb) %>% pull(Gene)
# 
# # compare
# venn.diagram(
#     x = list(H2vsLB_DEGs_tcga_counts, H2vsLB_DEGs_tcga_fpkm_fdr),
#     category.names = c("counts" , "FPKM"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/TCGA/count_vs_fpkm_H2vsLB_DEGs.png",
#     main = "TCGA: Comparison DEGs identified using count and FPKM data (Her2/LumB)")
# 
# venn.diagram(
#     x = list(H2vsLA_DEGs_tcga_counts, H2vsLA_DEGs_tcga_fpkm_fdr),
#     category.names = c("counts" , "FPKM"),
#     filename = "~/Desktop/MTP_project/Output/Plots/Transcriptomics/TCGA/count_vs_fpkm_H2vsLA_DEGs.png",
#     main = "TCGA: Comparison DEGs identified using count and FPKM data (Her2/LumA)")
