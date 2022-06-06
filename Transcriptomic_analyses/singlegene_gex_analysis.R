# Script: ERBB2 and ESR1 assessment in the subtypes lumHer2, lumA, LumB

# TODO next:

# empty environment
rm(list=ls())

# set working directory to where the data is
setwd("~/Desktop/MTP_project")

# packages
#library(matrixStats)
library(readxl)
library(biomaRt)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(janitor)

#######################################################################
# 1. set parameters
#######################################################################

# 1.1 input for which cohort the SA should be run
cohort <- "SCANB" # Metabric or SCANB 

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
    anno <- fu.all_v1.1 %>% filter(NCN.PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(!grepl('ERpHER2nLNn_None', evalGroupCM)) %>% filter(grepl('ERpHER2n', evalGroupCM)) %>% dplyr::rename(sampleID = rba, PAM50 = NCN.PAM50, fuV8 = flag_FUv8) %>% filter(fuV8==1)
    genematrix_noNeg <- as.data.frame(genematrix_noNeg)
    gex_data <- genematrix_noNeg[,names(genematrix_noNeg) %in% anno$sampleID]
    # remove the other files
    rm(fu.all_v1.1, genematrix_noNeg)
    
} else if (cohort=="Metabric") {
    # load annotation data
    load("Data/Metabric/Annotations/Merged_annotations.RData")
    anno <- as.data.frame(anno) %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(grepl('ERpHER2n', ClinGroup)) %>% dplyr::rename(sampleID=METABRIC_ID) %>% dplyr::select(sampleID,PAM50)
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
    anno <- anno %>% remove_rownames %>% column_to_rownames(var="sampleID")
    
    gex_data <- rownames_to_column(gex_data, "ensembl_gene_id") # name the column ensembl even though they arent so i can reuse the same functions as for the other cohorts
}

#######################################################################
# 3. Required data
#######################################################################

if (cohort!="Metabric") {
# 1. log transformed FPKM data
gex_data_log <- as.data.frame(log2(gex_data + 1))

# 2. sample annotation
sample_info <- anno[c("sampleID","PAM50")] %>% remove_rownames %>% column_to_rownames(var="sampleID")

# get the gex data in the correct format
gex_data_log <- rownames_to_column(gex_data_log, "ensembl_gene_id")
#gex_data_log$ensembl_gene_id <- gsub("\\..*","",gex_data_log$ensembl_gene_id) # commented this out because if not there are duplicate rows (when run for all of the genes, for the individual ones this has to be run) 
}

#######################################################################
# 4. ttest and gene score function definition
#######################################################################
# ttest function: should return 2 pvalues
sg_ttest <- function(id,gex_data_log,sample_info) {
    # sample ids
    Her2_cols <- sample_info %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="Her2") %>% pull(sampleID)
    LumA_cols <- sample_info %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumA") %>% pull(sampleID)
    LumB_cols <- sample_info %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumB") %>% pull(sampleID) 
    # extract the gex
    gex <- gex_data_log %>% filter(ensembl_gene_id == id) %>% column_to_rownames(var = "ensembl_gene_id")
    
    # vars
    hdata <- as.numeric(gex[1,Her2_cols])
    adata <- as.numeric(gex[1,LumA_cols])
    bdata <- as.numeric(gex[1,LumB_cols])
    
    # for Her2 vs LumA
    # equal variance check
    if (var.test(unlist(hdata),unlist(adata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = FALSE)
    } else {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = TRUE)
    }
    # for Her2 vs LumB
    # equal variance check
    if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
    } else {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
    }
    
    return(list("H2vsLA_pvalue" = H2vsLA_ttest_result$p.value,"H2vsLB_pvalue" = H2vsLB_ttest_result$p.value))
}

# gene score 
gene_score <- function(id,gex_data_log,method) {
    # extract the gex for each gene
    gex <- gex_data_log %>% filter(ensembl_gene_id == id) %>% column_to_rownames(var = "ensembl_gene_id")
    if (method=="scale") {
    # scale and calc. the score for each sample
    gex <- as.data.frame(t(scale(t(gex)))) 
    }
    result <- as.data.frame(apply(gex, 2, median)) %>% dplyr::rename(!!id:= "apply(gex, 2, median)")
    return(result)
}


#######################################################################
# 4. for selected single genes (ERBB2/GRB7/ESR1)
#######################################################################
# *<0.05, **<0.01 ***<0.001 ****< 0.0001   ns not significant

if (cohort!="Metabric") {
# save gex copy for the single gene analyses
sg_gex_data_log <- gex_data_log
# have to change the ids so i can match them to what i find in the databases
sg_gex_data_log$ensembl_gene_id <- gsub("\\..*","",gex_data_log$ensembl_gene_id)

# erbb2 #

erbb2_score <- gene_score("ENSG00000141736",sg_gex_data_log,method="scale")
erbb2_res <- sg_ttest("ENSG00000141736",sg_gex_data_log,sample_info)
erbb2_res$H2vsLA_pvalue # ns
erbb2_res$H2vsLB_pvalue # *

erbb2_score <- merge(erbb2_score,sample_info,by=0) %>% column_to_rownames(var = "Row.names")

#plot
ggplot(erbb2_score, aes(x=as.factor(PAM50),y=ENSG00000141736,fill=as.factor(PAM50))) +
    geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 subtype") +
    ylab("scaled log2 expression") +
    ylim(c(-5,10)) +
    ggtitle("ERBB2 expression (ENSG00000141736)")  +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="*", tip_length = 0.02, vjust=0.01, y_position = 8.8, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) + 
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")
#+geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4.5,nudge_x = 0.3,size=5) +

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ERBB2_expression.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm")


# ESR1 gene (as part of elucidating the steroid response metagene results) #

esr1_score <- gene_score("ENSG00000091831",sg_gex_data_log,method="scale")
#View(esr1_score)
esr1_score <- merge(esr1_score,sample_info,by=0) %>% column_to_rownames(var = "Row.names")
esr1_res <- sg_ttest("ENSG00000091831",sg_gex_data_log,sample_info)
esr1_res$H2vsLA_pvalue # ****
esr1_res$H2vsLB_pvalue # ****

#plot
ggplot(esr1_score, aes(x=as.factor(PAM50),y=ENSG00000091831,fill=as.factor(PAM50))) +
    geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 subtype") +
    ylab("scaled log2 expression") +
    ylim(c(-4,4.3))+
    ggtitle("ESR1 expression (ENSG00000091831)") +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="****", tip_length = 0.02, vjust=0.01, y_position = 3.5, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="****", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")
    #+ geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4,nudge_x = 0.3,size=5)

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ESR1_expression.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm")

} else if (cohort=="Metabric") {
    # erbb2 #
    # HGNC symbol: ERBB2 
    erbb2_score <- gene_score("ERBB2",gex_data,method="simple")
    erbb2_res <- sg_ttest("ERBB2",gex_data,anno)
    erbb2_res$H2vsLA_pvalue # ns
    erbb2_res$H2vsLB_pvalue # ****
    
    erbb2_score <- merge(erbb2_score,anno,by=0) %>% column_to_rownames(var = "Row.names")
    
    #plot
    ggplot(erbb2_score, aes(x=as.factor(PAM50),y=as.numeric(ERBB2),fill=as.factor(PAM50))) + 
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) + 
        xlab("PAM50 subtype") + 
        ylab("scaled log2 expression") + 
        ylim(c(-4,3.5)) + 
        ggtitle("ERBB2 expression")  + 
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations="****", tip_length = 0.02, vjust=0.01, y_position = 2.6, size = 2, textsize = 15) + 
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations="ns", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              legend.position = "none")
    #+geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4.5,nudge_x = 0.3,size=5) +
    
    ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ERBB2_expression.pdf", sep =""),
    width = 300,
    height = 300,
    units = "mm")
    
    # ESR1 gene (as part of elucidating the steroid response metagene results) #
    
    esr1_score <- gene_score("ESR1",gex_data,method="simple")
    #View(esr1_score)
    esr1_score <- merge(esr1_score,anno,by=0) %>% column_to_rownames(var = "Row.names")
    esr1_res <- sg_ttest("ESR1",gex_data,anno)
    esr1_res$H2vsLA_pvalue # ****
    esr1_res$H2vsLB_pvalue # ****
    
    #plot
    ggplot(esr1_score, aes(x=as.factor(PAM50),y=as.numeric(ESR1),fill=as.factor(PAM50))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 subtype") +
        ylab("scaled log2 expression") +
        ylim(c(-2.5,4.3))+
        ggtitle("ESR1 expression") +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations="****", tip_length = 0.02, vjust=0.01, y_position = 2.5, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations="****", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) +
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              legend.position="none")
    #+ geom_text(data=as.data.frame(dplyr::count(x=sample_info, group)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -4,nudge_x = 0.3,size=5)
    
    ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/",cohort,"/ESR1_expression.pdf", sep =""),
           width = 300,
           height = 300,
           units = "mm")
}







# #######################################################################
# # 5. for all genes (not finished)
# #######################################################################
# 
# # i have to rerun this because it crashed in the middle of it
# 
# for(i in 1:length(gex_data_log$ensembl_gene_id)) {
#     if(i==1) {
#         res <- gene_score(gex_data_log$ensembl_gene_id[i],gex_data_log)
#     } else {
#         temp_res <- gene_score(gex_data_log$ensembl_gene_id[i],gex_data_log)
#         res <- merge(res,temp_res,by=0) %>% column_to_rownames(var = "Row.names")
#     }
# }
# 
# View(res)
# 
# # add the pam50 annotations
# all_gene_scores <- merge(res,sample_info,by=0) %>% column_to_rownames(var = "Row.names")
# 
# # save the file
# save(all_gene_scores, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/all_gene_scores.RData", sep =""))
# 
# # create pdf with all the boxplots
# pdf("gene_plots.pdf", onefile = TRUE)
# for(i in 1:length(gex_data_log$ensembl_gene_id)-1) {
#     plot <- ggplot(all_gene_scores, aes(x=as.factor(group),y=all_gene_scores[,i])) +
#         geom_boxplot(fill="slateblue",alpha=0.2) +
#         xlab("PAM50 subtype") +
#         ylab("scaled log2 expression") +
#         ggtitle(colnames(all_gene_scores)[i])
#     print(plot)
# }
# dev.off()
# 
# 
# 
# ####
# 
# #reson <- some samples in anno but not in gex
# anno <- sample_info
# Her2_cols <- anno %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="Her2") %>% pull(sampleID)
# LumA_cols <- anno %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumA") %>% pull(sampleID)
# LumB_cols <- anno %>% rownames_to_column(var="sampleID") %>% filter(PAM50=="LumB") %>% pull(sampleID) 
# # extract the gex
# gex <- gex_data %>% filter(ensembl_gene_id == "ERBB2") %>% column_to_rownames(var = "ensembl_gene_id")
# 
# as.numeric(gex[1,Her2_cols])
# gex[1,LumA_cols]
# gex[1,LumB_cols]
# 
# # for Her2 vs LumA
# # equal variance check
# var.test(unlist(gex[1,Her2_cols]),unlist(gex[1,LumA_cols]), alternative = "two.sided")
# 
# H2vsLA_ttest_result <- t.test(as.numeric(gex[1,Her2_cols]),as.numeric(gex[1,LumA_cols]), var.equal = TRUE)
# t.test(gex[1,Her2_cols],gex[1,LumA_cols])
# 
# 
# 
# 
# #
# hdata <- as.numeric(gex[1,Her2_cols])
# adata <- as.numeric(gex[1,LumA_cols])
# bdata <- as.numeric(gex[1,LumB_cols])
# 
# var.test(unlist(hdata),unlist(adata), alternative = "two.sided")
# 
# H2vsLA_ttest_result <- t.test(as.numeric(gex[1,Her2_cols]),as.numeric(gex[1,LumA_cols]), var.equal = TRUE)
# t.test(gex[1,Her2_cols],gex[1,LumA_cols])