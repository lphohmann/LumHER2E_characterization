# heatmap script for the core set of DEGs using SCANB and Metabric gex data

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
library(janitor)

#######################################################################
# 1. set parameters
#######################################################################

# 1.1 input for which cohort the analysis should be run
cohort <- "Metabric" # Metabric or SCANB or TCGA

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
    
    #View(gex_data_log)
    
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
    
    # convert hugo symbols to ensembl
    gex_data <- gex_data %>% rownames_to_column(var="Hugo_Symbol") %>% filter(Hugo_Symbol %in% gene_anno$Hugo_Symbol)
    gex_data <- as.data.frame(merge(gex_data,gene_anno,by="Hugo_Symbol")) %>% filter(duplicated(ensembl_gene_id) == FALSE) %>%  column_to_rownames(var="ensembl_gene_id") %>% dplyr::select(-c("Hugo_Symbol","Entrez_Gene_Id"))
    
    gex_data_log <- gex_data
}


# 2. sample annotation
sample_info <- anno[c("sampleID","PAM50")]
#save(sample_info, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",output_dir,"/DE_anno.RData", sep = ""))

#######################################################################
# 4. Load DEGs of different cohorts and define the core set
#######################################################################

# load scan-b fdr DEGs
load(paste("~/Desktop/MTP_project/Output/Transcriptomics/","SCANB","/complete_DE_results.RData",sep = ""))
scanb_DEGs <- results
scanb_DEGs <- scanb_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,H2vsLB_padj_fdr = H2vsLB_padj) %>% rownames_to_column("Gene")

# load metabric fdr DEGs
load(paste("~/Desktop/MTP_project/Output/Transcriptomics/","Metabric","/complete_DE_results.RData",sep = ""))
mb_DEGs <- results
mb_DEGs <- mb_DEGs %>% dplyr::rename(H2vsLA_padj_fdr = H2vsLA_padj,H2vsLB_padj_fdr = H2vsLB_padj) %>% rownames_to_column("Gene")

# have to make the identifiers uniform so they can be compred between cohorts
# for metabric get the identifiers conform with scan be ids
load(file = paste("~/Desktop/MTP_project/Output/Transcriptomics/","Metabric","/gene_anno.RData",sep = ""))
gene_anno <- gene_anno %>% dplyr::rename(Gene = Hugo_Symbol) %>% dplyr::select(Gene,ensembl_gene_id)
mb_DEGs <- mb_DEGs %>% filter(Gene %in% gene_anno$Gene) 
mb_DEGs <- merge(mb_DEGs,gene_anno,by="Gene") 
mb_DEGs <- as.data.frame(mb_DEGs) %>% dplyr::select(-c(Gene)) %>% dplyr::rename(Gene = ensembl_gene_id)

scanb_DEGs$Gene <- gsub("\\..*","",scanb_DEGs$Gene)




# DEGs
H2vsLA_DEGs_mb_fdr <- mb_DEGs %>% dplyr::filter(H2vsLA_padj_fdr <= 0.01) #%>% dplyr::pull(Gene)
H2vsLB_DEGs_mb_fdr <- mb_DEGs %>% dplyr::filter(H2vsLB_padj_fdr <= 0.01) #%>% dplyr::pull(Gene)
H2vsLA_DEGs_scanb_fdr <- scanb_DEGs %>% dplyr::filter(H2vsLA_padj_fdr <= 0.01) #%>% dplyr::pull(Gene)
H2vsLB_DEGs_scanb_fdr <- scanb_DEGs %>% dplyr::filter(H2vsLB_padj_fdr <= 0.01) #%>% dplyr::pull(Gene)

length(H2vsLA_DEGs_mb_fdr$Gene)/17780*100
length(H2vsLB_DEGs_mb_fdr$Gene)/17780*100
length(H2vsLA_DEGs_scanb_fdr$Gene)
length(H2vsLB_DEGs_scanb_fdr$Gene)

# core genes
H2vsLA_core_DEGs <- H2vsLA_DEGs_scanb_fdr %>% filter(Gene %in% H2vsLA_DEGs_mb_fdr$Gene) #%>% pull(Gene)
H2vsLB_core_DEGs <- H2vsLB_DEGs_scanb_fdr %>% filter(Gene %in% H2vsLB_DEGs_mb_fdr$Gene) #%>% pull(Gene)
length(H2vsLA_core_DEGs$Gene)
length(H2vsLB_core_DEGs$Gene)

#######################################################################
# 4. Export the core set for functional enrichment analysis
#######################################################################

# export for functional enrichment analysis
#write.table(H2vsLA_core_DEGs$Gene, file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/H2vsLA_DEGs.txt",sep = ""), sep="\t", row.names = FALSE, quote = FALSE)
#write.table(H2vsLB_core_DEGs$Gene, file = paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/H2vsLB_DEGs.txt",sep = ""), sep="\t", row.names = FALSE, quote = FALSE)

#######################################################################
# 5. Heatmap parameters and functions
#######################################################################

n.tree <- 10  # nbr of branches/clusters -> based on what do I choose this? 3 for the pam50 subtypes?
nbr.c <- 2 # color palette size = max value (max value in top_gex is this)
breaksList = seq(-2, nbr.c, by = 0.1) # this defines the color palette range. 
my.link.method <- "ward.D" #
my.dist.method <- "euclidean" #"correlation" "euclidean"
my.meth <- paste("heatmap: COLS=",my.dist.method,"+",my.link.method,", ROWS=",my.dist.method,"+",my.link.method)

# colors
my_colour = list(Cluster=c("1"="#E41A1C","2"="#377EB8","3"="#4DAF4A"),
                 PAM50 = c("Her2" = "red" , "LumA" = "blue" , "LumB" = "yellow"),
                 NHG = c("missing"="white","1"="#E41A1C","2"="lightblue","3"="#4DAF4A"), # some no NHG (add stht o make white?)
                 Basal = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 Early_response = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 IR = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 Lipid = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 Mitotic_checkpoint = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 Mitotic_progression = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 SR = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"),
                 Stroma = c("<= -2"="#2e4053","-1 to -2"="#5d6d7e","-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1","0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00"))

#######

# bin continuous anno variables to color them correctly (e.g. the metagene scores)
mg_converter <- function(mg) {
    mg.class<-rep(NA,nrow(hm_anno))
    mg.class[which(mg <= -2)] <- "<= -2"
    mg.class[which((mg <= -1) & (mg > -2) )] <- "-1 to -2"
    mg.class[which((mg <= -0.5) & (mg > -1) )] <- "-0.5 to -1"
    mg.class[which((mg >= -0.5) & (mg < 0.5) )] <- "-0.5 to 0.5"
    mg.class[which((mg >= 0.5) & (mg < 1) )] <- "0.5 to 1"
    mg.class[which((mg >= 1) & (mg < 2) )] <- "1 to 2"
    mg.class[which(mg >= 2)] <- ">= 2"
    return(mg.class)
}

#######################################################################
# 6. Heatmap 1: Her2 vs LumA
#######################################################################

#### Part 1: select the gene subset and annotations to visualize ####

# 1. get the gex of only the relevant genes
top_gex <- gex_data_log %>% rownames_to_column("Gene") 
top_gex$Gene <- gsub("\\..*","",top_gex$Gene)
top_gex <- top_gex[!duplicated(top_gex$Gene), ] # remove duplicates

top_gex <- top_gex %>% filter(Gene %in% H2vsLA_core_DEGs$Gene) %>% column_to_rownames(var="Gene") 

# 2. annotations - include (pam50, grade, metagenes)
hm_anno <- anno[c("sampleID","PAM50","NHG")] %>% remove_rownames %>% column_to_rownames(var="sampleID") %>% dplyr::rename(PAM50_subtype=PAM50)

# 3. add the metagene scores (idea: add erbb2 also)
load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))
scaled_mg_scores$group <- NULL

# 4. handle missing NHG values
hm_anno$NHG <- as.character(hm_anno$NHG)
hm_anno$NHG[hm_anno$NHG == ""] <- "missing"

# 5. heatmap annotation data
scaled_mg_scores$NHG <- NULL
hm_anno <- merge(hm_anno,scaled_mg_scores,by=0) %>% column_to_rownames(var = "Row.names")

#### Part 2: create the heatmap ####

#install.packages("amap")
library(amap)
# cluster samples # old hm with dist() and euclidean not i try new one
c1 <- cutree(hclust(Dist(t(top_gex),method=my.dist.method),method=my.link.method),n.tree)

# merge
hm_anno <- merge(hm_anno,as.data.frame(c1),by=0) %>% column_to_rownames(var = "Row.names") %>% dplyr::rename(Cluster=c1)

#View(hm_anno)

# final hm annotation file
final_hm_anno <- data.frame(Cluster=factor(hm_anno$Cluster),
                            PAM50=factor(hm_anno$PAM50_subtype),
                            NHG=factor(hm_anno$NHG),
                            Basal=factor(mg_converter(hm_anno$Basal)),
                            IR=factor(mg_converter(hm_anno$IR)),
                            Lipid=factor(mg_converter(hm_anno$Lipid)),
                            Mitotic_checkpoint=factor(mg_converter(hm_anno$Mitotic_checkpoint)),
                            SR=factor(mg_converter(hm_anno$SR)),
                            Stroma=factor(mg_converter(hm_anno$Stroma)))
rownames(final_hm_anno) <- rownames(hm_anno)

#Early_response=factor(mg_converter(hm_anno$Early_response)),
#Mitotic_progression=factor(mg_converter(hm_anno$Mitotic_progression)),

#save(final_hm_anno, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_anno.RData",sep = ""))
#save(top_gex, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_top_gex.RData",sep = ""))
#load(file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_anno.RData",sep = ""))
#load(file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_top_gex.RData",sep = ""))

# lines if i want to do it for a subset
# final_hm_anno <- final_hm_anno %>% filter(PAM50 %in% c("Her2","LumB"))
# top_gex <- top_gex[,which(colnames(top_gex) %in% rownames(final_hm_anno))]

# create the heatmap
pheatmap(top_gex, cluster_rows=T, clustering_distance_rows = my.dist.method,
         clustering_distance_cols = my.dist.method, clustering_method = my.link.method,
         show_rownames=F, show_colnames=F, main=paste("Gene expression heatmap of core Her2/LumA DEGs in ",cohort,sep = ""), cutree_cols=n.tree,
         legend = T,
         color=colorRampPalette(c("blue", "white", "yellow"))(length(breaksList)),
         annotation_col=final_hm_anno[,c("Stroma","SR","Mitotic_checkpoint","Lipid","IR","Basal","NHG","PAM50")], annotation_colors=my_colour,breaks=breaksList,
         filename= paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/H2vsLA_core_DEGs_heatmap_",cohort,"_mbcore_001pval.pdf",sep = "")) # change name depending on scanb or metabric DELETE

#annotation_col=final_hm_anno[,c("Stroma","SR","Mitotic_progression","Mitotic_checkpoint","Lipid","IR","Early_response","Basal","NHG","PAM50")],

# try to make the heapmap nicer (only show annotation for one metagene) i.e. delete the other annotation legends 
#library(grid)
#grid.ls(grid.force())

#######################################################################
# 6. Heatmap 2: Her2 vs LumB
#######################################################################

#### Part 1: select the gene subset and annotations to visualize ####

# 1. get the gex of only the relevant genes
top_gex <- gex_data_log %>% rownames_to_column("Gene") 
top_gex$Gene <- gsub("\\..*","",top_gex$Gene)
top_gex <- top_gex[!duplicated(top_gex$Gene), ] # remove duplicates
top_gex <- top_gex %>% filter(Gene %in% H2vsLB_core_DEGs$Gene) %>% column_to_rownames(var="Gene") 

# 2. annotations - include (pam50, grade, metagenes)
hm_anno <- anno[c("sampleID","PAM50","NHG")] %>% remove_rownames %>% column_to_rownames(var="sampleID") %>% dplyr::rename(PAM50_subtype=PAM50)

# 3. add the metagene scores (idea: add erbb2 also)
load(paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/sample_metagene_scores.RData",sep = ""))
scaled_mg_scores$group <- NULL

# 4. handle missing NHG values
hm_anno$NHG <- as.character(hm_anno$NHG)
hm_anno$NHG[hm_anno$NHG == ""] <- "missing"

# 5. heatmap annotation data
scaled_mg_scores$NHG <- NULL
hm_anno <- merge(hm_anno,scaled_mg_scores,by=0) %>% column_to_rownames(var = "Row.names")

#### Part 2: create the heatmap ####

# cluster samples
c1 <- cutree(hclust(Dist(t(top_gex),method=my.dist.method),method=my.link.method),n.tree)

# merge the cluster anno
#hm_anno <- hm_anno %>% dplyr::select(-c("Cluster"))
hm_anno <- merge(hm_anno,as.data.frame(c1),by=0) %>% column_to_rownames(var = "Row.names") %>% dplyr::rename(Cluster=c1)

# final hm annotation file
final_hm_anno <- data.frame(Cluster=factor(hm_anno$Cluster),
                            PAM50=factor(hm_anno$PAM50_subtype),
                            NHG=factor(hm_anno$NHG),
                            Basal=factor(mg_converter(hm_anno$Basal)),
                            IR=factor(mg_converter(hm_anno$IR)),
                            Lipid=factor(mg_converter(hm_anno$Lipid)),
                            Mitotic_checkpoint=factor(mg_converter(hm_anno$Mitotic_checkpoint)),
                            SR=factor(mg_converter(hm_anno$SR)),
                            Stroma=factor(mg_converter(hm_anno$Stroma)))
rownames(final_hm_anno)<-rownames(hm_anno)

#Early_response=factor(mg_converter(hm_anno$Early_response)),
#Mitotic_progression=factor(mg_converter(hm_anno$Mitotic_progression)),

#save(final_hm_anno, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_anno.RData",sep = ""))
#save(top_gex, file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_top_gex.RData",sep = ""))
#load(file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_anno.RData",sep = ""))
#load(file = paste("~/Desktop/MTP_project/Output/Transcriptomics/",cohort,"/hm_top_gex.RData",sep = ""))



# # lines if i want to do it for a subset
# final_hm_anno <- final_hm_anno %>% filter(PAM50 %in% c("Her2","LumB"))
# top_gex <- top_gex[,which(colnames(top_gex) %in% rownames(final_hm_anno))]

# create the heatmap
pheatmap(top_gex, cluster_rows=T, clustering_distance_rows = my.dist.method,
         clustering_distance_cols = my.dist.method, clustering_method = my.link.method,
         show_rownames=F, show_colnames=F, main=paste("Gene expression heatmap of core Her2/LumB DEGs in ",cohort,sep = ""), cutree_cols=n.tree,
         legend = T,
         color=colorRampPalette(c("blue", "white", "yellow"))(length(breaksList)),
         annotation_col=final_hm_anno[,c("Stroma","SR","Mitotic_checkpoint","Lipid","IR","Basal","NHG","PAM50")], annotation_colors=my_colour,breaks=breaksList,
         filename= paste("~/Desktop/MTP_project/Output/Plots/Transcriptomics/intercohort/H2vsLB_core_DEGs_heatmap_",cohort,"_mbcore_001pval.pdf",sep = ""))

# try to make the heapmap nicer (only show annotation for one metagene) i.e. delete the other annotation legends 
#library(grid)
#grid.ls(grid.force())

