# Script: Copy number alteration analyses in the Metabric cohort

# TODO: 
# - check cytoband positions to see if they are co-localized

# empty environment
rm(list=ls())

# set working directory to where the data is
setwd("~/Desktop/MTP_project")

#packages
library(ggplot2)
library(tidyverse)
library(biomaRt)

library(regioneR)
library(karyoploteR)
library(CopyNumberPlots)
library(BSgenome.Hsapiens.UCSC.hg19)
#BiocManager::install("plyranges")
library(plyranges)

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################
cohort <- "Metabric"
# 2.1 load data
cn.data <- read.delim("Data/Metabric/data_CNA.txt",header = TRUE, sep = "\t", dec = ".")
#cn.meta <- read.delim("Data/Metabric/meta_CNA.txt",header = TRUE, sep = "\t", dec = ".")
# annotation data
load("Data/Metabric/Annotations/Merged_annotations.RData")
anno <- as.data.frame(anno) %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) %>% filter(grepl('ERpHER2n', ClinGroup)) %>% dplyr::rename(sampleID=METABRIC_ID) %>% dplyr::select(sampleID,PAM50)
anno$sampleID <- gsub("-",".",anno$sampleID) # have to get identifiers in same format
cn.data$Entrez_Gene_Id <- NULL
cn.data <- cn.data %>% column_to_rownames(var = "Hugo_Symbol") 
cn.data <- cn.data[,colnames(cn.data) %in% anno$sampleID]
anno <- anno %>% filter(sampleID %in% colnames(cn.data))

##########################################################################
# Part 1: general gain/loss, whole genome, freq. altered per group boxplot
##########################################################################

# alterations per sample 
res <- colSums(cn.data != 0) %>% as.data.frame() %>% rownames_to_column(var = "sampleID") %>% dplyr::rename(Nalt = 2)
#as percent
res$Nalt <- (res$Nalt/nrow(cn.data))*100

# add to anno
p1.anno <- merge(anno,res,by="sampleID")
p1.anno <- p1.anno %>% drop_na(Nalt)

# calc. significance
# luma
library(car)
library(ggsignif)
p1.anno$PAM50 <- as.factor(p1.anno$PAM50)
adata <- p1.anno %>% filter(PAM50 %in% c("LumA", "Her2"))
bdata <- p1.anno %>% filter(PAM50 %in% c("LumB", "Her2"))
# check equal variance
leveneTest(Nalt~PAM50, data=adata) # signif.
round(t.test(Nalt~PAM50, data=adata, var.equal = FALSE)$p.value,7) # 0.000001 -> signif. 

# lumb
leveneTest(Nalt~PAM50, data=bdata) # not signif.
round(t.test(Nalt~PAM50, data=bdata, var.equal = TRUE)$p.value,7) # 0.17-> not signif. 

# plot
#pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/overall_CN_alterations.pdf", sep =""))

plot <- ggplot(p1.anno, aes(x=as.factor(PAM50),y=Nalt,fill=PAM50)) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("PAM50 subtype") +
        ylab("Copy-number alterations") +
        scale_y_continuous(limits=c(0,28000),breaks = c(0,5000,10000,15000,20000,25000)) +
        ggtitle("Genome-wide copy-number alterations per luminal PAM50 subtype") +
        geom_text(data=as.data.frame(dplyr::count(x=p1.anno, PAM50)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -2,nudge_x = 0.3,size=5) +
        geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 26000, size = 2, textsize = 15) +
        geom_signif(comparisons=list(c("Her2", "LumA")), annotations="****", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) + 
        theme(axis.text.x = element_text(size = 30),
              axis.title.x = element_text(size = 35),
              axis.text.y = element_text(size = 30),
              axis.title.y = element_text(size = 35),
              legend.position = "none")
print(plot)
#dev.off()


# same plot but as fraction
#pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/perc_overall_CN_alterations.pdf", sep =""))
plot <- ggplot(p1.anno, aes(x=as.factor(PAM50),y=Nalt,fill=PAM50)) +
    geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 subtype") +
    ylab("Genome altered (%)") +
    scale_y_continuous(limits=c(0,0.0055),breaks = seq(0,0.005,0.001)) + #ggtitle("Genome-wide copy-number alterations per luminal PAM50 subtype") +
    geom_text(data=as.data.frame(dplyr::count(x=p1.anno, PAM50)), aes(y = 0, label = paste("n=",n,sep = "")),nudge_y = -2,nudge_x = 0.3,size=5) +
    geom_signif(comparisons=list(c("Her2", "LumB")), annotations="ns", tip_length = 0.02, vjust=0.01, y_position = 0.005, size = 2, textsize = 15) +
    geom_signif(comparisons=list(c("Her2", "LumA")), annotations="****", tip_length = 0.02, vjust=0.01, size = 2, textsize = 15) + 
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")
print(plot)
#dev.off()

#######################################################################
# Part 2: spatial gain/loss - get chromosomal position and cytoband from gene symbols 
#######################################################################

# apprach 1: - Biomart & Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://grch37.ensembl.org") # GRCh37

#listEnsemblArchives()

# Only use standard human chromosomes
normal.chroms <- c(1:22, "X", "Y", "M")

# Filter on HGNC symbol and chromosome, retrieve genomic location and band
my.symbols <- row.names(cn.data)

my.regions <- getBM(c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                    filters = c("hgnc_symbol", "chromosome_name"),
                    values = list(hgnc_symbol=my.symbols, chromosome_name=normal.chroms),
                    mart = ensembl)

#save(my.regions, file = paste("~/Desktop/MTP_project/Output/CopyNumber/",cohort,"/gene_regions.RData",sep = ""))
load(file = paste("~/Desktop/MTP_project/Output/CopyNumber/",cohort,"/gene_regions.RData",sep = ""))

#######################################################################
# Part 2: spatial gain/loss - plot each CNV of each sample 
#######################################################################
# manual approach
# get data into the right format
cn.data.meta <- as.data.frame(merge(my.regions,rownames_to_column(cn.data,var = "hgnc_symbol"), by="hgnc_symbol")) %>% dplyr::rename(chr=chromosome_name,start=start_position,end=end_position) %>% mutate(pos=end-((end-start)/2)) %>% relocate(pos, .after=band) %>% select(-c("start","end"))
# rename X chromosome to 23
cn.data.meta$chr[cn.data.meta$chr == "X"] <- "23" 
cn.data.meta$chr <- as.numeric(cn.data.meta$chr)
# sort
cn.data.meta <- cn.data.meta[with(cn.data.meta, order(chr, pos)),] 
# get chromosome lengths (i simply take the last pos on each)
chr.lengths <- cn.data.meta %>% group_by(chr) %>% summarise(length = max(pos)) %>% as.data.frame()
chr.lengths <- chr.lengths[with(chr.lengths, order(chr)),] %>% mutate(genome = cumsum(length))#%>% dplyr::rename(chromosome=chr)
# add genome pos 
cn.data.meta <- cn.data.meta %>% rowwise() %>% 
    mutate(genome = ifelse(chr==1,pos,pos + sum(chr.lengths$length[1:chr-1]))) %>% 
    relocate(genome, .after=pos)

# plot

pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/m_sample_CN_profiles.pdf", sep =""),onefile = TRUE, height = 21.0, width = 29.7)

pb = txtProgressBar(min = 0, max = ncol(cn.data.meta)-5, initial = 0, style = 3) 
# for each sample
for(i in c(1:(ncol(cn.data.meta)-5))) { 
    setTxtProgressBar(pb,i)
    # loop iteration example
    sampleID <- colnames(cn.data.meta)[i+5]
    samplePAM50 <- anno[which(anno$sampleID==sampleID),]$PAM50
    sample.data <- as.data.frame(cn.data.meta[,c(1:5,i+5)]) 
    
    plot <- ggplot(sample.data, aes(x=sample.data$genome, y=sample.data[,sampleID])) + #, group=1
        scale_y_continuous(limits=c(-2, 2), breaks=c(-2, -1, 0, 1, 2)) +
        geom_step() + #geom_path()
        geom_point() +
        geom_vline(xintercept = chr.lengths$genome, linetype="dotted") +
        ylab("gain/loss") +
        xlab("genome position (chromosome)") +
        scale_x_continuous(breaks=chr.lengths$genome,
                           labels=chr.lengths$chr)
    print(plot)
}
dev.off()


##################################################################
# Part 2: spatial gain/loss - plot CNV of each group
##################################################################

# for each group calc. the frequency of gain and loss for each gene
# sampleIDs
HER2.samples <- anno[which(anno$PAM50 == "Her2"),]$sampleID
LUMA.samples <- anno[which(anno$PAM50 == "LumA"),]$sampleID
LUMB.samples <- anno[which(anno$PAM50 == "LumB"),]$sampleID

# calc loss/gain freqs per group
cn.data.meta$freqloss.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x<0))/length(HER2.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x<0))/length(LUMA.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x<0))/length(LUMB.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqgain.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x>0))/length(HER2.samples))*100)
cn.data.meta$freqgain.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x>0))/length(LUMA.samples))*100)
cn.data.meta$freqgain.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x>0))/length(LUMB.samples))*100)

# plot

#pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/PAM50_CN_profiles.pdf", sep =""),onefile = TRUE, height = 21.0, width = 29.7)

#her2
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations in luminal HER2-enriched BC") +
    geom_line(aes(y = freqgain.H2), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.H2), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=c(1:22,"X"), #chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100)) +
    theme(text=element_text(size=27)) #+
    # annotate(x=80000, y=c(-25,25), label=c("Loss","Gain"), 
    #          geom="text", angle=90, hjust=0.5, size=10, colour=c("darkred","darkgreen")) 
    
print(plot)


#luma
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations in luminal A BC") +
    geom_line(aes(y = freqgain.LUMA), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.LUMA), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100)) +
    theme(text=element_text(size=27)) #+
    # annotate(x=95000, y=c(-25,25), label=c("Loss","Gain"), 
    #          geom="text", angle=90, hjust=0.5, size=10, colour=c("darkred","darkgreen")) 
    
print(plot)

#luma
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations in luminal B BC") +
    geom_line(aes(y = freqgain.LUMB), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.LUMB), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100)) +
    theme(text=element_text(size=27)) #+
    # annotate(x=-10000, y=c(-25,25), label=c("Loss","Gain"), 
    #          geom="text", angle=90, hjust=0.5, size=10, colour=c("darkred","darkgreen"))

print(plot)

# combined
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations in luminal BC") +
    geom_line(aes(y = freqgain.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqloss.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqgain.LUMA, color = "green"),size=4) + 
    geom_line(aes(y = freqloss.LUMA, color = "green"),size=4) + 
    geom_line(aes(y = freqgain.LUMB, color = "blue"),size=4) + 
    geom_line(aes(y = freqloss.LUMB, color = "blue"),size=4) + 
    scale_colour_manual(name="PAM50", values = c("#d72e2b", "#07a109", "#2176d5"),
                        labels = c("HER2", "LUMA", "LUMB")) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
    scale_x_continuous(name="Genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=c(as.character(1:22),"X"),
                       limits = c(0,max(cn.data.meta$genome)+50000000),
                       expand = c(0, 0)) +
    scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100),
                       expand = c(0, 0),
                       limits = c(-80,80)) +
    theme(text=element_text(size=30),
          legend.title = element_blank(),
          axis.title.y = element_text(vjust = 0.5),
          legend.position = c(0.97, 0.90)) +
    annotate(x=45000000,y=c(-50,50), label=c("Loss","Gain"), 
    geom="text", angle=90, hjust=0.5, size=9, colour=c("black","black")) 


# axis.title.y = element_text(angle=0, vjust = 0.5)

print(plot)

ggsave(filename="PAM50_CN_profiles.pdf", #_basal
       width = 720,
       height = 210,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/", sep =""))

#dev.off()

##################################################################
# Part 2: spatial gain/loss - calc. singificance of CN alterations between groups
##################################################################

# data

# genes
genes <- cn.data.meta$hgnc_symbol

# objects to store results
res.matrix <- data.frame(matrix(ncol = 9, nrow = length(cn.data.meta$hgnc_symbol))) %>% dplyr::rename(Gene=1, Loss_LUMA_pval=2, Loss_LUMB_pval=3 , Gain_LUMA_pval=4, Gain_LUMB_pval=5,DLoss_LUMA_pval=6, DLoss_LUMB_pval=7 , DGain_LUMA_pval=8, DGain_LUMB_pval=9)

# main loop testing all genes for gain/loss normal and driver
pb = txtProgressBar(min = 0, max = length(cn.data.meta$hgnc_symbol), initial = 0, style = 3) 
for (i in 1:length(cn.data.meta$hgnc_symbol)) {
    setTxtProgressBar(pb,i)
    res.matrix$Gene[i] <- cn.data.meta$hgnc_symbol[i]
    # transpose data # samples geneval pam50 columns
    cn.data <- cn.data.meta[i,] %>% column_to_rownames(var="hgnc_symbol") %>% select(-c("chr","band","pos","genome"))
    cn.data <- as.data.frame(t(cn.data)) %>% rownames_to_column(var="sampleID") 
    cn.data <- merge(cn.data,anno,by="sampleID") #add pam50anno as column
    
    # group data
    LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
    LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
    
    # normal
    # 1. gain
    LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
    LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
    LUMA.data[[2]][LUMA.data[[2]] <= 0] <- 0
    LUMA.data[[2]][LUMA.data[[2]] >= 1] <- 1
    LUMB.data[[2]][LUMB.data[[2]] <= 0] <- 0
    LUMB.data[[2]][LUMB.data[[2]] >= 1] <- 1
    
    # make count tables
    # lumA
    LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
    if(nrow(LUMA.cont.table)==2) {
        res.matrix$Gain_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
    } else {res.matrix$Gain_LUMA_pval[i] <- 1}
    # lumB
    LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
    if(nrow(LUMB.cont.table)==2) {
        res.matrix$Gain_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
    } else {res.matrix$Gain_LUMB_pval[i] <- 1}
    # 2. loss
    LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
    LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
    LUMA.data[[2]][LUMA.data[[2]] <= -1] <- -1
    LUMA.data[[2]][LUMA.data[[2]] >= 0] <- 0
    LUMB.data[[2]][LUMB.data[[2]] <= -1] <- -1
    LUMB.data[[2]][LUMB.data[[2]] >= 0] <- 0
    # lumA
    LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
    if(nrow(LUMA.cont.table)==2) {
        res.matrix$Loss_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
    } else {res.matrix$Loss_LUMA_pval[i] <- 1}
    # lumB
    LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
    if(nrow(LUMB.cont.table)==2) {
        res.matrix$Loss_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
    } else {res.matrix$Loss_LUMB_pval[i] <- 1}
    
    # 3. driver gain
    LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
    LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
    LUMA.data[[2]][LUMA.data[[2]] <= 1] <- 0
    LUMA.data[[2]][LUMA.data[[2]] >= 2] <- 2
    LUMB.data[[2]][LUMB.data[[2]] <= 1] <- 0
    LUMB.data[[2]][LUMB.data[[2]] >= 2] <- 2
    
    # make count tables
    # lumA
    LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
    if(nrow(LUMA.cont.table)==2) {
        res.matrix$DGain_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
    } else {res.matrix$DGain_LUMA_pval[i] <- 1}
    # lumB
    LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
    if(nrow(LUMB.cont.table)==2) {
        res.matrix$DGain_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
    } else {res.matrix$DGain_LUMB_pval[i] <- 1}
    # 4. driver loss
    LUMA.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumA"))
    LUMB.data <- cn.data %>% filter(PAM50 %in% c("Her2","LumB"))
    LUMA.data[[2]][LUMA.data[[2]] <= -2] <- -2
    LUMA.data[[2]][LUMA.data[[2]] >= -1] <- 0
    LUMB.data[[2]][LUMB.data[[2]] <= -2] <- -2
    LUMB.data[[2]][LUMB.data[[2]] >= -1] <- 0
    
    # make count tables
    # lumA
    LUMA.cont.table <- table(LUMA.data[[cn.data.meta$hgnc_symbol[i]]],LUMA.data$PAM50)
    if(nrow(LUMA.cont.table)==2) {
        res.matrix$DLoss_LUMA_pval[i] <- fisher.test(LUMA.cont.table)$p.value
    } else {res.matrix$DLoss_LUMA_pval[i] <- 1}
    # lumB
    LUMB.cont.table <- table(LUMB.data[[cn.data.meta$hgnc_symbol[i]]],LUMB.data$PAM50)
    if(nrow(LUMB.cont.table)==2) {
        res.matrix$DLoss_LUMB_pval[i] <- fisher.test(LUMB.cont.table)$p.value
    } else {res.matrix$DLoss_LUMB_pval[i] <- 1}
}

View(res.matrix)

# significant genes
# adjust
res.matrix$Loss_LUMA_pval <- p.adjust(res.matrix$Loss_LUMA_pval,method = "fdr")
res.matrix$Gain_LUMA_pval <- p.adjust(res.matrix$Gain_LUMA_pval,method = "fdr")
res.matrix$Loss_LUMB_pval <- p.adjust(res.matrix$Loss_LUMB_pval,method = "fdr")
res.matrix$Gain_LUMB_pval <- p.adjust(res.matrix$Gain_LUMB_pval,method = "fdr")
res.matrix$DLoss_LUMA_pval <- p.adjust(res.matrix$DLoss_LUMA_pval,method = "fdr")
res.matrix$DGain_LUMA_pval <- p.adjust(res.matrix$DGain_LUMA_pval,method = "fdr")
res.matrix$DLoss_LUMB_pval <- p.adjust(res.matrix$DLoss_LUMB_pval,method = "fdr")
res.matrix$DGain_LUMB_pval <- p.adjust(res.matrix$DGain_LUMB_pval,method = "fdr")

# save
#save(res.matrix, file = paste("~/Desktop/MTP_project/Output/CopyNumber/",cohort,"/cn_enrichment.RData",sep = ""))

load(paste("~/Desktop/MTP_project/Output/CopyNumber/",cohort,"/cn_enrichment.RData",sep = ""))

#modify my.regions to include chromosome name in the cytoband column
my.regions <- my.regions %>% mutate(band = paste(chromosome_name,band,sep=""))
my.regions <- dplyr::rename(my.regions,Gene=hgnc_symbol)
res.matrix <- merge(res.matrix,my.regions,by="Gene") %>% as.data.frame()

# gain
Gain_LUMA <- res.matrix %>% filter(Gain_LUMA_pval  <= 0.05) %>% dplyr::select(c(Gene,band))
Gain_LUMA_CB <- unique(Gain_LUMA$band)
Gain_LUMB <- res.matrix %>% filter(Gain_LUMB_pval <= 0.05) %>% select(c(Gene,band))
Gain_LUMB_CB <- unique(Gain_LUMB$band)
# loss
Loss_LUMA <- res.matrix %>% filter(Loss_LUMA_pval  <= 0.05) %>% select(c(Gene,band))
Loss_LUMA_CB <- unique(Loss_LUMA$band)
Loss_LUMB <- res.matrix %>% filter(Loss_LUMB_pval  <= 0.05) %>% select(c(Gene,band))
Loss_LUMB_CB <- unique(Loss_LUMB$band)
# gain
DGain_LUMA <- res.matrix %>% filter(DGain_LUMA_pval <= 0.05) %>% select(c(Gene,band))
DGain_LUMA_CB <- unique(DGain_LUMA$band)
DGain_LUMB <- res.matrix %>% filter(DGain_LUMB_pval  <= 0.05) %>% select(c(Gene,band))
DGain_LUMB_CB <- unique(DGain_LUMB$band)
# loss
DLoss_LUMA <- res.matrix %>% filter(DLoss_LUMA_pval <= 0.05) %>% select(c(Gene,band))
DLoss_LUMA_CB <- unique(DLoss_LUMA$band)
DLoss_LUMB <- res.matrix %>% filter(DLoss_LUMB_pval <= 0.05) %>% select(c(Gene,band))
DLoss_LUMB_CB <- unique(DLoss_LUMB$band)

length(Gain_LUMA$Gene) #1385
length(Gain_LUMB$Gene) #0
length(Loss_LUMA$Gene) #15757
length(Loss_LUMB$Gene) #133
length(DGain_LUMA$Gene) #1240
length(DGain_LUMB$Gene) #0
length(DLoss_LUMA$Gene) #0
length(DLoss_LUMB$Gene) #0

sort(Gain_LUMA_CB,decreasing = TRUE) # 8q21-24;6q21-25; etc.
sort(Gain_LUMB_CB,decreasing = TRUE) #0
sort(Loss_LUMA_CB,decreasing = TRUE) # big list
sort(Loss_LUMB_CB,decreasing = TRUE) # 5q21-23; 12q13
sort(DGain_LUMA_CB,decreasing = TRUE) #8,19,17,12
sort(DGain_LUMB_CB,decreasing = TRUE) #0
sort(DLoss_LUMA_CB,decreasing = TRUE)
sort(DLoss_LUMB_CB,decreasing = TRUE)


##########################################################################
# barplots to visualize the reustls furthetr
##########################################################################

length(Gain_LUMA$Gene) #1385
length(Gain_LUMB$Gene) #0
length(Loss_LUMA$Gene) #15757
length(Loss_LUMB$Gene) #133
# Gains vs Losses
plot.data <- data.frame(matrix(ncol = 3, nrow = 4))
colnames(plot.data) <- c("Cat","pam50","N")
plot.data$Cat <- c("Loss","Loss","Gain","Gain")
plot.data$pam50 <- c("LumA","LumB","LumA","LumB")
plot.data$N <- c(length(Loss_LUMA$Gene),length(Loss_LUMB$Gene),length(Gain_LUMA$Gene),length(Gain_LUMB$Gene))

# plot losses 
ggplot(plot.data[which(plot.data$Cat=="Loss"),], aes(x=as.factor(pam50),y=N,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Significant loss alterations") +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/loss_barplot.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm",)

# plot gains 
ggplot(plot.data[which(plot.data$Cat=="Gain"),], aes(x=as.factor(pam50),y=N,fill=as.factor(pam50))) +
    geom_bar(stat="identity") +
    scale_x_discrete(name="PAM50 subtype") +
    scale_y_continuous(name="Significant gain alterations") +
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")

ggsave(filename = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/gain_barplot.pdf", sep =""),
       width = 300,
       height = 300,
       units = "mm",)


# boxplots to plot the # of significant genes per sample
# count number of significant genes 
length(Gain_LUMA$Gene) #1385
length(Gain_LUMB$Gene) #0
length(Loss_LUMA$Gene) #15757
length(Loss_LUMB$Gene) #133


bp.cn.data.meta <- cn.data.meta %>% dplyr::select(-c("chr","band","pos","genome"))

# objects to store results
bp.res.matrix <- data.frame(matrix(ncol = 4, nrow = length(cn.data.meta$hgnc_symbol))) %>% dplyr::rename(sampleID=1,PAM50=2,Nloss=3,Ngain=4)

# main loop testing all genes for gain/loss normal and driver
pb = txtProgressBar(min = 0, max = ncol(bp.cn.data.meta)-1, initial = 0, style = 3) 
for (i in 1:(ncol(bp.cn.data.meta)-1)) {
    setTxtProgressBar(pb,i)
    sample <- colnames(bp.cn.data.meta)[i+1]
    bp.res.matrix$sampleID[i] <- sample
    if (sample %in% LUMA.samples) {
        bp.res.matrix$PAM50[i] <- "LumA"
        # get sample data
        #gain
        s.data.gain <- bp.cn.data.meta %>% dplyr::select(c("hgnc_symbol",sample)) %>% 
            filter(hgnc_symbol %in% Gain_LUMA$Gene)
        Freq <- as.data.frame(table(s.data.gain[,2])) %>% filter(Var1 == 1) %>% pull(Freq)
        if (identical(Freq, integer(0))) { 
            bp.res.matrix$Ngain[i] <- 0 } else {(bp.res.matrix$Ngain[i] <- Freq) 
        }
        #loss
        s.data.loss <- bp.cn.data.meta %>% dplyr::select(c("hgnc_symbol",sample)) %>% 
            filter(hgnc_symbol %in% Loss_LUMA$Gene)
        Freq <- as.data.frame(table(s.data.loss[,2])) %>% filter(Var1 == -1) %>% pull(Freq)
        if (identical(Freq, integer(0))) { 
            bp.res.matrix$Nloss[i] <- 0 } else {(bp.res.matrix$Nloss[i] <- Freq) 
        }
    
    } else if (sample %in% LUMB.samples) {
        bp.res.matrix$PAM50[i] <- "LumB"
        # get sample data
        #gain: there are no significant genes int his group
        # s.data.gain <- bp.cn.data.meta %>% dplyr::select(c("hgnc_symbol",sample)) %>% 
        #     filter(hgnc_symbol %in% Gain_LUMB$Gene)
        # bp.res.matrix$Ngain[i] <- as.data.frame(table(s.data.gain[,2])) %>% filter(Var1 == 1) %>% pull(Freq)
        bp.res.matrix$Ngain[i] <- 0
        #loss
        s.data.loss <- bp.cn.data.meta %>% dplyr::select(c("hgnc_symbol",sample)) %>% 
            filter(hgnc_symbol %in% Loss_LUMB$Gene)
        Freq <- as.data.frame(table(s.data.loss[,2])) %>% filter(Var1 == -1) %>% pull(Freq)
        if (identical(Freq, integer(0))) { 
            bp.res.matrix$Nloss[i] <- 0 } else {(bp.res.matrix$Nloss[i] <- Freq)
        }
    }
}

save(bp.res.matrix, file = paste("~/Desktop/MTP_project/Output/CopyNumber/",cohort,"/cn_gainloss_counts.RData",sep = ""))

# plot

#pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/perc_overall_CN_alterations.pdf", sep =""))
bp.res.matrix$PAM50 <- droplevels(as.factor(bp.res.matrix$PAM50))

plot <- ggplot(bp.res.matrix, aes(x=as.factor(PAM50),y=Ngain,fill=PAM50)) +
    geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
    xlab("PAM50 subtype") +
    ylab("Number of singificant genes") + #scale_y_continuous(limits=c(0,0.0055),breaks = seq(0,0.005,0.001)) + 
    theme(axis.text.x = element_text(size = 30),
          axis.title.x = element_text(size = 35),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 35),
          legend.position = "none")
print(plot)
#dev.off()

##################################################################
# Part 3: spatial driver gain/loss - plot CN driver variation of each group
##################################################################

# calc loss/gain freqs per group
cn.data.meta$freqloss.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x<(-1)))/length(HER2.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x<(-1)))/length(LUMA.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqloss.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x<(-1)))/length(LUMB.samples))*-100) # i add a minus to make it easer for plotting
cn.data.meta$freqgain.H2 <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% HER2.samples)], 1, function(x) (length(which(x>1))/length(HER2.samples))*100)
cn.data.meta$freqgain.LUMA <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMA.samples)], 1, function(x) (length(which(x>1))/length(LUMA.samples))*100)
cn.data.meta$freqgain.LUMB <- apply(cn.data.meta[,(colnames(cn.data.meta) %in% LUMB.samples)], 1, function(x) (length(which(x>1))/length(LUMB.samples))*100)

# plot

pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/PAM50_driver_CN_profiles.pdf", sep =""),onefile = TRUE, height = 21.0, width = 29.7)

#her2
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of driver gain/loss CN alterations in luminal HER2-enriched BC") +
    geom_line(aes(y = freqgain.H2), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.H2), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,5)),
                       labels=c(seq(-100,100,5))) +
    theme(text=element_text(size=27)) #+
# annotate(x=80000, y=c(-25,25), label=c("Loss","Gain"), 
#          geom="text", angle=90, hjust=0.5, size=10, colour=c("darkred","darkgreen")) 

print(plot)


#luma
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of driver gain/loss CN alterations in luminal A BC") +
    geom_line(aes(y = freqgain.LUMA), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.LUMA), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,5)),
                       labels=c(seq(-100,100,5))) +
    theme(text=element_text(size=27)) #+
# annotate(x=95000, y=c(-25,25), label=c("Loss","Gain"), 
#          geom="text", angle=90, hjust=0.5, size=10, colour=c("darkred","darkgreen")) 

print(plot)
#lumb
plot <- ggplot(cn.data.meta, aes(x=genome)) + 
    ggtitle("Frequency of driver gain/loss CN alterations in luminal B BC") +
    geom_line(aes(y = freqgain.LUMB), color = "darkgreen",size=3) + 
    geom_line(aes(y = freqloss.LUMB), color = "darkred",size=3) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=chr.lengths$chr,
                       limits = c(min(cn.data.meta$pos),NA)) +
    scale_y_continuous(name="Frequency (%) \n Loss (red) & Gain (green)",
                       breaks=c(seq(-100,100,5)),
                       labels=c(seq(-100,100,5))) +
    theme(text=element_text(size=27)) #+
# annotate(x=-10000, y=c(-25,25), label=c("Loss","Gain"), 
#          geom="text", angle=90, hjust=0.5, size=10, colour=c("darkred","darkgreen"))

print(plot)

# combined
plot <- ggplot(cn.data.meta, aes(x=genome)) + #ggtitle("Frequency of driver gain/loss CN alterations in luminal BC") +
    geom_line(aes(y = freqgain.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqloss.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqgain.LUMA, color = "green"),size=4) + 
    geom_line(aes(y = freqloss.LUMA, color = "green"),size=4) + 
    geom_line(aes(y = freqgain.LUMB, color = "blue"),size=4) + 
    geom_line(aes(y = freqloss.LUMB, color = "blue"),size=4) + 
    scale_colour_manual(name="PAM50", values = c("#d72e2b", "#07a109", "#2176d5"),
                        labels = c("HER2", "LUMA", "LUMB")) + 
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=1) +
    scale_x_continuous(name="Genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=c(as.character(1:22),"X"),
                       limits = c(0,max(cn.data.meta$genome)+50000000),
                       expand = c(0, 0)) +
    scale_y_continuous(name="Alteration frequency (%)", 
                       breaks=c(seq(-20,50,10)),
                       labels=c(20,10,0,10,20,30,40,50),
                       expand = c(0, 0),
                       limits = c(-15,45)) +
    theme(text=element_text(size=30),
          legend.title = element_blank(),
          axis.title.y = element_text(vjust = 0.5),
          legend.position = c(0.97, 0.90)) +
    annotate(x=45000000,y=c(-10,10), label=c("Loss","Gain"), 
             geom="text", angle=90, hjust=0.5, size=9, colour=c("black","black")) 

print(plot)

ggsave(filename="PAM50_driver_CN_profiles.pdf", #_basal
       width = 720,
       height = 210,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/", sep =""))





##################################################################
# Part 4: check selected gene list
##################################################################
library("readxl")

data <- read_excel("Data/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx", sheet = 4)
length(unique(data$Gene)) #93

important.genes <- unique(data$Gene)

important.genes.data <- res.matrix %>% filter(Gene %in% important.genes)

View(important.genes.data)


# how to indicate in the plot which genes are significant


#################################################
# # Karyoploter approach
# # get data into the right format
# cn.data.meta <- as.data.frame(merge(my.regions,rownames_to_column(cn.data,var = "hgnc_symbol"), by="hgnc_symbol")) %>% dplyr::rename(chr=chromosome_name,start=start_position,end=end_position) %>% relocate(hgnc_symbol, .after=end) 
# # add "chr" to the chromosome positions to solve an error later on when matching with a genome
# cn.data.meta$chr <- sub("^","chr",cn.data.meta$chr)
# 
# # some genes are lost due to no available anno
# nrow(cn.data.meta)
# nrow(cn.data)
# 
# #convert to GRanges obj
# cn.data.meta.gr <- toGRanges(cn.data.meta)
# 
# # get genome
# hg19 <- filterChromosomes(getGenome("hg19"), organism = "hg", chr.type = "canonical")
# 
# # this checks the overlap between genome and cn data regions and sorts them in
# cn.data.meta.gr <- sort(c(cn.data.meta.gr, subtractRegions(hg19, cn.data.meta.gr))) 
# 
# # plot
# pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/sample_CN_profiles.pdf", sep =""),onefile = TRUE, height = 21.0, width = 29.7)
# 
# pb = txtProgressBar(min = 0, max = length(names(mcols(cn.data.meta.gr)))-2, initial = 0, style = 3) 
# # for each sample
# for(i in c(1:(length(names(mcols(cn.data.meta.gr)))-2))) { 
#     setTxtProgressBar(pb,i)
#     if(i == 2) { # test
#         break
#     } 
#     # loop iteration example
#     sampleID <- names(mcols(cn.data.meta.gr))[i+2]
#     samplePAM50 <- anno[which(anno$sampleID==sampleID),]$PAM50
#     sample.data <- cn.data.meta.gr[,c(1,2,i+2)]
#     # rename the metadata column cn for the plotting function
#     names(mcols(sample.data))[3] <- "cn"
#     sample.data <- sample.data[,-c(1,2)] # exclude the other metadata
#     
#     # plot
#     kp <- plotKaryotype("hg19", plot.type = 4, labels.plotter = NULL, main=paste("CN data for sample: ",sampleID,"(",samplePAM50,")",sep=""), cex=3.2)
#     kpAddChromosomeNames(kp, srt=45, cex=2.2)
#     plotCopyNumberCallsAsLines(kp, cn.calls = sample.data, ymin= min(cn.data,na.rm = TRUE)-1, ymax= max(cn.data,na.rm = TRUE)+1, lwd=3, add.axis=FALSE, labels = NA, style = "segments")
#     kpAddChromosomeSeparators(kp, lwd=2, col = "#666666")
#     kpAxis(kp, ymin = min(cn.data,na.rm = TRUE)-1, ymax= max(cn.data,na.rm = TRUE)+1, tick.pos = (min(cn.data,na.rm = TRUE)-1):(max(cn.data,na.rm = TRUE)+1), cex=2.2) # fix axis position
#     kpAddLabels(kp, labels = "log2", cex=3, srt=90, pos=3, label.margin = 0.025)
# }
# 
# dev.off()
# 
# ##################################################################
# # Part 2: spatial gain/loss - plot CNV of each group
# ##################################################################
# 
# # use # plotCopyNumberSummary
# 
# # make a list with all samples for each subgroup
# H2.list <- list()
# LA.list <- list()
# LB.list <- list()
# 
# pb = txtProgressBar(min = 0, max = length(names(mcols(cn.data.meta.gr)))-2, initial = 0, style = 3) 
# # loop to add all relevant data to the lists
# for(i in c(1:(length(names(mcols(cn.data.meta.gr)))-2))) { 
#     setTxtProgressBar(pb,i)
#     # loop iteration example
#     sampleID <- names(mcols(cn.data.meta.gr))[i+2]
#     samplePAM50 <- anno[which(anno$sampleID==sampleID),]$PAM50
#     sample.data <- cn.data.meta.gr[,c(1,2,i+2)]
#     # rename the metadata column cn for the plotting function
#     names(mcols(sample.data))[3] <- "cn"
#     sample.data <- sample.data[,-c(1,2)] # exclude the other metadata
#     # filter out rows with NA due to plotCopyNumberSummary() not being able to handle NAs
#     sample.data <- sample.data[(!is.na(elementMetadata(sample.data)[,1]))]
#     # append to the correct list
#     if(samplePAM50 == "Her2") {
#         H2.list[length(H2.list)+1] <- sample.data
#     } else if(samplePAM50 == "LumA") {
#         LA.list[length(LA.list)+1] <- sample.data
#     } else if(samplePAM50 == "LumB") {
#         LB.list[length(LB.list)+1] <- sample.data
#     }
# }
# 
# # plot
# pdf(file = paste("~/Desktop/MTP_project/Output/Plots/CN/",cohort,"/PAM50_CN_profiles.pdf", sep =""),onefile = TRUE, height = 21.0, width = 29.7)
# 
# # Her2
# kp <- plotKaryotype("hg19", plot.type = 4, labels.plotter = NULL, main="CN data for luminal HER2-enriched BC", cex=3.2)
# kpAddChromosomeNames(kp, srt=45, cex=2.2)
# plotCopyNumberSummary(kp, H2.list, r1=0.25, direction = "out") # cant handle NAs
# # LumA
# kp <- plotKaryotype("hg19", plot.type = 4, labels.plotter = NULL, main="CN data for luminal A BC", cex=3.2)
# kpAddChromosomeNames(kp, srt=45, cex=2.2)
# plotCopyNumberSummary(kp, LA.list, r1=0.25, direction = "out") # cant handle NAs
# # LumB
# kp <- plotKaryotype("hg19", plot.type = 4, labels.plotter = NULL, main="CN data for luminal B BC", cex=3.2)
# kpAddChromosomeNames(kp, srt=45, cex=2.2)
# plotCopyNumberSummary(kp, LB.list, r1=0.25, direction = "out") # cant handle NAs
# 
# dev.off()
