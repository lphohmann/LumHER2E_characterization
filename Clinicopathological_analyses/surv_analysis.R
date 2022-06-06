# Script: Survival analyses in the Metabric and clinical SCANB cohort 

# TODO: fix the output directory paste()

# empty environment
rm(list=ls())

# set working directory to where the data is
setwd("~/Desktop/MTP_project")

#packages
library(ggplot2)
library(ggfortify)
library(survival)
library(tidyverse)
library(survminer)
library(grid)

#######################################################################
# 1. set parameters
#######################################################################

# 1.1 input for which cohort the SA should be run
cohort <- "SCANB" # Metabric or SCANB

# 1.2 input desired outcome measure for SA
# available OMs: 
# SCANB -> OS, RFI
# Metabric -> OS, DSS, DRFI, RFI, IDFS

OM <- "RFI" 
OM_event <- paste(OM, "bin", sep = "")

#######################################################################
# 2. Cohort-specific data preprocessing including selection of  
# the clinical ER+Her2- subtyped samples
#######################################################################

# for metabric cohort
if (cohort=="Metabric") {
    # 2.1 load data
    load("Data/Metabric/Annotations/Merged_annotations.RData")
    # 2.2 extract relevant variables
    survdata <- anno %>% select(METABRIC_ID, Age, Grade, TumSize,
                                lymph_nodes_positive, ClinGroup, 
                                PAM50, 
                                OS, OSbin, DSS, DSSbin, 
                                DRFI, DRFIbin, RFI, RFIbin, 
                                IDFS, IDFSbin, 
                                Chemotherapy, Endocrine)
    # 2.3 add treatment column
    survdata <- survdata %>% mutate(Treatment = case_when(Chemotherapy==1 & Endocrine==1 ~ "EC",
                                                  Endocrine==1 ~ "E"))
    # 2.4 change lymph node variable to status
    survdata <- survdata %>% mutate(LN = ifelse(lymph_nodes_positive > 0, "N+", "N0"))
    survdata$lymph_nodes_positive <- NULL
    # 2.5 change outcome data to days instead of years
    survdata <- survdata %>% mutate(across(c(OS,DSS,DRFI,RFI,IDFS), (function(years) return(years*365))))
    # 2.6 filter to only include cERpHER2n subjects
    survdata <- survdata %>% filter(grepl('ERpHER2n', ClinGroup))
    # 2.7 getting correct strucutre
    survdata$DSS <- as.numeric(survdata$DSS)
    survdata$DSSbin <- as.numeric(survdata$DSSbin)
    survdata$DRFI <- as.numeric(survdata$DRFI)
    survdata$DRFIbin <- as.numeric(survdata$DRFIbin)
    survdata$IDFS <- as.numeric(survdata$IDFS)
    survdata$IDFSbin <- as.numeric(survdata$IDFSbin)
    
    # for SCANB cohort
    } else if (cohort=="SCANB") {
    # 2.1 load data
    load("Data/SCAN_B/Summarized_SCAN_B_rel4_with_ExternalReview_Bosch_data.RData")
    # 2.2 extract relevant variables
    survdata <- pam50.frame %>% dplyr::select(rba_rel4, fuV8,
                                       PAM50_NCN_ProSigna_rel4,
                                       treatment_Bosch, ERpHER2n_Bosch,
                                       Chemo_Bosch, OS, OSbin,
                                       TumSize_Bosch, Age, NHG_Bosch,
                                       LNstatus_Bosch, relapse_Bosch, 
                                       relapseTime_Bosch) 
    # 2.3 rename columns to match metabric annotation
    survdata <- survdata %>% dplyr::rename(RFI = relapseTime_Bosch,
                                    RFIbin = relapse_Bosch, 
                                    TumSize = TumSize_Bosch,
                                    PAM50 = PAM50_NCN_ProSigna_rel4,
                                    Grade = NHG_Bosch,
                                    LN = LNstatus_Bosch)
    # 2.4 filter to only include ERpHER2n subjects suited for outcome analysis
    survdata <- survdata %>% filter(ERpHER2n_Bosch==1 & fuV8==1)
    # 2.5 add treatment column (ASK: IS IT EC OR ONLY C??)
    survdata <- survdata %>% mutate(Treatment = case_when(Chemo_Bosch==1 ~ "EC",
                                                          Chemo_Bosch==0 & treatment_Bosch==1 ~ "E"))
} 

# 2.6 getting correct structure for common variables
survdata$PAM50 <- as.factor(survdata$PAM50)
survdata$Age <- as.numeric(survdata$Age)
survdata$TumSize <- as.numeric(survdata$TumSize)
survdata$Grade <- as.factor(survdata$Grade) 
survdata$LN <- as.factor(survdata$LN) 
survdata$LN <- relevel(survdata$LN, ref = "N0")
# outcome measures
survdata$OS <- as.numeric(survdata$OS)
survdata$OSbin <- as.numeric(survdata$OSbin)
survdata$RFI <- as.numeric(survdata$RFI)
survdata$RFIbin <- as.numeric(survdata$RFIbin)

#######################################################################
# 2. Defining the PAM50 subtypes of interest
#######################################################################

# 2.1 filter to only include subjects that are PAM50 == Her2 | LumA | LumB
survdata <- survdata %>% filter(PAM50 %in% c("LumA", "LumB", "Her2")) # remove basal
survdata$PAM50 <- droplevels(survdata$PAM50) # drop empty levels
barplot(table(survdata$PAM50))

#######################################################################
# 3. general KM 
#######################################################################

# 3.1 survival object and model
survdata.surv <- Surv(survdata[[OM]], survdata[[OM_event]])
survdata.fit <- survfit(survdata.surv~1, data=survdata, conf.type="log-log")

# 3.2 Check assumptions
# Censoring is random -> not fulfilled, problem?
t.test(survdata[[OM]]~survdata$Treatment)
fit1 <- lm(survdata[[OM_event]]~survdata$PAM50) 
anova(fit1)
summary(fit1)

# 3.3 plot 
autoplot(survdata.fit) + 
    xlab(paste(OM,"days",sep=" ")) + # not days but what?
    ylab(OM) +
    theme_bw()

#######################################################################
# 4. Effect of treatment (log-rank model) 
#######################################################################

# 4.1 Construct model
survdata.treatment <- survdiff(survdata.surv~Treatment, data=survdata)

# 4.2 results
print(survdata.treatment) # depends on which OM is used

# 4.3 plot
autoplot(
    survfit(survdata.surv~Treatment, data=survdata, conf.type="log-log")) + 
    xlab(paste(OM,"days",sep=" ")) +
    ylab(OM) +
    scale_colour_manual(values=c("green", "purple")) +
    scale_fill_manual(values=c("green", "purple")) +
    theme_bw()

# 4.4 save plot
ggsave(filename=paste(cohort,OM,"treatment_KM.pdf",sep="_"),
       width = 297,
       height = 210,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/SA/",cohort,"/", sep =""))


# alternative 
fit <- survfit(survdata.surv~Treatment, data=survdata, conf.type="log-log")
ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE,
            risk.table = FALSE, # Add risk table
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            surv.median.line = "hv", # Specify median survival
            ggtheme = theme_bw(), # Change ggplot2 theme
            palette = c("#E7B800", "#2E9FDF"))

#######################################################################
# 4. Relevel to use HER2 as base for comparison
#######################################################################

# relevel and check
levels(survdata$PAM50)
survdata$PAM50 <- relevel(survdata$PAM50, ref = "Her2")
levels(survdata$PAM50)

#######################################################################
# 5. Investigate the EC treatment group
#######################################################################

# define the group
EC_group <- survdata %>% filter(Treatment == "EC")
EC_group.surv <- Surv(EC_group[[OM]], EC_group[[OM_event]])
table(EC_group$PAM50) 

##########################

# 5.1 Univariate Cox proportional hazards model

# 5.1.1 Model construction
EC_main.pam50 <- coxph(EC_group.surv~PAM50, data=EC_group)

# 5.1.2 Checking assumptions
# mention that prop. Hazard ratio is not fulfilled
cox.zph(EC_main.pam50, transform="km", global=TRUE)
plot(cox.zph(EC_main.pam50, transform="km", global=TRUE))

# 5.1.3 Result
cres <- summary(EC_main.pam50)
round(cres$coefficients[9],5) #lumA
round(cres$coefficients[10],5) #lumB
ggforest(EC_main.pam50,fontsize = 1)
# 5.1.4 Plot

# add count
# add column n with counts for each group
EC_group <- EC_group %>% add_count(PAM50) %>% mutate(PAM50_count = paste0(PAM50, ' (', n, ')'))

# KM 
# autoplot(survfit(EC_group.surv~PAM50_count, data=EC_group, conf.type="log-log"),conf.int = FALSE,censor.size = 11,surv.size = 8,main = "Survival analysis in the chemotherapy + endocrine therapy treatment group") + #censor.size = 8,surv.size = 5
#     xlab("Relapse-free interval (days)") +
#     ylab("Relapse-free interval probability") + 
#     labs(color='PAM50 subtype') +
#     theme(plot.title = element_text(size=22),
#           legend.position= "none", #c(0.85,0.90), # include or exlcude legend
#           legend.title = element_text(size=20), #20
#           legend.key.size = unit(0.5,"cm"), 
#           legend.text = element_text(size = 20), #20
#           axis.text.x = element_text(size = 30), #20
#           axis.title.x = element_text(size = 35), #25
#           axis.text.y = element_text(size = 30), #20
#           axis.title.y = element_text(size = 35)) #25

 
fit <- survfit(EC_group.surv~PAM50, data=EC_group, conf.type="log-log")

survdiff(EC_group.surv ~ PAM50, data = EC_group) # p=0.05 scnab, p=0.004 metabric

plot <- ggsurvplot(
    fit,
    censor.size = 0,
    censor.color = "black",
    size = 8,
    risk.table = FALSE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = FALSE,         # show confidence intervals for 
    xlim = c(0,3500),         # present narrower X axis
    xlab = "Relapse-free interval (days)",
    ylab = "Relapse-free interval probability",
    ylim = c(0.4,1),
    legend= "none",
    break.time.by = 500,     # break X axis in time intervals by 500.
    ggtheme = theme(legend.position= "none", #c(0.85,0.90), # include or exlcude legend
                    legend.title = element_text(size=20), #20
                    legend.key.size = unit(0.5,"cm"), 
                    legend.text = element_text(size = 20), #20
                    axis.text.x = element_text(size = 30), #20
                    axis.title.x = element_text(size = 35), #25
                    axis.text.y = element_text(size = 30), #20
                    axis.title.y = element_text(size = 35)),
    title="Survival analysis in the chemotherapy + endocrine therapy treatment group") 

plot %++% scale_y_continuous(name="Relapse-free interval probability",
                breaks=seq(0.4,1,0.1),
                limits = c(0.4,1))
    
# annotation_custom(grobTree(textGrob("Her2 - LumA: p = 0.004 ", x=0.1,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))) +
    # annotation_custom(grobTree(textGrob("Her2 - LumB: p = 0.08 ", x=0.1,  y=0.05, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic"))))

# save
ggsave(filename=paste(cohort,OM,"EC_PAM50_KM.pdf",sep="_"), #_basal
       width = 325,
       height = 210,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/SA/",cohort,"/", sep =""))

# forest 
ggforest(EC_main.pam50,fontsize = 1)

##########################

# 5.2 Multivariate Cox proportional hazards model

# 5.2.1 Model construction 
# parameters to incl: PAM50, Age, Grade, TumSize
EC_main.all <- coxph(EC_group.surv~PAM50+Age+LN+TumSize+Grade, data=EC_group) 

# 5.2.2 Checking assumptions
# Check for violation of proportional hazard (constant HR over time)
cox.zph(EC_main.all, transform="km", global=TRUE) 
plot(cox.zph(EC_main.all, transform="km", global=TRUE))

# 5.2.3 Result
summary(EC_main.all) 

# 5.2.4 Plot forest 
ggforest(EC_main.all,fontsize = 3,cpositions = c(0.01,0.13,0.35))

ggsave(filename=paste(cohort,OM,"EC_forest.pdf",sep="_"),
       width = 560,
       height = 480,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/SA/",cohort,"/", sep =""))

#######################################################################
# 6. Investigate the EndoCyto treatment group using LumA+B as the reference
#######################################################################
# 
# # 6.1 Show no significant difference between LumA and B in Endo group
# summary(EC_main.all) # depends on chosen endpoint
# 
# ##########################
# 
# # 6.2 Defining LumA+B
# EC_group <- EC_group %>% mutate(PAM50_grouped = if_else(PAM50 == "Her2", "Her2", "LumA+B"))
# 
# ##########################
# 
# # 6.3 Setting LumA+B as reference
# EC_group$PAM50_grouped <- as.factor(EC_group$PAM50_grouped)
# EC_group$PAM50_grouped <- relevel(EC_group$PAM50_grouped, ref = "LumA+B")
# 
# ##########################
# 
# # 6.4 Cox analyses
# 
# # 6.4.1 Univariate
# 
# # 6.4.1.1 model and result
# EC_main.AB.pam50 <- coxph(EC_group.surv~PAM50_grouped, data=EC_group)
# summary(EC_main.AB.pam50)
# 
# # 6.4.1.2 plots
# autoplot(survfit(EC_group.surv~PAM50_grouped, data=EC_group, conf.type="log-log")) + # conf.int = FALSE
#     xlab(paste(OM,"days",sep=" ")) +
#     ylab(OM) +
#     theme_bw()
# # save
# ggsave(filename=paste(cohort,OM,"EC_PAM50_LumAB_KM.pdf",sep="_"),
#        width = 297,
#        height = 210,
#        units = "mm",
#        path = "~/Desktop/MTP_project/Output/Plots/SA/SCANB/")
# 
# ggforest(EC_main.AB.pam50)
# 
# ##########################
# 
# # 6.4.2 Multivariate
# # 6.4.2.1 model and result
# EC_main.AB.all <- coxph(EC_group.surv~PAM50_grouped+Age+LN+TumSize+Grade, data=EC_group)
# summary(EC_main.AB.all)
# 
# # 6.4.2.2 plot 
# ggforest(EC_main.AB.all)
# # save
# ggsave(filename=paste(cohort,OM,"EC_LumAB_forest.pdf",sep="_"),
#        width = 297,
#        height = 210,
#        units = "mm",
#        path = paste("./Plots/SA/",cohort,sep=""))

#######################################################################
# 7. Investigate the Endo treatment group
#######################################################################

# define the group
E_group <- survdata %>% filter(Treatment == "E")
E_group.surv <- Surv(E_group[[OM]], E_group[[OM_event]])
table(E_group$PAM50)

##########################

# 7.1 Univariate Cox proportional hazards model

# 7.1.1 Model construction
E_main.pam50 <- coxph(E_group.surv~PAM50, data=E_group)

# 7.1.2 Checking assumptions
cox.zph(E_main.pam50, transform="km", global=TRUE)

# 7.1.3 Result
cres <- summary(E_main.pam50)
cres
round(cres$coefficients[9],5) #luma
round(cres$coefficients[10],5) #lumb

# 7.1.4 Plots

# add column n with counts for each group
E_group <- E_group %>% add_count(PAM50) %>% mutate(PAM50_count = paste0(PAM50, ' (', n, ')'))

# KM
# autoplot(survfit(E_group.surv~PAM50_count, data=E_group, conf.type="log-log"),conf.int = FALSE,censor.size = 11,surv.size = 8,main = "Survival analysis in the endocrine therapy treatment group")  +
#     labs(color='PAM50 subtype') +
#     ylab("Relapse-free interval probability") +
#     xlab("Relapse-free interval (days)") +
#     theme(plot.title = element_text(size=22),
#           legend.position= "none", #c(0.85,0.90), # include or exlcude legend
#           legend.title = element_text(size=20), #20
#           legend.key.size = unit(0.5,"cm"), 
#           legend.text = element_text(size = 20), #20
#           axis.text.x = element_text(size = 30), #20
#           axis.title.x = element_text(size = 35), #25
#           axis.text.y = element_text(size = 30), #20
#           axis.title.y = element_text(size = 35)) +  #25 
#     scale_y_continuous(breaks=c(seq(40,100,10)), labels=c(seq(40,100,10)))




#+
# annotation_custom(grobTree(textGrob("Her2 - LumA: p = 0.006 ", x=0.1,  y=0.1, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic")))) +
#     annotation_custom(grobTree(textGrob("Her2 - LumB: p = 0.3 ", x=0.1,  y=0.05, hjust=0, gp=gpar(col="black", fontsize=20, fontface="italic"))))

fit <- survfit(E_group.surv~PAM50, data=E_group, conf.type="log-log")
survdiff(E_group.surv ~ PAM50, data = E_group) # p<0.00001 scanb, p= 0.0004 metabric

plot <- ggsurvplot(
    fit,
    censor.size = 0,
    size = 8,
    risk.table = FALSE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = FALSE,         # show confidence intervals for 
    xlim = c(0,3500),         # present narrower X axis
    xlab = "Relapse-free interval (days)",
    ylab = "Relapse-free interval probability",
    ylim = c(0.4,1),
    legend= "none",
    break.time.by = 500,     # break X axis in time intervals by 500.
    ggtheme = theme(legend.position= "none", #c(0.85,0.90), # include or exlcude legend
                    legend.title = element_text(size=20), #20
                    legend.key.size = unit(0.5,"cm"), 
                    legend.text = element_text(size = 20), #20
                    axis.text.x = element_text(size = 30), #20
                    axis.title.x = element_text(size = 35), #25
                    axis.text.y = element_text(size = 30), #20
                    axis.title.y = element_text(size = 35)),
    title="Survival analysis in the endocrine therapy treatment group") 

plot %++% scale_y_continuous(name="Relapse-free interval probability",
                             breaks=seq(0.4,1,0.1),
                             limits = c(0.4,1))


# save
ggsave(filename=paste(cohort,OM,"E_PAM50_KM.pdf",sep="_"), #_basal
       width = 325,
       height = 210,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/SA/",cohort,"/", sep =""))

# forest 
ggforest(E_main.pam50,fontsize = 1)

##########################

# 7.2 Multivariate Cox proportional hazards model
# 7.2.1 Model construction
# parameters to incl: PAM50, Age, NHG, LN, Size
E_main.all <- coxph(E_group.surv~PAM50+Age+LN+TumSize+Grade, data=E_group) 

# 7.2.2 Checking assumptions
# Check for violation of proportional hazard (constant HR over time)
cox.zph(E_main.all, transform="km", global=TRUE) 
plot(cox.zph(E_main.all, transform="km", global=TRUE))

# 7.2.3 Result
summary(E_main.all) 

# 7.2.4 Plot forest 
ggforest(E_main.all,fontsize = 3,cpositions = c(0.01,0.13,0.35))
# save
ggsave(filename=paste(cohort,OM,"E_forest.pdf",sep="_"),
       width = 560,
       height = 480,
       units = "mm",
       path = paste("~/Desktop/MTP_project/Output/Plots/SA/",cohort,"/", sep =""))

#######################################################################
#######################################################################





