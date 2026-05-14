#*****************
#*Script for the analysis of TCGA data
#*HER2 specific genes on macrophage polarization/infiltration
#*Rosy Li
#*2023-01-16
#*****************
#*

rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2) #boxplot
library(ggpubr) #boxplot
# library(gplots) #heatmap.2
library(dplyr) #%>%
library(limma)
library(survival)
library(dplyr)
### data: breast cancer tcga
{
  #txt file from cBioPortal original RNA data
  data1<-read.table("/Users/rosyli/Desktop/R_Rosy/ECM_compression_drug_resistance/brca_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt",
                    header = TRUE,fill=TRUE)%>%as.data.frame #20531*1084
  
  ## remove rows with duplicated Hugo Symbols
  data1<-data1[!duplicated(data1$Hugo_Symbol), ]%>%as.data.frame #20530*1084
  # # rename rows so that all information are kept
  # #problem: new names are not recognized when doing enrichment analysis and error reports
  # data1$Hugo_Symbol=rename.duplicate(data1$Hugo_Symbol)
  
  
  #remove blank spaces and NAs
  data1[data1==""] <- NA 
  data1<-na.omit(data1,cols=seq_along(data1)) #20501*1084
  
  #expression normalization
  library(preprocessCore)
  data1_normalized<-normalize.quantiles(data1[,-c(1:2)]%>%as.matrix,copy=TRUE)
}


###Selection of HER2+ patients only
{
  clinicAnnotation <- read.table("/Users/rosyli/Desktop/R_Rosy/ECM_compression_drug_resistance/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
                                 header = TRUE,fill=TRUE,sep="\t")%>%as.data.frame
  
  #check the order of Clinical annotation and RNAdata
  newcolnames <- gsub("-01","",gsub("\\.","-",colnames(data1)[-c(1:2)]))
  which((newcolnames==clinicAnnotation$PATIENT_ID)==FALSE)
  colnames(data1) <- c(colnames(data1)[c(1:2)],newcolnames)
  
  #not in the same order, re-arrange
  rownames(clinicAnnotation) <- clinicAnnotation$PATIENT_ID
  clinicAnnotation <- clinicAnnotation[newcolnames,]
  which((newcolnames==clinicAnnotation$PATIENT_ID)==FALSE)
  #should return integer(0)
  
  #select HER2+ patients
  idx_her2patient <- which(clinicAnnotation$SUBTYPE=="BRCA_Her2")
  
  #78 patients selected (patient RNAseq starts from the 3rd column)
  # data1 <- data1[,c(1,2,(idx_her2patient+2))]
  data_her2 <- data1[,c(1,2,(idx_her2patient+2))]
  
}

{#select patients with ER+ or TNBC
  #select luminal A patients
  idx_LumApatient <- which(clinicAnnotation$SUBTYPE=="BRCA_LumA")
  #499 patients
  data_LumA <- data1[,c(1,2,(idx_LumApatient+2))]
  
  #select luminbal B patients
  idx_LumBpatient <- which(clinicAnnotation$SUBTYPE=="BRCA_LumB")
  #197 patients
  data_LumB <- data1[,c(1,2,(idx_LumBpatient+2))]
  
  #select basal-like patients
  idx_Basalpatient <- which(clinicAnnotation$SUBTYPE=="BRCA_Basal")
  #171 patients
  data_Basal <- data1[,c(1,2,(idx_Basalpatient+2))]
  
}

res2_quantiseq<-read.csv("quantiseq_BRCA_TCGA_1084patients_11immuneCellSignature.csv")

#run this for HER2, lumA, lumB, basal patients only
{
  # res2_quantiseq <- res2_quantiseq[,c(1,2,(idx_her2patient+2))]
  res2_quantiseq_her2 <- res2_quantiseq[,c(1,2,(idx_her2patient+2))]
  res2_quantiseq_lumA <- res2_quantiseq[,c(1,2,(idx_LumApatient+2))]
  res2_quantiseq_lumB <- res2_quantiseq[,c(1,2,(idx_LumBpatient+2))]
  res2_quantiseq_basal <- res2_quantiseq[,c(1,2,(idx_Basalpatient+2))]
}

{
#correlations between cell types
  library(corrplot)
  corrdata <- t(res2_quantiseq_basal[,-c(1:2)])
  colnames(corrdata) <- res2_quantiseq_her2$cell_type
  corr_plot <- cor(corrdata)
  corrplot(corr_plot, type = 'lower')
  }

# #run this for ER+ patients only
# {
#   idx_her2patient <- c(which(clinicAnnotation$SUBTYPE=="BRCA_LumA"),
#                        which(clinicAnnotation$SUBTYPE=="BRCA_LumB"))
#   
#   data_her2 <- data1[,c(1,2,(idx_her2patient+2))]
#   
#   res2_quantiseq_her2 <- res2_quantiseq[,c(1,2,(idx_her2patient+2))]
# }
# 
# #run this for TNBC patients only
# {
#   idx_her2patient <- c(which(clinicAnnotation$SUBTYPE=="BRCA_Basal"))
#   
#   data_her2 <- data1[,c(1,2,(idx_her2patient+2))]
#   
#   res2_quantiseq_her2 <- res2_quantiseq[,c(1,2,(idx_her2patient+2))]
#   
# }

#grouping
subtype_names <- c()
subtype_names[idx_her2patient] <- "her2"
subtype_names[idx_LumApatient] <- "lumA"
subtype_names[idx_LumBpatient] <- "lumB"
subtype_names[idx_Basalpatient] <- "basal"

#identify the genes that are high/low in HER2+ tumors specifically
#identify the genes that showed significance in ANOVA
# df_anova <- c()
# df_anova.p <- c()
# df_anova.mean_her2 <- c()
# df_anova.mean_lumA <- c()
# df_anova.mean_lumB <- c()
# df_anova.mean_basal <- c()

means <- data.frame(her2=c(0),lumA=c(0),lumB=c(0),basal=c(0))
dunnetdf <- data.frame(lumAHer2Diff=c(0),lumBHer2Diff=c(0),basalHer2Diff=c(0),
                       lumAHer2P=c(0),lumBHer2P=c(0),basalHer2P=c(0))
data2 <- data1[,-c(1:2)]
data2 <- data2[which(rowSums(data2)>100),]#non-normalized
rownames(data2) <- data1$Hugo_Symbol[which(rowSums(data2)>100)]
library(DescTools)
for (i in 1:nrow(data2)) {
  data=c()
  data=data.frame(Expr=as.numeric(data2[i,]),group=subtype_names)
  data$group <- factor(data$group,levels=c("her2","lumA","lumB","basal"))
  # fit <- lm(Expr~group,data)
  res <- DunnettTest(data$Expr, data$group,control="her2")$her2
  dunnetdf[i,1:6] <- c(res[,1],res[,4])
  #first column is diff, 4th column is p values
  means[i,1:4] <- unlist(aggregate(data$Expr, list(data$group), FUN=mean)[2])%>%as.numeric
}

her2_anova_df <- cbind(dunnetdf,means)
rownames(her2_anova_df) <- rownames(data2)
write.csv(her2_anova_df,"her2_anova_df_nonnrmalized.csv")

diff_sum <- rowSums(abs(her2_anova_df[,1:3]))
order(diff_sum,decreasing=TRUE)

#what are the genes that are differentially expressed in her2 tumors?
idx_her2diff <- which(her2_anova_df[,4]<0.01 | her2_anova_df[,5]<0.01 | her2_anova_df[,6]<0.01)
names_her2_diff <- rownames(her2_anova_df)[idx_her2diff]
idx_her2_low <- which(her2_anova_df$lumAHer2Diff>0 & her2_anova_df$lumBHer2Diff>0 & her2_anova_df$basalHer2Diff>0)
idx_her2_high <- which(her2_anova_df$lumAHer2Diff<0 & her2_anova_df$lumBHer2Diff<0 & her2_anova_df$basalHer2Diff<0)

intersect(rownames(her2_anova_df)[intersect(idx_her2diff,idx_her2_low)],DDIR_signature)
intersect(idx_her2diff,idx_her2_low)

#intersect between DNA damage signatures and her2 different genes
DDIR <- read.csv("44-gene DNA Damage Immune Response DDIR signature.csv",header=TRUE)
DDIR_signature <- DDIR$Gene.Symbol
intersect(names_her2_diff,DDIR_signature)

M2total_high_cor_HER2 <- read.csv("/Users/rosyli/Desktop/R_Rosy/OVCAmicrophageInfiltration/M2total_high_correlation_genelist_HER2.csv")
intersect(names_her2_diff,M2total_high_cor_HER2$names)

monocyte_high_cor_HER2 <- read.csv("/Users/rosyli/Desktop/R_Rosy/OVCAmicrophageInfiltration/mono_high_correlation_genelist_HER2.csv")
intersect(names_her2_diff,monocyte_high_cor_HER2$names)

M1M2_high_cor_HER2 <- read.csv("/Users/rosyli/Desktop/R_Rosy/OVCAmicrophageInfiltration/M1M2_high_correlation_genelist_HER2.csv")
intersect(names_her2_diff,M1M2_high_cor_HER2$names)


#load genelists and plot in boxplots or histograms
{
  #wnt genelist
  wnt_genelist <- c(data1$Hugo_Symbol[grep("WNT",data_her2$Hugo_Symbol)],
                    "PORCN","TMED2","GPR177","NOTUM","NDP","EVR2","RSPO1","CER1",
                    "SFRP5","FRP1B","SARP3","WIF","SOST","DKK4","IGFBP4"
  )
  # write.csv(wnt_genelist,"wnt_genelist.csv")
  wnt_genelist <- wnt_genelist[wnt_genelist%in%data1$Hugo_Symbol]
  rownames(data_her2) <- data_her2$Hugo_Symbol
  library(tidyverse)
  library(hrbrthemes)
  library(viridis)
  
  p <- list()
  for (i in 1:length(wnt_genelist)){
    df <- c()
    df <- data.frame(expr=data1[data1$Hugo_Symbol==wnt_genelist[i],-c(1:2)]%>%as.numeric,group=subtype_names)
    df <- na.omit(df)
    df$group <- factor(df$group,levels = c("her2","lumA","lumB","basal"))
    p[[i]] <- df %>%
      ggplot( aes(x=group, y=expr, fill=group)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.4) +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme_ipsum() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle(wnt_genelist[i]) +
      xlab("")+
      stat_compare_means(method="t.test",ref.group="her2",hide.ns=TRUE,label="p.signif")
  }
  library(gridExtra)
  do.call(grid.arrange,p)
}

{#plot the genelist levels for 4 sybtypes
  genelist <- c("ITGAL","NFKB1")
  
  subtype_frame <- data.frame(Expr = as.numeric(data1[data1$Hugo_Symbol==genelist[2],-c(1:2)]),subtype=subtyepe_names)
  subtype_frame <- na.omit(subtype_frame)
  subtype_frame$subtype <- factor(subtype_frame$subtype,levels=c("her2","lumA","lumB","basal"))
  
  ggboxplot(data=subtype_frame, x="subtype", y="Expr", ylab=genelist[2]) +
    ylim(c(0,10000))+
    stat_compare_means(method="t.test",ref.group="her2",hide.ns=TRUE,label="p.signif")
  # ggboxplot(data=subtype_frame, x="subtype", y="IGFBP2", select = c("lumA", "her2")) +
  #   stat_compare_means(method = "t.test")+
  #   ylim(c(0,30000))
  
}


#**********
#*Analysis of metabric data
#*

#***
#*Load metabric data
# untar("brca_metabric.tar.gz")
library(readr)
library(dplyr)
data_raw <- read_tsv("brca_metabric/data_mrna_agilent_microarray.txt")

data_raw[1:10,1:10]
hugo_symbol <- data_raw$Hugo_Symbol

entrezid <- data_raw$Entrez_Gene_Id

library(limma)
library(preprocessCore) 
data_metabric<-normalize.quantiles(as.matrix(data_raw[,-c(1:2)]),copy=TRUE)
ncol(data_metabric)#1980 patients
colnames(data_metabric) <- colnames(data_raw)[-c(1:2)]
rownames(data_metabric) <- hugo_symbol
remove(data_raw)

#***
#*read in the clinical data and filter out HER2+ patients

data_clinical <- read_tsv("brca_metabric/data_clinical_patient.txt")
data_clinical <- data_clinical[-c(1:4),]
data_clinical$`Overall Survival Status` <- gsub("\\:.*","",data_clinical$`Overall Survival Status`)%>%as.numeric()

#all RNAseq patients have clinical annotations
which((colnames(data_metabric)%in%data_clinical$`#Patient Identifier`)==FALSE)

#identify her2 patients' ID and select the survival, RNAseq data exclusively
idx_her2 <- data_clinical$`#Patient Identifier`[which(data_clinical$`Pam50 + Claudin-low subtype`=="Her2")]
# MUST check whether the following line work or not.
rownames(data_clinical) <- data_clinical$`#Patient Identifier`

#subset her2 data exclusively
{
  data_her2 <- data_metabric[,idx_her2]
  data_clinical <- data_clinical[idx_her2,]
}

# survival analysis genelist: wnt ligands
wnt_genelist <- c("SERPINE1")
wnt_genelist <- read.csv("wnt_genelist.csv",header = TRUE, row.names = 1)$x
wnt_genelist <- wnt_genelist[wnt_genelist%in%rownames(data_her2)]

#for loop for survival analysis
p_thresh <- 0.05
library(tidyr)
{
  survival_pval_logRank <- c()
  for (i in 1:length(wnt_genelist)){
    curr_gene <- wnt_genelist[i]
    
    if (curr_gene%in%rownames(data_her2)==TRUE){
      curr_gene_expr <- data_her2[curr_gene,]%>%data.matrix #convert to numerical array
      group_expr <- c(replicate(ncol(data_her2),"NA"))
      group_expr[which(curr_gene_expr>mean(curr_gene_expr))] <- "high"
      group_expr[which(curr_gene_expr<=mean(curr_gene_expr))] <- "low"
      
      # #Threashold with quantile instead of mean
      # group_expr[which(curr_gene_expr>=quantile(curr_gene_expr)[4])] <- "high"
      # group_expr[which(curr_gene_expr<=quantile(curr_gene_expr)[2])] <- "low"
      
      
      survdata <- data.frame(os_status = data_clinical$`Overall Survival Status`,
                             os_months = data_clinical$`Overall Survival (Months)`%>%as.numeric,
                             group = group_expr)
      survdata <- survdata[survdata$group!="NA",]
      
      # survdata <- survdata[c(idx_LumApatient,idx_LumBpatient),]#for ER+ only
      # survdata <- survdata[c(idx_Basalpatient),]#for TNBC only
      # survdata <- survdata #for all patients
      if (length(table(survdata$group))>1){
        survdata <- na.omit(survdata)
        survival_pval_logRank[i] <- survdiff(Surv(os_months, os_status) ~ group, data = survdata[,])$p
      }
      else
        survival_pval_logRank[i] <- 1
    }
    else
      survival_pval_logRank[i] <- 1
    
    # plotting block
    if(survival_pval_logRank[i]<0.05){
      sfit <- survfit(Surv(os_months, os_status)~group, data=survdata)
      sfit
      summary(sfit)

      range(survdata$survmonth)
      seqtimes <- seq(0, 185, 10) #seq(from, to, by)
      seqtimes

      summary(sfit, times=seqtimes)
      ggsurvplot(sfit,data=survdata, conf.int=TRUE, pval=TRUE, #risk.table=TRUE,
                 legend.title=curr_gene,
                 palette=c("dodgerblue2", "orchid2"),
                 title="Overall Survival",
                 risk.table.height=.15)%>%print
    }
  }
  os_her2_surv_wnt <- wnt_genelist[which(survival_pval_logRank<p_thresh)] #os for HER2
}

#***
#*Stroma score calculation
#*
library(immunedeconv)
data_immunedeconv <- as.matrix(data_raw[,-c(1:2)])
rownames(data_immunedeconv) <- make.names(data_raw$Hugo_Symbol)
res_estimate <- immunedeconv::deconvolute(data_immunedeconv,"estimate")

write.csv(res_estimate,"ESTIMATE_metabric_all_patients_Immune_Stromal_Signature.csv")

res_estimate<-read.csv("ESTIMATE_metabric_all_patients_Immune_Stromal_Signature.csv")
rownames(res_estimate) <- res_estimate$cell_type
res_estimate <- res_estimate[,-c(1,2)]


#***
#*Stroma score correlations to PLK1 and SERPINE1 expression
#*
gene <- "SERPINE1"
library(ggplot2)
library(ggpubr)
#in all patients that are HER2+:
idx_her2 <- data_clinical$`#Patient Identifier`[which(data_clinical$`Pam50 + Claudin-low subtype`=="Her2")]

res_estimate_her2 <- res_estimate[,gsub("-",".",idx_her2)]
data_her2 <- data_metabric[,idx_her2]

df_gene_stromaclass <- data.frame(
  gene_value = as.numeric(data_her2[gene,]),#pre therapy gene expression levels
  stroma_value = as.numeric(res_estimate_her2["stroma score",])
)

stroma_value_class <- ifelse(df_gene_stromaclass$stroma_value>median(df_gene_stromaclass$stroma_value),"high","low")
df_gene_stromaclass <- cbind(df_gene_stromaclass, stroma_value_class)

df_gene_stromaclass %>%
  ggplot( aes(x=stroma_value_class, y=gene_value, fill=stroma_value_class)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=stroma_value_class))+
  ylab(paste("Expression levels of",gene))+
  stat_compare_means(method="wilcox.test")

ggplot(df_gene_stromaclass, aes(x=stroma_value, y=gene_value)) +
  geom_point()+ 
  geom_smooth(method = "lm", fill = NA)+
  stat_cor()+
  ylab(paste("Expression levels of",gene))

# write.csv(df_gene_stromaclass,"METABRIC_Stroma_score_PLK1_expression_HER2 patients only.csv")

# temp <- read.csv("METABRIC_Stroma_score_SERPINE1_expression_HER2 patients only.csv")





#***
#*Filter the survival-associated genes for HER2+ patients
#* All genes associated with survival

p_thresh <- 0.05
{
  survival_pval_logRank <- c()
  survival_obs_high <- c()
  survival_obs_low <- c()
  for (i in 1:nrow(data_her2)){
    curr_gene <- rownames(data_her2)[i]
    
    if (curr_gene%in%rownames(data_her2)==TRUE){
      curr_gene_expr <- data_her2[curr_gene,]%>%data.matrix #convert to numerical array
      group_expr <- c(replicate(ncol(data_her2),"NA"))
      group_expr[which(curr_gene_expr>mean(curr_gene_expr))] <- "high"
      group_expr[which(curr_gene_expr<=mean(curr_gene_expr))] <- "low"
      
      survdata <- data.frame(os_status = data_clinical$`Overall Survival Status`,
                             os_months = data_clinical$`Overall Survival (Months)`%>%as.numeric,
                             group = group_expr)
      
      # survdata <- survdata[c(idx_LumApatient,idx_LumBpatient),]#for ER+ only
      # survdata <- survdata[c(idx_Basalpatient),]#for TNBC only
      # survdata <- survdata #for all patients
      if (length(table(survdata$group))>1){
        survdata <- na.omit(survdata)
        surv_test <- survdiff(Surv(os_months, os_status) ~ group, data = survdata[,])
        survival_pval_logRank[i] <- surv_test$p
        survival_obs_high[i] <- surv_test$obs[1]
        survival_obs_low[i] <- surv_test$obs[2]
      }
      else
      {
        survival_pval_logRank[i] <- 1
        survival_obs_high[i] <- 0
        survival_obs_low[i] <- 0
      }
        
    }
    else
    {
      survival_pval_logRank[i] <- 1
      survival_obs_high[i] <- 0
      survival_obs_low[i] <- 0
    }
    
    #plotting block
    # if(survival_pval_logRank[i]<0.05){
    #   sfit <- survfit(Surv(os_months, os_status)~group, data=survdata)
    #   sfit
    #   summary(sfit)
    #   
    #   range(survdata$survmonth)
    #   seqtimes <- seq(0, 185, 10) #seq(from, to, by)
    #   seqtimes
    #   
    #   summary(sfit, times=seqtimes)
    #   ggsurvplot(sfit,data=survdata, conf.int=TRUE, pval=TRUE, #risk.table=TRUE, 
    #              legend.title=curr_gene,  
    #              palette=c("dodgerblue2", "orchid2"), 
    #              title="Overall Survival", 
    #              risk.table.height=.15)%>%print
    # }
  }
  os_her2_surv_allgenes <- rownames(data_her2)[which(survival_pval_logRank<p_thresh)] #os for HER2
}

df_os_her2_allgenes <- data.frame(
  geneName = rownames(data_her2),
  survival_pval_logRank,
  survival_highExpr = survival_obs_high, #the length of observed survival for high expression patients
  survival_lowExpr = survival_obs_low #the length of observed survival for low expression patients
)

df_os_hers_allgenes_sig <- df_os_her2_allgenes[which(survival_pval_logRank<p_thresh),]

write.csv(df_os_her2_allgenes,"metabric_df_os_her2_allgenes.csv")
write.csv(df_os_hers_allgenes_sig,"metabric_df_os_hers_allgenes_sig.csv")


#***
#*Immune deconvolution
#*
library(immunedeconv)
library(quantiseqr)

library(EnsDb.Hsapiens.v86)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys= hugo_symbol, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
dupli = duplicated(geneIDs$SYMBOL)
geneIDs = geneIDs[!dupli,]
ind = intersect(hugo_symbol, geneIDs$SYMBOL)
data_quanTIseq = data_metabric[ind,]
rownames(data_quanTIseq) = geneIDs$SYMBOL

data_quanTIseq <- as.matrix(data_quanTIseq)
data_quanTIseq <- na.omit(data_quanTIseq)

{#Run for immune deconvolution
  # res_quanTIseq <- quantiseqr::run_quantiseq(
  #   expression_data = data_quanTIseq,
  #   is_arraydata = TRUE,
  #   is_tumordata = TRUE,
  #   scale_mRNA = TRUE
  # )
  # write.csv(res_quanTIseq,"res_quanTIseq_metabric_brca_1980samples.csv")
  # rownames(res_quanTIseq) <- res_quanTIseq$Sample
  
}
res_quanTIseq <- read.csv("res_quanTIseq_metabric_brca_1980samples.csv",header=TRUE,row.names = 1)
res_quanTIseq <- res_quanTIseq[,-1]



#***
#*Identify the genes correlated to M2 macrophages
quanTIseq_her2 <- res_quanTIseq[idx_her2,]
M2_her2 <- quanTIseq_her2$Macrophages.M2
hist(M2_her2,breaks = 100)
cor_Pearson_M2 <- c()
cor_pval <- c()
for (i in 1:nrow(data_her2)){
  curr_test <- cor.test(M2_her2,data_her2[i,])
  cor_Pearson_M2[i] <- curr_test$estimate[1]
  cor_pval[i] <- curr_test$p.value[1]
}

df_her2_m2 <- data.frame(
  geneName = rownames(data_her2),
  cor_coeff_Pearson = cor_Pearson_M2,
  pval = cor_pval
)
write.csv(df_her2_m2,"metabric_her2_m2_correlations_all_genes.csv")

df_her2_m2_sig <- df_her2_m2[which(df_her2_m2$pval<0.05),]
write.csv(df_her2_m2_sig,"metabric_her2_m2_correlations_sig_genes.csv")


#***
#*Identify the genes correlated to CD8 T cells
quanTIseq_her2 <- res_quanTIseq[idx_her2,]
CD8_her2 <- quanTIseq_her2$T.cells.CD8
hist(CD8_her2,breaks=100)
cor_Pearson <- c()
cor_pval <- c()
for (i in 1:nrow(data_her2)){
  curr_test <- cor.test(CD8_her2,data_her2[i,])
  cor_Pearson[i] <- curr_test$estimate[1]
  cor_pval[i] <- curr_test$p.value[1]
}

df_her2_cd8 <- data.frame(
  geneName = rownames(data_her2),
  cor_coeff_Pearson = cor_Pearson,
  pval = cor_pval
)
write.csv(df_her2_cd8,"metabric_her2_cd8_correlations_all_genes.csv")

df_her2_cd8_sig <- df_her2_cd8[which(df_her2_cd8$pval<0.05),]
write.csv(df_her2_cd8_sig,"metabric_her2_cd8_correlations_sig_genes.csv")


#***
#*Transcription factor analysis
#*The Human Transcription Factor downloaded from:
#*https://www.sciencedirect.com/science/article/pii/S0092867418301065?via%3Dihub
#*Lambert et al. 2018 Cell
#*Supplementary table S1
#*

#Load the transcription factor genelist

TF_csv <- read.csv("human_Transcription_Factor_list_Lambert_2018_Cell.csv",
                        header=TRUE)
TF_genelist <- TF_csv$Name[which(TF_csv$Is.TF.=="Yes")] # 1639 TFs

TcellActivation_TF_genelist <- 

CD8_correlated_TF <- intersect(df_her2_cd8_sig$geneName,TF_genelist)
#the BCL family transcription factors
grep("BCL",CD8_correlated_TF,value=TRUE)
grep("BCL",TF_genelist,value=TRUE)
#Yah!!!! ALL BCL family transcription factors are associated with CD8 infiltration!!
#AND all BCL POSITIVELY correlated with CD8 infiltration!!!

grep("NFKB",CD8_correlated_TF,value=TRUE)
grep("NFKB",TF_genelist,value=TRUE)
#NFKB1 also correlated to CD8 infiltration!!!

CD8_correlated_TF <- df_her2_cd8_sig[df_her2_cd8_sig$geneName%in%TF_genelist,]
write.csv(CD8_correlated_TF,"CD8_correlated_TF_metabric.csv")

#intersect of CD8 correlated TFs and survival correlated genes
df_os_hers_allgenes_sig <- read.csv("metabric_df_os_hers_allgenes_sig.csv",header=TRUE,row.names=1)
intersect(df_os_hers_allgenes_sig$geneName,CD8_correlated_TF$geneName)

#intersect of CD8 correlated TFs and survival correlated genes
chemokines <- read.csv("C:/Users/TME Lab User/Desktop/RL R projects/R_Rosy_macBackUp/OVCAmicrophageInfiltration/geneset_KEGG_cytokine_cytokine_receptor_interaction.txt",
                       header=TRUE)[-1,]
intersect(df_os_hers_allgenes_sig$geneName,chemokines)

#identify genes correlated to M2 macrophages
df_her2_m2_sig <- read.csv("metabric_her2_m2_correlations_sig_genes.csv",header = TRUE,
                           row.names = 1)
m2_correlated_TF <- df_her2_m2_sig[df_her2_m2_sig$geneName%in%TF_genelist,]
grep("BCL",m2_correlated_TF$geneName,value=TRUE)

write.csv(m2_correlated_TF,"m2_correlated_TF_metabric.csv")
















