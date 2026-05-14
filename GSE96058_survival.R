#**************
#*Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
#*Analysis of the GSE96058 survival data
#*
#*Rosy Li 3.3.2024
#**************
#*


# library(GEOquery)
# 
# gse <- getGEO("GSE96058")
# head(Meta(gse))
# gsedata <- gse[[1]]
# data_RNA <- exprs(gsedata)#no RNAseq data in the database. Check the supplementary data later
# data_clinical <- pData(gsedata)
# write.csv(data_clinical,"GSE96058_clinical_annotation.csv")

data_clinical <- read.csv("GSE96058_clinical_annotation.csv",header = TRUE, row.names = 1)

#download the supplementary data
# eList <- getGEOSuppFiles("GSE96058")
# gunzip("GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
library(data.table)
data_RNA <- fread("GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv") #30665 3410
# data_RNA <- fread("GSE96058/GSE96058_transcript_expression_3273_samples_and_136_replicates.csv")

data_RNA[1:10,1:10]

data_RNA$V1
"PLK1"%in%data_RNA$V1
"SERPINE1"%in%data_RNA$V1


row.names(data_RNA) <- make.names(data_RNA$V1, unique=TRUE)
allGenes <- row.names(data_RNA)
data_RNA <- data_RNA[,-1]
{###stratify by patient PAM50 values
  # HER2_patient_ID <- data_clinical$title[which(data_clinical$pam50.subtype.ch1=="Her2")]
# # match the order of RNAseq data to clinical data
# data_clinical_HER2 <- subset(data_clinical,pam50.subtype.ch1=="Her2")
#   library(dplyr)
#   data_RNA_HER2 <- select(data_RNA, HER2_patient_ID)
#   data_clinical_HER2$title==colnames(data_RNA_HER2)#need to be all TRUE
  
  ###stratify by patient HER2 status
  HER2_patient_ID <- data_clinical$title[which(data_clinical$her2.status.ch1==1)]
  # match the order of RNAseq data to clinical data
  data_clinical_HER2 <- subset(data_clinical,her2.status.ch1==1)
  library(dplyr)
  data_RNA_HER2 <- select(data_RNA, HER2_patient_ID)
  data_clinical_HER2$title==colnames(data_RNA_HER2)#need to be all TRUE
}


# #run this for all patients
# data_RNA_HER2 <- data_RNA
# data_clinical_HER2 <- data_clinical

#***
#*Survival analysis
library(survival)
library(survminer)
row.names(data_RNA_HER2) <- allGenes

curr_gene <- "SERPINE1"
curr_gene %in% row.names(data_RNA_HER2)
curr_gene_expr <- data_RNA_HER2[which(row.names(data_RNA_HER2)==curr_gene),]%>%as.numeric
hist(curr_gene_expr, breaks=300)
median_data = median(curr_gene_expr)
# idx_high <- which(curr_gene_expr>quantile(curr_gene_expr)[4])
# idx_low <- which(curr_gene_expr<=quantile(curr_gene_expr)[2])
# idx_high <- which(curr_gene_expr>median_data)
# idx_low <- which(curr_gene_expr<=median_data)
cutoff<-3.94
idx_high <- which(curr_gene_expr>cutoff)
idx_low <- which(curr_gene_expr<=cutoff)

{
  id <- c(idx_high,idx_low)
  survstatus <- data_clinical_HER2$overall.survival.event.ch1[id]
  # survstatus %<>%gsub(":LIVING","",.)%>%gsub(":DECEASED","",.)%>%as.numeric
  length(which(survstatus%in%c(0,1)))==length(id)#check if all samples have status annotation
  survmonth <- data_clinical_HER2$overall.survival.days.ch1[id]
  hist(survmonth,breaks = 100)
  gene_expr_group <- c(replicate(length(idx_high),"high"),replicate(length(idx_low),"low"))
  survdata <- data.frame(id,survstatus,survmonth,gene_expr_group)
}


{ 
  sfit <- survfit(Surv(survmonth, survstatus)~gene_expr_group, data=survdata)
  
  summary(sfit)
  plot(sfit)
  
  ggsurvplot(sfit,data=survdata, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
             legend.labs=c("high", "low"), legend.title="Sample",  
             palette=c("dodgerblue2", "orchid2"), 
             title=paste("Overall Survival",curr_gene), 
             risk.table.height=.25)
}
















