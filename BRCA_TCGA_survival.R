#****************
#*Script to analyze TCGA data patient survivl correlations to gene expression
#*Fibroblast drug resistance paper
#*Rosy Li
#*
#*Adapted from CD8TrastuzumabHer2BreastCancer/BRCA_RPPA_macrophge_IGFBP2
#*
#*[updated 02/27/24]
#****************

# install.packages("remotes")
# remotes::install_github("grst/immunedeconv")

rm(list=ls())
options(stringsAsFactors = FALSE)

# library(matrixStats) #rowMaxs function

library(immunedeconv)
# library(cgdsr) #load only when requesting cbioportal data
library(ggplot2) #boxplot
library(ggpubr) #boxplot
library(gplots) #heatmap.2
library(dplyr) #%>%
library(corrplot) 
# library(genefu) #load only when rename.duplicate function is required
# library(MCPcounter) #load only when MCPcounter is used
library(clusterProfiler) #enrichGO function
library(survival)

#************
#*Get HER2+ patient ID
#************
data_clinical <- read.table("E:/OneDrive - University of Pittsburgh/Lab_OneDrive/R/TMElab_projects_MacBackup/R_Rosy/ECM_compression_drug_resistance/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt",header = TRUE,fill=TRUE,sep="\t")
rownames(data_clinical) <- make.names(data_clinical$PATIENT_ID)
#filter out all HER2+ patients and match with expression data
#get the HER2+ patient ID
table(data_clinical$SUBTYPE) #78 HER2+
idx_HER2 <- which(data_clinical$SUBTYPE == "BRCA_Her2")
patientID_HER2 <- data_clinical$PATIENT_ID[idx_HER2]
#change patient ID to rnaseq data patient ID format
patientID_HER2 <- gsub("-", ".", patientID_HER2)
patientID_HER2 <- paste0(patientID_HER2,".01")

### data: breast cancer tcga
#txt file from cBioPortal original RNAseq data
#need to select the correct file path
data1<-read.table("E:/OneDrive - University of Pittsburgh/Lab_OneDrive/R/TMElab_projects_MacBackup/R_Rosy/ECM_compression_drug_resistance/brca_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt",
                  header = TRUE,fill=TRUE)%>%as.data.frame #20531*1084

# write.csv(data1,"TCGA_Data_Fibroblast_Paper.csv")


#remove rows with duplicated Hugo Symbols
data1<-data1[!duplicated(data1$Hugo_Symbol), ]%>%as.data.frame #20500*1084
# # rename rows so that all information are kept
# #problem: new names are not recognized when doing enrichment analysis and error reports
# data1$Hugo_Symbol=rename.duplicate(data1$Hugo_Symbol)

#remove blank spaces and NAs
data1[data1==""] <- NA 
data1<-na.omit(data1,cols=seq_along(data1)) #20501*1084
hugoSymbols <- data1$Hugo_Symbol

#prepare data for quanTIseq
data1_quantiseq<-data1[,-(1:2)]%>%as.data.frame
# data1_quantiseq<-data1_quantiseq[na.omit(data1_quantiseq)]

#check row names unique
length(unique(hugoSymbols))==length(hugoSymbols) 

# data2_quantiseq <- data1_quantiseq
# rownames(data2_quantiseq) <- make.names(data1$Hugo_Symbol) 
#use make.names to convert "_" into "." for names.
#colnames cannot include special characters

#convert to numerival matrix and name columns with patient ID and rows with Hugo symbol
data2_quantiseq<-matrix(as.numeric(unlist(data1_quantiseq)),nrow=nrow(data1_quantiseq))
dim(data2_quantiseq)  #20501*1084
colnames(data2_quantiseq)=colnames(data1_quantiseq)
rownames(data2_quantiseq)=make.names(hugoSymbols)


#*********
#*Immune Deconvolution
#*********

### run quanTIseq
{print("Run the following 2 lines for a new immunedeconvolution. Otherwise, open saved csv.")
res2_estimate<-immunedeconv::deconvolute(data2_quantiseq, "estimate")
write.csv(res2_estimate,"ESTIMATE_BRCA_TCGA_1084patients_Immune_Stromal_Signature.csv",row.names = TRUE,col.names = TRUE)
}
### read ESTIMATE results
res2_estimate<-read.csv("ESTIMATE_BRCA_TCGA_1084patients_Immune_Stromal_Signature.csv")

corr_data<-as.matrix(res2_estimate[,-c(1:2)])
rownames(corr_data)<-res2_estimate$cell_type

corr_data_t <- t(corr_data)%>%as.data.frame
hist(corr_data_t$`stroma score`,breaks=100)


#***
#*run this for HER2 patients only
#*
which(c(rownames(corr_data_t)==colnames(data2_quantiseq))==TRUE)
idx_HER2 <- which(colnames(data2_quantiseq)%in%patientID_HER2)
data2_quantiseq <- data2_quantiseq[,idx_HER2]
corr_data_t <- corr_data_t[idx_HER2,]

# write.csv(data2_quantiseq,"TCGA_HER2_patients_only_fibroblast_paper.csv")

#check correlation between each pair of immune cells
# mypar(1,1)
corrplot(cor(corr_data_t),
         method="color",type="upper",
         addCoef.col = "black",
         tl.col="black", tl.srt=45,#text rotation
         diag=FALSE,
         tl.cex=1,#change text font size
         number.cex=1)#change corr coef font size

#***
#*Correlations between gene expression and stroma scores
#*
which(c(rownames(corr_data_t)==colnames(data2_quantiseq))==TRUE)#everyone needs to be true
plot(corr_data_t$`stroma score`,log2(data2_quantiseq["PLK1",]))
library(ggplot2)
library(hrbrthemes)
library(ggpubr)

# Create dummy data
gene <- "PLK1"
data_scatter <- data.frame(
  my_x = corr_data_t$`stroma score`, 
  my_y = log2(data2_quantiseq[gene,])
)

# write.csv(corr_data_t,"ESTIMATE_stroma_score_HER2_patients_only_TCGA_fibroblast_paper.csv")


ggplot(data_scatter, aes(x=my_x, y=my_y)) +
  geom_point() +
  geom_smooth(method=lm , color="red", se=FALSE)+
  xlab("ESTIMATE stroma score")+
  ylab(paste("log2",gene))+
  stat_cor(method="pearson")
  



#******************
#*Survival analysis for all patients or HER2+ patients
#******************
#*data_clinical_patient.txt
#ALL PATIENTS

#Separate by fibroblast score, median
#get the median of M2/M1 ratio
median_data = median(corr_data_t$`stroma score`)
idx_high <- which(corr_data_t$`stroma score`>median_data)
idx_low <- which(corr_data_t$`stroma score`<=median_data)

names_high <- colnames(data2_quantiseq)[idx_high]
names_low <- colnames(data2_quantiseq)[idx_low]


#Match survival data to RNAseq data/grouping info
names_data_clinical <- make.names(data_clinical$PATIENT_ID)
names_high_match <- gsub(".01","",x=names_high)
names_low_match <- gsub(".01","",x=names_low)

#return patient id (not index)
id_high <- intersect(names_high_match,names_data_clinical)
id_low <- intersect(names_low_match,names_data_clinical)

#generate a data frame
#data.frame() including: patient ID, Survival status, Survival months, grouping (high/low)
id <- c(id_high,id_low)
survstatus <- data_clinical[id,"OS_STATUS"]
survstatus %<>%gsub(":LIVING","",.)%>%gsub(":DECEASED","",.)%>%as.numeric
length(which(survstatus%in%c(0,1)))==length(id)#check if all samples have status annotation

survmonth <- data_clinical[id,"OS_MONTHS"]%>%as.numeric()
hist(survmonth,breaks = 100)

stromalgroup <- c(replicate(length(id_high),"high"),replicate(length(id_low),"low"))

survdata <- data.frame(id,survstatus,survmonth,stromalgroup)
survdata <- na.omit(survdata)

#example of survival package: 
{ #separate by M2/M1, mean
  s <- Surv(survdata$survmonth, survdata$survstatus)
  class(s)
  s
  
  #check models
  survfit(s~1)
  survfit(Surv(survmonth, survstatus)~1, data=survdata)
  
  #generate plots
  sfit <- survfit(Surv(survmonth, survstatus)~1, data=survdata)
  summary(sfit)
  sfit <- survfit(Surv(survmonth, survstatus)~stromalgroup, data=survdata)
  sfit
  summary(sfit)
  
  range(survdata$survmonth)
  seqtimes <- seq(0, 300, 20) #seq(from, to, by)
  seqtimes
  
  summary(sfit, times=seqtimes)
  plot(sfit)
  library(survminer)
  ggsurvplot(sfit,data=survdata)
  
  ggsurvplot(sfit,data=survdata, conf.int=TRUE, pval=TRUE, #risk.table=TRUE, 
             legend.labs=c("high ESTIMATE stroma score", "low ESTIMATE stroma score"), legend.title="Sample",  
             palette=c("dodgerblue2", "orchid2"), 
             title="Overall Survival", 
             risk.table.height=.15)
}


###

#Separate by single gene, median
#get the median of M2/M1 ratio
"PLK1" %in% rownames(data2_quantiseq)
"SERPINE1" %in% rownames(data2_quantiseq) #PAI1
gene <- "SERPINE1"
median_data <- median(data2_quantiseq[gene,])
idx_high <- which(data2_quantiseq[gene,]>median_data)
idx_low <- which(data2_quantiseq[gene,]<=median_data)

names_high <- colnames(data2_quantiseq)[idx_high]
names_low <- colnames(data2_quantiseq)[idx_low]

#Match survival data to RNAseq data/grouping info
names_data_clinical <- make.names(data_clinical$PATIENT_ID)
names_high_match <- gsub(".01","",x=names_high)
names_low_match <- gsub(".01","",x=names_low)

#return patient id (not index)
id_high <- intersect(names_high_match,names_data_clinical)#144
id_low <- intersect(names_low_match,names_data_clinical)#148

#generate a data frame
#data.frame() including: patient ID, Survival status, Survival months, grouping (high/low)
id <- c(id_high,id_low)
survstatus <- data_clinical[id,"OS_STATUS"]
survstatus %<>%gsub(":LIVING","",.)%>%gsub(":DECEASED","",.)%>%as.numeric
length(which(survstatus%in%c(0,1)))==length(id)#check if all samples have status annotation

survmonth <- data_clinical[id,"OS_MONTHS"]%>%as.numeric
hist(survmonth,breaks = 100)

survivalgroup <- c(replicate(length(id_high),"high"),replicate(length(id_low),"low"))

survdata <- data.frame(id,survstatus,survmonth,survivalgroup)
survdata <- na.omit(survdata)

{ 
  #generate plots
  sfit <- survfit(Surv(survmonth, survstatus)~survivalgroup, data=survdata)
  sfit
  summary(sfit)
  
  range(survdata$survmonth)
  seqtimes <- seq(0, 300, 10) #seq(from, to, by)
  seqtimes
  
  summary(sfit, times=seqtimes)
  plot(sfit)
  library(survminer)
  ggsurvplot(sfit,data=survdata)
  
  ggsurvplot(sfit,data=survdata, conf.int=TRUE, pval=TRUE, #risk.table=TRUE, 
             legend.labs=c(paste("high",gene), paste("low",gene)), legend.title="Sample",  
             palette=c("dodgerblue2", "orchid2"), 
             title="Overall Survival", 
             risk.table.height=.15)
}























###*******
#*Survival plots for multiple groups
#*key: change the "grouping" 

###

#Separate by IGFBP2 and M2/M1 ratio
#grouping: samples that have high or low IGFBP2 and high or low M2/M1
idx_highhigh <- intersect(idx_highIGFBP2,idx_highm2m1)
idx_highlow <- intersect(idx_highIGFBP2,idx_lowm2m1)
idx_lowhigh <- intersect(idx_lowIGFBP2,idx_highm2m1)
idx_lowlow <- intersect(idx_lowIGFBP2,idx_lowm2m1)

names_highhigh <- colnames(data2_quantiseq)[idx_highhigh]
names_highlow <- colnames(data2_quantiseq)[idx_highlow]
names_lowhigh <- colnames(data2_quantiseq)[idx_lowhigh]
names_lowlow <- colnames(data2_quantiseq)[idx_lowlow]

#Match survival data to RNAseq data/grouping info
names_data_clinical <- make.names(data_clinical$PATIENT_ID)

names_highhigh_match <- gsub(".01","",x=names_highhigh)
names_highlow_match <- gsub(".01","",x=names_highlow)
names_lowhigh_match <- gsub(".01","",x=names_lowhigh)
names_lowlow_match <- gsub(".01","",x=names_lowlow)


#return patient id (not index)
id_highhigh <- intersect(names_highhigh_match,names_data_clinical)
id_highlow <- intersect(names_highlow_match,names_data_clinical)
id_lowhigh <- intersect(names_lowhigh_match,names_data_clinical)
id_lowlow <- intersect(names_lowlow_match,names_data_clinical)


#generate a data frame
#data.frame() including: patient ID, Survival status, Survival months, grouping (high/low)
id <- c(id_highhigh,id_highlow,id_lowhigh,id_lowlow)
survstatus <- data_clinical[id,"OS_STATUS"]
survstatus %<>%gsub(":LIVING","",.)%>%gsub(":DECEASED","",.)%>%as.numeric
length(which(survstatus%in%c(0,1)))==length(id)#check if all samples have status annotation

survmonth <- data_clinical[id,"OS_MONTHS"]%>%as.numeric
hist(survmonth,breaks = 100)

m2m1group <- c(replicate(length(id_highhigh),"highhigh"),replicate(length(id_highlow),"highlow"),
               replicate(length(id_lowhigh),"lowhigh"),replicate(length(id_lowlow),"lowlow"))

survdata <- data.frame(id,survstatus,survmonth,m2m1group)
survdata <- na.omit(survdata)

{ 
  s <- Surv(survdata$survmonth, survdata$survstatus)
  class(s)
  s
  
  #check models
  survfit(s~1)
  survfit(Surv(survmonth, survstatus)~1, data=survdata)
  
  #generate plots
  sfit <- survfit(Surv(survmonth, survstatus)~1, data=survdata)
  summary(sfit)
  sfit <- survfit(Surv(survmonth, survstatus)~m2m1group, data=survdata)
  sfit
  summary(sfit)
  
  range(survdata$survmonth)
  seqtimes <- seq(0, 185, 10) #seq(from, to, by)
  seqtimes
  
  summary(sfit, times=seqtimes)
  plot(sfit)
  library(survminer)
  ggsurvplot(sfit,data=survdata)
  
  ggsurvplot(sfit,data=survdata, conf.int=TRUE, pval=TRUE, #risk.table=TRUE, 
             legend="right",
             legend.labs=c("high IGFBP2, high M2/M1","high IGFBP2, low M2/M1", "low IGFBP2, high M2/M1","low IGFBP2, low M2/M1"), legend.title="Sample",  
             palette=c("dodgerblue2", "orchid2","green","purple"), 
             title="Overall Survival", 
             risk.table.height=.15)
}




