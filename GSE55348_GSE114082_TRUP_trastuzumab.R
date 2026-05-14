#*************
#*Script for CD8 T cell infiltration and trastuzumab responses
#*Deconvolution, Regression, Classification, DEG, Network analysis
#*Rosy Li
#*2021-03-08
#*************

###
rm(list=ls())
options(stringsAsFactors = FALSE) # do not convert strings to factors in dataframes

##T cell GEO dataset
##GSE130788

library(GEOquery)
library(dplyr)
library(glmnet)
library(plotmo)
library(caret)
library(illuminaHumanv4.db)
library(doParallel)
library(e1071)
library(data.table)
library(heatmaply)
library(immunedeconv)
library(ggplot2)
library(reshape2)#melt
library(tidyr)#gather
library(tibble)
library(corrplot)
library(ggpubr)#ggboxplot
library(gprofiler2)#enrichment

#***************************
###test the GSE130788 dataset 
#(phaseII clinical trial for trastuzumab lapatinib combination paper)

# gse <- getGEO("GSE130788", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
# 
# # untar( "GSE130788_RAW.tar")
# test <- read.table("/Users/rosyli/Desktop/R_Rosy/BRCAfibroblastDrugResistance/untarGSE/GSM3753582_US22502571_251485068472_S01_GE2_107_Sep09_2_OFF_1_1.txt.gz",
#                    fill = TRUE , header = FALSE)

#*****************************
###test the GSE123845 dataset
#(neoadjuvant chemotherapy)
# gse1 <- getGEO("GSE123845", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
# 
# GSE123845_exp_tpm_matrix.csv.gz
# readin1 <- fread("GSE123845_exp_tpm_matrix.csv.gz")
# gene_names <- readin1$V1
# readin1 <- readin1[,-1]
# mydata1 <- t(as.matrix(readin1))
# colnames(mydata1) <- gene_names
# log_data1 <- log2(mydata1+1)
# meta_data1 <- pData(gse1)
# 
# sampnames1 <- make.names(meta_data1$title)
# meta_data1 <- meta_data1[match(rownames(mydata1), sampnames1),]# reorder metadata with RNAseq patient order
# 
# ### meta data
# #therapy method: AC(adjuvant chemotherapy), T(Taxol), H(Trastuxumab)
# table(meta_data1$characteristics_ch1.16)
# 
# #response: pCR(pathology complete response) 0:complete response 1:non-complete response
# table(meta_data1$characteristics_ch1.35)

#*****************************
###Validate the chemokine results with the GSE55348 dataset
gse_tras_val <- getGEO("GSE55348",GSEMatrix = TRUE,AnnotGPL = TRUE)[[1]]
meta_data_val <- pData(gse_tras_val)

readin_tras_val <- read.table("E:/OneDrive - University of Pittsburgh/Lab_OneDrive/R/TMElab_projects_MacBackup/R_Rosy/BRCAfibroblastDrugResistance/GSE55348_non-normalized.txt", sep='\t', header=T)[-c(1:2),]
nonnormalizeddata_tras_val <- readin_tras_val[-1,-1]
colnames(nonnormalizeddata_tras_val) <- readin_tras_val[1,-1]
rownames(nonnormalizeddata_tras_val) <- readin_tras_val[-1,1]
nonnormalizeddata_tras_val <- nonnormalizeddata_tras_val[,grep("Signal",colnames(nonnormalizeddata_tras_val))]
colnames(nonnormalizeddata_tras_val)<-gsub(".AVG_Signal","",colnames(nonnormalizeddata_tras_val))

which(colnames(nonnormalizeddata_tras_val)!=meta_data_val$description)

#discard rows with NAs or blanks
nonnormalizeddata_tras_val1<-na.omit(mutate_all(nonnormalizeddata_tras_val, ~ifelse(. %in% c("N/A", "null", ""),  NA, .)))
nonnormalizeddata_tras_val1<-matrix(as.numeric(as.matrix(nonnormalizeddata_tras_val1)),ncol=53)
rownames(nonnormalizeddata_tras_val1) <- rownames(nonnormalizeddata_tras_val)
colnames(nonnormalizeddata_tras_val1) <- colnames(nonnormalizeddata_tras_val)
#nonnormalizeddata_tras_val1 is numeric

#convert ILMN (Illumina probe annotations) to gene sybmols
library("illuminaHumanv4.db")
data.frame(Gene=unlist(mget(x=rownames(nonnormalizeddata_tras_val1),envir=illuminaHumanv4SYMBOL,ifnotfound = NA)))

rownames(nonnormalizeddata_tras_val1) <- unlist(mget(x=rownames(nonnormalizeddata_tras_val1),envir=illuminaHumanv4SYMBOL,ifnotfound = NA))

###data normalization
mydata_nonnormalized_tras_t_val <- matrix(as.numeric(t(nonnormalizeddata_tras_val)),nrow =53)
colnames(mydata_nonnormalized_tras_t_val) <- colnames(mydata_nonnormalized_tras_t)
mydata_log2_tras_val <- log2(mydata_nonnormalized_tras_t_val+1)
colnames(mydata_log2_tras_val) <- rownames(nonnormalizeddata_tras)
rownames(mydata_log2_tras_val) <- colnames(nonnormalizeddata_tras_val)
mydata_normalized_tras_val <- percentize(mydata_log2_tras_val)%>%as.matrix


###Deconvolution
res_estimate_val <- deconvolute(nonnormalizeddata_tras_val1, "estimate")
rownames_to_replace <- res_estimate_val$cell_type
res_estimate_val<- t(res_estimate_val[,-1])%>%as.data.frame
colnames(res_estimate_val) <- rownames_to_replace

hist(res_estimate_val$`stroma score`)

#classify by relapse
#define "relapse" as relapse happends within 600 days
time_to_relspse_val <- gsub("time to relapse: ","",meta_data_val$characteristics_ch1.8)%>%as.numeric
class_relapse_val <- ifelse((time_to_relspse_val>700),"NonRelapse","Relapse")
table(class_relapse_val)

#classify by estimate stroma score: high/low (compared to mean)
medianE_score <- median(res_estimate_val$`stroma score`)
class_E_score <- ifelse(res_estimate_val$`stroma score`>medianE_score,"High","Low")

table(class_relapse_val,class_E_score)

#Boxplots of gene expression in relapsed and non-relapsed patients
genelist <- c("PLK1","SERPINE1")
genelist <- genelist[genelist%in%rownames(nonnormalizeddata_tras_val1)]
library(tidyverse)
library(hrbrthemes)
library(viridis)
# par(mfrow=c(5,6)) 
p <- list()
for (i in 1:length(genelist)){
  df <- c()
  df <- data.frame(expr=nonnormalizeddata_tras_val1[genelist[i],]%>%as.numeric%>%log2,group=class_relapse_val)
  p[[i]] <- df %>%
    ggplot( aes(x=group, y=expr, fill=group)) +
    geom_boxplot() +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(genelist[i]) +
    xlab("")+
    stat_compare_means(method="wilcox.test")
}
library(gridExtra)
do.call(grid.arrange,p)


#boxplot of stroma signature in relapse and non-relapse groups
df <- data.frame(escore=res_estimate_val$`stroma score`,group=class_relapse_val)
ggplot(df, aes(x=group, y=escore, fill=group)) +
  geom_boxplot()+
  stat_compare_means(method="wilcox.test")+
  theme_ipsum()+
  ylab("ESTIMATE stroma score")+
  xlab("")+
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  )




###Correlations
cor2<-cor(res_quantiseq1_val$`T cell CD8+`,res_cibersort1_val$`T cell CD8+`)
plot(res_quantiseq1_val$`T cell CD8+`,res_cibersort1_val$`T cell CD8+`,xlab="quanTIseq CD8 val",ylab="CIBERSORT CD8 val",main=paste("cor=",round(cor2,3)))

#IGFBP2 levels
#Hypoxia updated 1.1.2023
gene = "HIF1A"
# gene="IGFBP2"
repapse = nonnormalizeddata_tras_val1[gene,class_relapse_val=="Relapse"]
nonrelapse = nonnormalizeddata_tras_val1[gene,class_relapse_val=="NonRelapse"]
hist(repapse,main=gene)
hist(nonrelapse, main=gene)
t.test(repapse,nonrelapse)

hypoxia_genelist <- c("HIF1A","EPAS1",
                      "VEGFA","SLC2A1","PGAM1","ENO1","LDHA",
                      "TPI1","P4HA1","MRPS17","CDKN3","ADM",
                      "NDRG1","TUBB6","ALDOA","MIF","ACOT7")
hypoxia_genelist <- hypoxia_genelist[hypoxia_genelist%in%rownames(nonnormalizeddata_tras_val1)]
library(tidyverse)
library(hrbrthemes)
library(viridis)
# par(mfrow=c(5,6)) 
p <- list()
for (i in 1:length(hypoxia_genelist)){
  df <- c()
  df <- data.frame(expr=nonnormalizeddata_tras_val1[hypoxia_genelist[i],]%>%as.numeric,group=class_relapse_val)
  p[[i]] <- df %>%
    ggplot( aes(x=group, y=expr, fill=group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(hypoxia_genelist[i]) +
    xlab("")+
    stat_compare_means(method="wilcox.test")
}
library(gridExtra)
do.call(grid.arrange,p)


#genelist <- c("CCL5")
genelist <- c("CASP3","CCL5")
gene_curr <- c()
genedata_curr <- c()
cor_curr <- c()
p_corr_CD8_cytokine <- c()

#multiple correlations
for(i in 1:length(genelist)){
  gene_curr <- genelist[i]
  genedata_curr <- c(mydata_log2_tras_val[,gene_curr])
  cor_curr <- round(cor(res_quantiseq1_val$`T cell CD8+`, genedata_curr),3)
  p_corr_CD8_cytokine[i] <- plot(res_quantiseq1_val$`T cell CD8+`,genedata_curr, 
                              xlab="T cell fraction quanTIseq",ylab=genelist,
                              main=paste("validation data, cor=",cor_curr))
  
}
#single correlation plot
genedata_val <- c(mydata_normalized_tras_val[,genelist])
cor(res_quantiseq1_val$`T cell CD8+`, genedata_val)
plot(res_quantiseq1_val$`T cell CD8+`,genedata_val, 
     xlab="T cell fraction quanTIseq",ylab=genelist,
     main=paste("validation data, cor=",round(cor(res_quantiseq1_val$`T cell CD8+`,genedata_val),3)))


###Organize data
mydata_T_forLasso_val <- data.frame(condition=as.factor(class_Tcell_val), t(nonnormalizeddata_tras_val1))

#*******************************
#* Plot figures for cytokines
#*******************************
#*
#check cytokine expression
#boxplot_genelist_val <- c("CCR1","CCR2","CCR3","CCR4","CCR5","CXCR3","CCL2","CCL5","CCL7","CCL8","CCL20","CXCL9","CXCL10","CXCL11","CXCL12","XCL1","IL10","EBI3","TGFB1","TGFB2","PDCD1","PTGES")
boxplot_genelist_val <- c("CCR2","CCR4","CCR5","CCL2","CCL5","CCL8","CXCL9","CXCL11")
boxplot_genelist_val%in%colnames(mydata_T_forLasso_val)

###Boxplots for POST-TREATMENT cytokine levels
# ONLY POST-treatment data (row 18:34) are used!!!
data_cytokine_boxplot_val <- c()
for(i in 1:length(boxplot_genelist_val)){
  data_cytokine_boxplot_val <- rbind(data_cytokine_boxplot_val , 
                                 data.frame(
                                   Label=class_Tfraction_val,#Label is the class of the T cell level.
                                   variable=c(replicate(53,boxplot_genelist_val[i])),
                                   value=mydata_T_forLasso_val[,boxplot_genelist_val[i]])
  )
}

#horizontal version
p_validation <- ggplot(data = data_cytokine_boxplot_val, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=Label))+
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  stat_compare_means(aes(group = Label), label = "p.format")#p.format,p.signif
p_validation

#facet version
p_validation + facet_wrap( ~ variable, scales="free") #facet version


#*********
#*DEG for validation dataset
#*********

#Identify DEGs
fc_all_val <- c()
pval_all_val <- c()

for (i in 1:ncol(mydata_nonnormalized_tras_t_val)){
  fc_all_val[i] <- sum(mydata_nonnormalized_tras_t_val[idx_Tpositive_val,i])/sum(mydata_nonnormalized_tras_t_val[idx_Tnegative_val,i]) -1
  pval_all_val[i] <- t.test(mydata_nonnormalized_tras_t_val[idx_Tpositive_val,i],mydata_nonnormalized_tras_t_val[idx_Tnegative_val,i])$p.val
}
length(which(pval_all_val<0.05))#2441

pval_all_adj_val <- p.adjust(pval_all_val,method="BH")

names(pval_all_val) <- colnames(mydata_nonnormalized_tras_t_val)
names(pval_all_adj_val) <- colnames(mydata_nonnormalized_tras_t_val)
names(fc_all_val) <- colnames(mydata_nonnormalized_tras_t_val)

length(which(pval_all_adj_val<0.05))#0
length(which(pval_all_val<0.05 & abs(fc_all_val)>0.2))#1611
length(which(pval_all_val<0.05 & abs(fc_all_val)>0.5))#935
length(which(pval_all_val<0.05 & abs(fc_all_val)>1))#101

#Volcano plot for genes changed 50% and p<0.05

volcanodata_val<-cbind(fc_all_val,pval_all_val)%>%as.data.frame
colnames(volcanodata_val)<-c("fc","pval")

volcanodata_val$delabel <- NA
volcanodata_val$delabel[abs(volcanodata_val$fc)>0.5 & volcanodata_val$pval<0.05] <- rownames(volcanodata_val)[abs(volcanodata_val$fc)>0.5 & volcanodata_val$pval<0.05]

volcanodata_val$trend_val<-replicate(length(fc_all_val),"no")
volcanodata_val$trend_val[volcanodata_val$fc>0.5 & volcanodata_val$pval<0.05]<-"up"
volcanodata_val$trend_val[volcanodata_val$fc<(-0.5) & volcanodata_val$pval<0.05]<-"down"

colnames(volcanodata_val)<-c("fc","pval","lab","trend")

library(ggrepel)

ggplot(volcanodata_val, 
       aes(x=log2(fc+1), y=-log10(pval),label=lab,color=trend))+
  geom_point(cex=1)+
  theme_minimal()+
  #geom_text_repel(cex=0) +
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_vline(xintercept=c(log2(-0.5+1), log2(0.5+1)), col="grey60") +
  geom_hline(yintercept=-log10(0.05), col="grey60")+
  ggtitle(label="Post-trastuzumab therapy, high vs. low")




###******************
#Enrichment analysis
###*****************

###Perform gene ontology enrichment analysis

sig_and_fc_thres <- names(which(pval_all_val<0.05 & fc_all_val>0.2))

# library('org.Hs.eg.db')
# columns(org.Hs.eg.db)
# symbols_all <- rownames(all_gene_symbol)
# mapIds(org.Hs.eg.db, symbols_all, 'ENTREZID', 'UNIGENE')

library("illuminaHumanv4.db")
data.frame(Gene=unlist(mget(x = rownames(all_gene_symbol),envir = illuminaHumanv4SYMBOL)))

ego <- enrichGO(gene = sig_and_fc_thres,                  #differentially expressed genes
                universe = all_gene_symbol,                 #list of all genes
                keyType = "GENENAME",                                #ID used to identify genes
                OrgDb = org.Hs.eg.db, 
                ont = "BP",                                          #select GO term to use
                pAdjustMethod = "BH", 
                qvalueCutoff = 1, 
                pvalueCutoff = 1,
                readable = TRUE)

#need to figure out how to convert symbol names

###Ghose plot for multiple datasets
length(which(abs(fc_all_val)>0.3))
deg_idx_val <- order(abs(fc_all_val))[1:5000]
deg_list_val <- colnames(mydata_normalized_tras_val)[deg_idx_val]

res_val <- gost(deg_list_val, organism="hsapiens")
gostplot(res_val,interactive=TRUE)



mydata_T_forLasso_val <- data.frame(condition=as.factor(class_Tcell_val), mydata_nonnormalized_tras_t_val[c(idx_Tnegative_val,idx_Tpositive_val),])

#*************
# Feature Selection
#*************
#
###Lasso regression for parameter selection
colnames(mydata_normalized_tras_val) <- all_gene_symbol$SYMBOL
fit_val <- glmnet(mydata_normalized_tras_val, class_Tcell_val == "Positive", alpha = 1, normalize=FALSE ) # 1 for lasso
plot(fit_val) # default plotting
plot_glmnet(fit_val, xvar = "lambda") # beta > lambda are shrunken to 0

lasso_cv_tras_val <- cv.glmnet(mydata_normalized_tras_val,
                           class_Tcell_val == "Positive", 
                           alpha = 1, 
                           nfolds = 10)
plot(lasso_cv_tras_val)

#Binary classify with CIBERSORT result: presence or absence of T cell
lasso_model_val <- glmnet(mydata_normalized_tras_val,
                        class_Tcell_val == "Positive", 
                        lambda = lasso_cv_tras_val$lambda.min, 
                        alpha = 1)
lasso_coef_val <- as.matrix(coef(lasso_model_val))
nonzero_lasso_val <- lasso_coef_val[lasso_coef_val != 0,]
nonzero_lasso_val
#write.csv(nonzero_lasso_1, "nonzero_coeff_lasso_cibersort.csv")
# write.csv(nonzero_lasso_1, "test_nonzero_lasso.csv")
# length(nonzero_lasso_1)

#Linear regression with quanTIseq CD8 T cell fraction
fit2_val <- glmnet(mydata_normalized_tras_val, res_quantiseq1_val$`T cell CD8+`, alpha = 1, normalize=FALSE ) # 1 for lasso
plot(fit2_val) # default plotting
plot_glmnet(fit2_val, xvar = "lambda") # beta > lambda are shrunken to 0

set.seed(1)
lasso_cv_tras_2_val <- cv.glmnet(mydata_normalized_tras_val,
                             res_quantiseq1_val$`T cell CD8+`, 
                             alpha = 1, 
                             nfolds = 10)
plot(lasso_cv_tras_2_val)

lasso_model_2_val <- glmnet(mydata_normalized_tras_val,
                        res_quantiseq1_val$`T cell CD8+`, 
                        lambda = lasso_cv_tras_2_val$lambda.min, 
                        alpha = 1)
lasso_coef_2_val <- as.matrix(coef(lasso_model_2_val))
nonzero_lasso_2_val <- lasso_coef_2_val[lasso_coef_2_val != 0,]
nonzero_lasso_2_val
# write.csv(cbind(nonzero_lasso_2_val,names(nonzero_lasso_2_val)), "nonzero_coeff_lasso_quantiseq_val.csv", row.names = TRUE)
length(nonzero_lasso_2_val)

intersect(names(nonzero_lasso_1),names(nonzero_lasso_2))
#no intersect

###correlation with CD8
genelist_val <- c("CASP3","CCR2","CCR5","CXCR3","CCL2","CCL5","CCL7","CCL8","CCL20","CXCL9","CXCL10","CXCL11","CXCL12","XCL1","TGFB1","PDCD1","PTGES")
#genelist_val <- c("CXCL9","CXCL10","CXCR3","CCL2","CCL8")
cd8a_data_val <- c(mydata_log2_tras_val[,"CD8A"])
hist(c(mydata_log2_tras_val[,"CD8A"]),xlab="log2(CD8A) (RPKM)",main="CD8A Gene Expression")#Reads Per Kilobase per Million mapped reads
plot(cd8a_data_val,genedata_val, 
     xlab="CD8A",ylab=genelist_val,
     main=paste("validation data, cor=",round(cor(cd8a_data_val,genedata_val),3)))
#multiple correlations
p_corr_CD8_cytokine_val <- c()
for(i in 1:length(genelist_val)){
  gene_curr <- genelist_val[i]
  genedata_curr <- c(mydata_log2_tras_val[,gene_curr])
  cor_curr <- round(cor(cd8a_data_val, genedata_curr),3)
  p_corr_CD8_cytokine_val[i] <- plot(cd8a_data_val,genedata_curr, 
                                 xlab="CD8A",ylab=genelist_val[i],
                                 main=paste("validation data, cor=",cor_curr))
  
}















#******************************
#*
#use the GSE114082 dataset, which has 17 sample RNAseq before and after trastuzumab therapy
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114082
#*
#*
#*

ges_tras <- getGEO("GSE114082", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
meta_data <- pData(ges_tras)

readin_tras <- read.table("E:/OneDrive - University of Pittsburgh/Lab_OneDrive/R/TMElab_projects_MacBackup/R_Rosy/BRCAfibroblastDrugResistance/GSE114082_non-normalized.txt", sep='\t', header=T,
                          row.names = 1)
nonnormalizeddata_tras <- readin_tras[,grep("Signal",colnames(readin_tras))]
all_gene_symbol <- readin_tras["SYMBOL"]
# write.csv(all_gene_symbol,"GSE_gene_symbol.csv")
rownames(nonnormalizeddata_tras) <- make.names(all_gene_symbol$SYMBOL, unique = TRUE)
nonnormalizeddata_tras <- nonnormalizeddata_tras[!duplicated(rownames(nonnormalizeddata_tras)), ]

#another way to make names if symbols are not provided: 
# probeID <- rownames(readin_tras)
# symbol_name_to_match <- mget(x = probeID,envir = illuminaHumanv4SYMBOL,ifnotfound = NA)
# rownames(nonnormalizeddata_tras) <- make.names(
#   ifelse(symbol_name_to_match =="NA",probeID,symbol_name_to_match),
#   unique = TRUE)


###Data normalization
#TPM original data
mydata_nonnormalized_tras_t <- t(as.matrix(nonnormalizeddata_tras))
colnames(mydata_nonnormalized_tras_t) <- make.names(colnames(mydata_nonnormalized_tras_t))
rownames(mydata_nonnormalized_tras_t) <- make.names(rownames(mydata_nonnormalized_tras_t))
#Log2 TPM
mydata_log2_tras <- log2(mydata_nonnormalized_tras_t+1)
#normalize log2 data with percentage
mydata_normalized_tras <- percentize(mydata_log2_tras)%>%as.matrix

heatmap(mydata_nonnormalized_tras_t[1:20,1:20])
heatmap(mydata_log2_tras[1:20,1:20])
heatmap(mydata_normalized_tras[1:20,1:20])
{
  #OR ue the preprocessCore package to normalize data
  library(limma)
  library(preprocessCore) 
  mydata_normalized_tras<-normalize.quantiles(mydata_nonnormalized_tras_t,copy=TRUE)
  rownames(mydata_normalized_tras) <- rownames(mydata_nonnormalized_tras_t)%>%make.names()
  colnames(mydata_normalized_tras) <- colnames(mydata_nonnormalized_tras_t)%>%make.names()
  heatmap(mydata_normalized_tras[1:20,1:20])
}


###

#*************
# Deconvolution: identification of immune cell fraction
#*************
#
library(immunedeconv)
{
  res_estimate <- deconvolute(nonnormalizeddata_tras, "estimate")
  write.csv(res_estimate,"ESTIMATE_GSE114082_Immune_Stromal_Signature.csv",row.names = TRUE,col.names = TRUE)
}
res_estimate <- read.csv("ESTIMATE_GSE114082_Immune_Stromal_Signature.csv")
rownames_to_be_added <- res_estimate$cell_type
res_estimate <- t(res_estimate[,-c(1:2)])%>%as.data.frame
colnames(res_estimate) <- rownames_to_be_added


# devtools::install_github("psyteachr/introdataviz")
#create df for split violin plot
df <- data.frame(
  group = c(rep("stroma score", 34), 
            rep("immune score", 34), 
            rep("estimate score", 34), 
            rep("tumor purity", 34)),
  therapy = c(rep(c("pre","post"), 68)),
  value = c(res_estimate$`stroma score`,
            res_estimate$`immune score`,
            res_estimate$`estimate score`,
            res_estimate$`tumor purity`
            )
)

#use ggplot and self-written function to create the split violin plot
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

df %>%
  ggplot(aes(x=group,y=value,fill=therapy))+
  geom_split_violin()+
  ylab("ESTIMATE results")+
  stat_compare_means(method = "t.test", label.y=0.4,paired=TRUE)+ 
  stat_compare_means(label = "p.signif", paired=TRUE,method = "t.test", hide.ns = TRUE)

#***
#*paired box plot for stroma score, pre vs. post-therapy
df_stroma <- data.frame(
  therapy = c(rep(c("pre","post"), 17)),
  value = c(res_estimate$`stroma score`)
)

ggpaired(df_stroma, x = "therapy", y = "value",line.color = "gray", line.size = 0.4,
         color="therapy",palette = "npg")+
  ylab("ESTIMATE stroma score")+
  stat_compare_means(method = "t.test",paired=TRUE)


#***
#*box plot for post minus pre PLK1 and SERPINE1 changes, high vs. low pre-therapy stroma scores
gene <- "SERPINE1"
idx_pre <- which(meta_data$characteristics_ch1.1=="time: Pre")
idx_post <- which(meta_data$characteristics_ch1.1=="time: Post")
#in all patients that are paired:
df_deltagene_pre_stroma <- data.frame(
  gene_value = as.numeric(nonnormalizeddata_tras[gene,idx_post]-nonnormalizeddata_tras[gene,idx_pre]),#post minus pre therapy gene expression levels
  stroma_value_pre = as.numeric(res_estimate$`stroma score`[idx_pre]),
  stroma_value_post = as.numeric(res_estimate$`stroma score`[idx_post])
)

stroma_value_pre_class <- ifelse(df_deltagene_pre_stroma$stroma_value_pre>median(df_deltagene_pre_stroma$stroma_value_pre),"high","low")
stroma_value_post_class <- ifelse(df_deltagene_pre_stroma$stroma_value_post>median(df_deltagene_pre_stroma$stroma_value_post),"high","low")
df_deltagene_pre_stroma <- cbind(df_deltagene_pre_stroma, stroma_value_pre_class, stroma_value_post_class)

df_deltagene_pre_stroma %>%
  ggplot( aes(x=stroma_value_pre_class, y=gene_value, fill=stroma_value_pre_class)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=stroma_value_pre_class))+
  ylab(paste("Post minus pre therapy",gene))+
  stat_compare_means(method="wilcox.test")


t.test(df_deltagene_pre_stroma$gene_value[which(df_deltagene_pre_stroma$stroma_value_pre_class=="low")], mu = 0)

#***
#*Correlation plot for pre minus post therapy PLK1 vs. stroma scores
ggplot(df_deltagene_pre_stroma, aes(x=stroma_value_pre, y=gene_value)) +
  geom_point()+ 
  geom_smooth(method = "lm", fill = NA)+
  stat_cor()+
  ylab(paste("Post minus pre therapy",gene))


#***
#*box plot for pre or post PLK1 and SERPINE1 expression, high vs. low pre or post-therapy stroma scores
gene <- "SERPINE1"
idx_pre <- which(meta_data$characteristics_ch1.1=="time: Pre")
idx_post <- which(meta_data$characteristics_ch1.1=="time: Post")
#in all patients that are paired:
df_gene_pre_post_stromaclass <- data.frame(
  gene_value_pre = as.numeric(nonnormalizeddata_tras[gene,idx_pre]),#pre therapy gene expression levels
  gene_value_post = as.numeric(nonnormalizeddata_tras[gene,idx_post]),#post therapy gene expression levels
  stroma_value_pre = as.numeric(res_estimate$`stroma score`[idx_pre]),
  stroma_value_post = as.numeric(res_estimate$`stroma score`[idx_post])
)

stroma_value_pre_class <- ifelse(df_gene_pre_post_stromaclass$stroma_value_pre>median(df_gene_pre_post_stromaclass$stroma_value_pre),"high","low")
stroma_value_post_class <- ifelse(df_gene_pre_post_stromaclass$stroma_value_post>median(df_gene_pre_post_stromaclass$stroma_value_post),"high","low")
df_deltagene_pre_stroma <- cbind(df_gene_pre_post_stromaclass, stroma_value_pre_class, stroma_value_post_class)

df_gene_pre_post_stromaclass %>%
  ggplot( aes(x=stroma_value_pre_class, y=gene_value_pre, fill=stroma_value_pre_class)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=stroma_value_pre_class))+
  ylab(paste("Pre therapy",gene))+
  stat_compare_means(method="wilcox.test")

df_gene_pre_post_stromaclass %>%
  ggplot( aes(x=stroma_value_pre_class, y=gene_value_post, fill=stroma_value_pre_class)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=stroma_value_pre_class))+
  ylab(paste("Post therapy",gene))+
  stat_compare_means(method="wilcox.test")



#***
#*paired box plot and correlations
"PLK1" %in% rownames(nonnormalizeddata_tras)
"SERPINE1" %in% rownames(nonnormalizeddata_tras)
gene <- "SERPINE1"
df_gene <- data.frame(
  therapy = c(rep(c("pre","post"), 17)),
  value = as.numeric(log2(nonnormalizeddata_tras[gene,]))
)

ggpaired(df_gene, x = "therapy", y = "value",line.color = "gray", line.size = 0.4,
         color="therapy",palette = "npg")+
  ylab(paste("log2",gene))+
  stat_compare_means(method = "t.test",paired=TRUE)


#***
#*correlations between pre- and post-therapy gene expression and stroma scores
gene <- "SERPINE1"
library(forcats)
df_gene_stroma <- data.frame(
  therapy = c(rep(c("pre","post"), 17)),
  gene_value = as.numeric(log2(nonnormalizeddata_tras[gene,])),
  stroma_value = c(res_estimate$`stroma score`)
)%>%mutate(therapy = fct_relevel(therapy, 
                            "pre","post")) #mutate forces the levels of therapy to follow a "pre" "post" order
ggplot(df_gene_stroma, aes(x=gene_value, y=stroma_value,color=therapy)) +
  geom_point()+ 
  geom_smooth(method = "lm", fill = NA)+
  stat_cor()+
  xlab(paste("log2",gene))+
  ylab("ESTIMATE stroma score")


















# #genelist 6.1: genes correlated to M1; M2; M2/total changes of pre- post-trastuzumab
# pval <- c()
# corr <- c()
# row_odd <- seq_len(nrow(mydata_normalized_tras)) %% 2    
# groups <- unique(df$group)
# mydata_normalized_diff <- mydata_nonnormalized_tras_t[!row_odd,]-mydata_nonnormalized_tras_t[row_odd,]
# rownames(mydata_normalized_diff) <- gsub("Post.AVG_Signal","Diff",rownames(mydata_normalized_diff))
# #each row is a patient and each column is a gene
# library(doParallel)
# library(foreach)
# detectCores()
# cl <- makeCluster(7)
# registerDoParallel(cl)
# 
# for (i in 1:length(groups)){
#   curr_group <- groups[i]
#   fractiondiff <- df$value[!row_odd & df$group==curr_group]-df$value[row_odd & df$group==curr_group]
#   corr_coeff <- c()
#   # for (j in 1:ncol(mydata_normalized_diff)){
#   #   corr_coeff[j] <- cor(fractiondiff,mydata_normalized_diff[,j])
#   # }
#   corr_coeff <- foreach(j=1:ncol(mydata_normalized_diff), .combine=cbind) %dopar% {
#     temp_corr = cor(fractiondiff,mydata_normalized_diff[,j])
#     temp_corr
#     # temp_corr  #Equivalent to corr_coeff = cbind(corr_coeff, temp_corr_coeff)
#   }
#   print(curr_group)
#   print(sum(abs(corr_coeff)>0.5))
#   idx_sig <- which(abs(corr_coeff)>0.5)
#   sig_df <- data.frame(geneName = colnames(mydata_normalized_diff)[idx_sig],
#                        spearman = corr_coeff[idx_sig])
#   sig_KEGG <- sig_df[sig_df$geneName==intersect(geneset_chemokine_cytokine,colnames(mydata_normalized_diff)[idx_sig]),]
#   write.csv(sig_df,paste(curr_group,"all genes Pearson.csv"))
#   write.csv(sig_KEGG,paste(curr_group,"KEGG chemokine receptor Pearson.csv"))
# }
# 
# 
# 
# 
# 
# 
# 
# ###QuanTIseq
# res_quantiseq <- deconvolute(nonnormalizeddata_tras, "quantiseq", tumor = TRUE)
# res_quantiseq1 <- t(res_quantiseq[,-1])%>%as.data.frame
# colnames(res_quantiseq1) <- res_quantiseq$cell_type
# rownames(res_quantiseq1) <- gsub(".AVG_Signal","",rownames(res_quantiseq1))
# colnames(res_quantiseq) <- gsub(".AVG_Signal","",colnames(res_quantiseq))
# 
# hist(res_quantiseq1$`T cell CD8+`,breaks = 40,xlab="quanTIseq CD8 Test",main="")
# 
# res_quantiseq %>%
#   gather(sample, fraction, -cell_type) %>%
#   # plot as stacked bar chart
#   ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
#   geom_bar(stat='identity') +
#   coord_flip() +
#   scale_fill_brewer(palette="Paired") +
#   scale_x_discrete(limits = rev(levels(res_quantiseq)))
# 
# 
# ###Classify
# #classify by pre-treatment T cell presence: positive/negative (results from CIBERSORT)
# idx_Tpositive <- which(res_cibersort1$`T cell CD8+`!=0)
# idx_Tnegative <- which(res_cibersort1$`T cell CD8+`==0)
# class_Tcell <- ifelse(res_cibersort1$`T cell CD8+`!=0,"Positive","Negative")
# 
# #classify by trastuzumab treatment condition: before/after
# idx_pre <- grep("Pre",rownames(mydata_normalized_tras))
# idx_post <- grep("Post",rownames(mydata_normalized_tras))
# class_therapy <- ifelse(grepl("Pre",rownames(mydata_normalized_tras)),"Pre","Post")
# 
# table(class_Tcell,class_therapy)
# 
# #classify by post-treatment T cell level: increased/decreased comparing to pre-treatment
# #(results from quanTIseq) 
# idx_increase <- which(res_quantiseq1$`T cell CD8+`[idx_pre]>res_quantiseq1$`T cell CD8+`[idx_post])
# idx_decrease <- which(res_quantiseq1$`T cell CD8+`[idx_pre]<=res_quantiseq1$`T cell CD8+`[idx_post])
# class_increase <- ifelse(c(1:length(idx_pre))%in%idx_increase,"Increase","Decrease")
# 
# 
# #DEgenes between pre and post samples
# pval <- c()
# fc <- c()
# row_odd <- seq_len(nrow(mydata_normalized_tras)) %% 2    
# for (i in 1:ncol(mydata_normalized_tras)){
#   pre_data <- mydata_normalized_tras[row_odd==1,i]
#   post_data <- mydata_normalized_tras[row_odd==0,i]
#   fc[i] <- sum(post_data)/sum(pre_data)-1
#   pval[i] <- t.test(pre_data, post_data, paired = TRUE, alternative = "two.sided")$p.val
# }
# 
# padj <- p.adjust(pval)
# 
# #threashold DEGs between pre and post samples and save as .csv
# idx_deg <- which(abs(fc)>0.5 & pval<0.0001)
# DEG_pre_post_trastuzumab_GSE114082 <- data.frame(
#   geneName = colnames(mydata_normalized_tras)[idx_deg],
#   p = pval[idx_deg],
#   fc = fc[idx_deg]
# )
# write.csv(DEG_pre_post_trastuzumab_GSE114082,"DEG_pre_post_trastuzumab_GSE114082.csv")
# 
# #which are hypoxia markers
# DEG_pre_post_trastuzumab_GSE114082 <- read.csv("~/Desktop/R_Rosy/BRCAfibroblastDrugResistance/DEG_pre_post_trastuzumab_GSE114082.csv",header=TRUE)[,c(-1)]
# intersect(DEG_pre_post_trastuzumab_GSE114082$geneName,hypoxia_genelist)
# 
# #which are ligand/receptors
# lig_rec_genelist <- read.csv("ligand_receptor_genelist.csv",header = TRUE)[2]%>%as.matrix
# intersect(lig_rec_genelist,DEG_pre_post_trastuzumab_GSE114082$geneName)
# write.csv(DEG_pre_post_trastuzumab_GSE114082[intersect(lig_rec_genelist,DEG_pre_post_trastuzumab_GSE114082$geneName),],
#           "DEG_pre_post_trastuzumab_GSE114082_ligand_receptor.csv"
# )
# geneset_chemokine_cytokine <- read.table("~/Desktop/R_Rosy/OVCAmicrophageInfiltration/geneset_KEGG_cytokine_cytokine_receptor_interaction.txt",
#                                          sep="\t",
#                                          header=TRUE)[-1,]
# intersect(geneset_chemokine_cytokine,DEG_pre_post_trastuzumab_GSE114082$geneName)
# write.csv(DEG_pre_post_trastuzumab_GSE114082[intersect(geneset_chemokine_cytokine,DEG_pre_post_trastuzumab_GSE114082$geneName),],
#           "DEG_pre_post_trastuzumab_GSE114082_KEGG_chemokine.csv"
# )
# #further narrow down by threasholding adjusted p values
# colnames(mydata_normalized_tras)[which(padj<0.05)]
# #intersection with the BT474 pre post trastuzumab data
# geneset_chemokine_cytokine <- read.table("~/Desktop/R_Rosy/OVCAmicrophageInfiltration/geneset_KEGG_cytokine_cytokine_receptor_interaction.txt",
#                                          sep="\t",
#                                          header=TRUE)[-1,]
# lig_rec_genelist <- read.csv("ligand_receptor_genelist.csv",header = TRUE)[2]%>%as.matrix
# 
# intersect(geneset_chemokine_cytokine,lig_rec_genelist)
# 
# genesBT474 <- read.csv("~/Desktop/R_Rosy/CD8TrastuzumabHer2BreastCancer/BT474 DEGs pre- post-trastuzumab.csv",
#                        sep=",")
# intersectBT474 <- which(genesBT474$GeneName%in%DEG_pre_post_trastuzumab_GSE114082$geneName)
# intersectBT474 <- intersect(genesBT474$GeneName,DEG_pre_post_trastuzumab_GSE114082$geneName)[
#   which(genesBT474[genesBT474$GeneName%in%intersectBT474,4]*DEG_pre_post_trastuzumab_GSE114082[DEG_pre_post_trastuzumab_GSE114082$geneName%in%intersectBT474,3]>0)
# ]
# 
# 
# 
# 
# ###Correlations
# cor1 <- cor(res_quantiseq1$`T cell CD8+`,res_cibersort1$`T cell CD8+`)
# plot(res_quantiseq1$`T cell CD8+`,res_cibersort1$`T cell CD8+`,main=paste("cor=",round(cor1,3)),xlab="quanTIseq CD8 test",ylab="CIBERSORT CD8 test")
# 
# #genelist <- c("TGFB2")
# genelist <- c("CASP3")
# 
# genedata <- c(mydata_normalized_tras[,genelist])
# cor3 <- round(cor(res_quantiseq1$`T cell CD8+`, genedata),3)
# plot(res_quantiseq1$`T cell CD8+`,genedata, 
#      xlab="T cell fraction quanTIseq",ylab=genelist,
#      main=paste("test data, cor=",cor3))
# 
# mypar(3,3)
# ###correlation with CD8
# genelist <- c("CASP3","CCR2","CCR5","CXCR3","CCL2","CCL5","CCL7","CCL8","CCL20","CXCL9","CXCL10","CXCL11","CXCL12","XCL1","TGFB1","PDCD1","PTGES")
# genelist <- c("CXCL9","CXCL10","CXCR3","CCL2","CCL8","CCR5")
# cd8a_data <- c(mydata_log2_tras[,"CD8A"])
# hist(c(mydata_log2_tras[,"CD8A"]),xlab="log2(CD8A), RPKM",main="CD8A Gene Expression")#Reads Per Kilobase per Million mapped reads
# plot(cd8a_data,genedata, 
#      xlab="CD8A",ylab=genelist,
#      main=paste("test data, cor=",round(cor(cd8a_data,genedata),3)))
# #multiple correlations
# mypar(3,3)
# for(i in 1:length(genelist)){
#   gene_curr <- genelist[i]
#   genedata_curr <- c(mydata_log2_tras[,gene_curr])
#   cor_curr <- round(cor(cd8a_data, genedata_curr),3)
#   p_corr_CD8_cytokine[i] <- plot(cd8a_data,genedata_curr, 
#                                  xlab="CD8A",ylab=genelist[i],
#                                  main=paste("cor =",cor_curr))
#   abline(lm(genedata_curr ~ cd8a_data))
# }
# mypar(1,1)
# 
# 
# 
# # cor_positive <- round(cor(res_quantiseq1$`T cell CD8+`[idx_Tpositive],genedata[idx_Tpositive]),3)
# # cor_negative <- round(cor(res_quantiseq1$`T cell CD8+`[idx_Tnegative],genedata[idx_Tnegative]),3)
# # 
# # plot(res_quantiseq1$`T cell CD8+`[idx_Tpositive],genedata[idx_Tpositive], 
# #      xlab="T cell fraction quanTIseq",ylab=genelist,
# #      main=paste("Test T positive, cor=",cor_positive))
# # 
# # plot(res_quantiseq1$`T cell CD8+`[idx_Tnegative],genedata[idx_Tnegative], 
# #      xlab="T cell fraction quanTIseq",ylab=genelist,
# #      main=paste("Test T positive, cor=",cor_negative))
# 
# cor_pre <- round(cor(res_quantiseq1$`T cell CD8+`[idx_pre],genedata[idx_pre]),3)
# cor_post <- round(cor(res_quantiseq1$`T cell CD8+`[idx_post],genedata[idx_post]),3)
# 
# cor_inc <- round(cor(res_quantiseq1$`T cell CD8+`[2*idx_increase],genedata[2*idx_increase]),3)
# cor_dec <- round(cor(res_quantiseq1$`T cell CD8+`[2*idx_decrease],genedata[2*idx_decrease]),3)
# 
# plot(res_quantiseq1$`T cell CD8+`[2*idx_increase],genedata[2*idx_increase],
#      xlab="T cell fraction quanTIseq",ylab=genelist,
#      main=paste("Test T increase, cor=",cor_inc))
# 
# plot(res_quantiseq1$`T cell CD8+`[2*idx_decrease],genedata[2*idx_decrease],
#      xlab="T cell fraction quanTIseq",ylab=genelist,
#      main=paste("Test T decrease, cor=",cor_dec))
# 
# 
# ###Organize data
# #mydata_T_forLasso <- data.frame(condition=as.factor(class_Tcell), mydata_normalized_tras)
# mydata_T_forLasso <- data.frame(condition=as.factor(class_increase), mydata_nonnormalized_tras_t[c(idx_pre,idx_post),])
# 
# 
# 
# #*******************************
# #* Plot figures for cytokines
# #*******************************
# #*
# #check cytokine expression
# boxplot_genelist <- c()
# boxplot_genelist%in%colnames(mydata_T_forLasso)
# 
# ###Boxplots for POST-TREATMENT cytokine levels
# # ONLY POST-treatment data (row 18:34) are used!!!
# data_cytokine_boxplot <- c()
# for(i in 1:length(boxplot_genelist)){
#   data_cytokine_boxplot <- rbind(data_cytokine_boxplot , 
#                                  data.frame(
#                                    Label=class_increase,#Label is the post-treatment T cell level.
#                                    variable=c(replicate(17,boxplot_genelist[i])),
#                                    value=mydata_T_forLasso[18:34,boxplot_genelist[i]])
#                                  )
# }
# 
# #horizontal version
# p <- ggplot(data = data_cytokine_boxplot, aes(x=variable, y=value)) + 
#   geom_boxplot(aes(fill=Label))+
#   scale_color_manual(values = c("#00AFBB", "#E7B800")) +
#   stat_compare_means(aes(group = Label), label = "p.format")#p.format,p.signif
# p 
# 
# #facet version
# p + facet_wrap( ~ variable, scales="free") #facet version
# 
# 
# # ggboxplot(mydata_T_forLasso[18:34,], x="condition", y="CCL5", color="condition") +
# #   geom_signif(comparisons = list(c("Increase","Decrease")), 
# #               map_signif_level=FALSE)+
# #   ggtitle("CCL5")
# # 
# 
# ggplot(data = mydata_T_forLasso, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Label))
# 
# 
# ggboxplot(mydata_T_forLasso, x="condition", y="PTGES", color="condition")+
#   geom_signif(comparisons = list(c("Positive","Negative")), 
#               map_signif_level=FALSE)
# 
# ggboxplot(mydata_T_forLasso, x="condition", y="CCL5", color="condition") +
#   geom_signif(comparisons = list(c("Positive","Negative")), 
#               map_signif_level=FALSE)+
#   ggtitle("CCL5")
# 
# ggboxplot(mydata_T_forLasso, x="condition", y="CXCL12", color="condition") +
#   geom_signif(comparisons = list(c("Positive","Negative")), 
#               map_signif_level=FALSE)+
#   ggtitle("CXCL12")
# 
# ggboxplot(mydata_T_forLasso, x="condition", y="CXCL10", color="condition") +
#   geom_signif(comparisons = list(c("Positive","Negative")), 
#               map_signif_level=FALSE)+
#   ggtitle("CXCL10")
# 
# ggboxplot(mydata_T_forLasso, x="condition", y="CXCL11", color="condition") +
#   geom_signif(comparisons = list(c("Positive","Negative")), 
#               map_signif_level=FALSE)+
#   ggtitle("CXCL11")
# 
# ggboxplot(mydata_T_forLasso, x="condition", y="CXCR3", color="condition") +
#   geom_signif(comparisons = list(c("Positive","Negative")), 
#               map_signif_level=FALSE)+
#   ggtitle("CXCR3")
# 
# 
# ###Paired boxplot with significance level
# # paired_data <- cbind(res_quantiseq1, class_therapy)
# paired_data <- cbind(res_cibersort1, class_therapy)
# ggpaired(paired_data, x ="class_therapy", y = "T cell CD8+",
#          color = "class_therapy", line.color = "gray", line.size = 0.4,
#          palette = "npg")+
#   geom_signif(comparisons = list(c("Pre","Post")), 
#               map_signif_level=FALSE)
# 
#  
# ###
# 
# #*************
# # Feature Selection
# #*************
# #
# ###Lasso regression for parameter selection
# fit <- glmnet(mydata_normalized_tras, class_Tcell == "Positive", alpha = 1, normalize=FALSE ) # 1 for lasso
# plot(fit) # default plotting
# plot_glmnet(fit, xvar = "lambda") # beta > lambda are shrunken to 0
# 
# lasso_cv_tras <- cv.glmnet(mydata_normalized_tras,
#                            class_Tcell == "Positive", 
#                           alpha = 1, 
#                           nfolds = 10)
# plot(lasso_cv_tras)
# 
# #Binary classify with CIBERSORT result: presence or absence of T cell
# lasso_model_1 <- glmnet(mydata_normalized_tras,
#                         class_Tcell == "Positive", 
#                             lambda = lasso_cv_tras$lambda.min, 
#                             alpha = 1)
# lasso_coef_1 <- as.matrix(coef(lasso_model_1))
# nonzero_lasso_1 <- lasso_coef_1[lasso_coef_1 != 0,]
# nonzero_lasso_1
# #write.csv(nonzero_lasso_1, "nonzero_coeff_lasso_cibersort.csv")
# write.csv(nonzero_lasso_1, "test_nonzero_lasso.csv")
# length(nonzero_lasso_1)
# 
# #Linear regression with quanTIseq CD8 T cell fraction
# fit2 <- glmnet(mydata_normalized_tras, res_quantiseq1$`T cell CD8+`, alpha = 1, normalize=FALSE ) # 1 for lasso
# plot(fit2) # default plotting
# plot_glmnet(fit2, xvar = "lambda") # beta > lambda are shrunken to 0
# 
# lasso_cv_tras_2 <- cv.glmnet(mydata_normalized_tras,
#                            res_quantiseq1$`T cell CD8+`, 
#                            alpha = 1, 
#                            nfolds = 10)
# plot(lasso_cv_tras_2)
# 
# lasso_model_2 <- glmnet(mydata_normalized_tras,
#                         res_quantiseq1$`T cell CD8+`, 
#                         lambda = lasso_cv_tras_2$lambda.min, 
#                         alpha = 1)
# lasso_coef_2 <- as.matrix(coef(lasso_model_2))
# nonzero_lasso_2 <- lasso_coef_2[lasso_coef_2 != 0,]
# nonzero_lasso_2
# write.csv(nonzero_lasso_2, "nonzero_coeff_lasso_quantiseq.csv")
# length(nonzero_lasso_2)
# 
# intersect(names(nonzero_lasso_1),names(nonzero_lasso_2))
# #no intersect
# 
# 
# ###
# 
# #*****************
# # DEG analysis and Enrichment analysis
# #*****************
# #*
# ###Identify Differentially expressed genes
# 
# #CD8 T positive vs. negative
# #fold change of normalized data
# #@Rosy need to check if using raw non-normalized TPM data is different
# fc_all <- c()
# pval_all <- c()
# for (i in 1:ncol(mydata_normalized_tras)){
#   fc_all[i] <- sum(mydata_normalized_tras[idx_Tpositive,i])/sum(mydata_normalized_tras[idx_Tnegative,i]) -1
#   #pval_all[i] <- t.test(mydata_normalized_tras[idx_Tpositive,i],mydata_normalized_tras[idx_Tnegative,i])$p.val
#   }
# #adj.p_all <- p.adjust(pval_all,method = "BH") #default method is bonferroni
# 
# 
# #upregulated gene ratio
# mean(fc_all>0)
# #88% are upregulated
# 
# #define DEG by pval and fc threshold
# # fc_threshold = 0.5
# #length(which(abs(fc_all)>fc_threshold))#4282
# 
# #define DEG by abs foldchange
# length(which(abs(fc_all)>0.3))
# deg_idx <- order(abs(fc_all))[1:8000]
# deg_list <- colnames(mydata_normalized_tras)[deg_idx]
# 
# res <- gost(deg_list, organism="hsapiens")
# gostplot(res,interactive=TRUE)
# 
# 
# #************************
# #*Ligand-receptor analysis
# #*ligand-receptor list form this paper:
# #*https://www.nature.com/articles/ncomms8866
# #*data: validation set
# 
# lig_rec_genelist <- read.csv("ligand_receptor_genelist.csv",header = TRUE)[2]%>%as.matrix
# 
# #correlate every gene to CD8A, CD8B and CD8+ T cell fraction(quanTIseq)
# intersect(lig_rec_genelist,colnames(mydata_log2_tras_val))
# intersect_lig_rec_test <- c(0)
# for (i in 1:length(lig_rec_genelist)){
#   if(lig_rec_genelist[i]%in%colnames(mydata_log2_tras_val)){
#     intersect_lig_rec_test[i] <- 1
#   }
# }
# new_idx <- which(is.na(intersect_lig_rec_test)==TRUE)
# 
# new_lig_rec_genelist <- lig_rec_genelist[-new_idx]
# 
# cor_CD8A <- c()
# for(i in 1:length(new_lig_rec_genelist)){
#   cor_CD8A[i] <- cor(mydata_log2_tras_val[,new_lig_rec_genelist[i]],mydata_log2_tras_val[,"CD8A"])
# }
# names(cor_CD8A)<-new_lig_rec_genelist
# hist(cor_CD8A)
# hi_cor_CD8A <- cor_CD8A[which(abs(cor_CD8A)>0.5)]
# 
# cor_CD8B <- c()
# for(i in 1:length(new_lig_rec_genelist)){
#   cor_CD8B[i] <- cor(mydata_log2_tras_val[,new_lig_rec_genelist[i]],mydata_log2_tras_val[,"CD8B"])
# }
# names(cor_CD8B)<-new_lig_rec_genelist
# hist(cor_CD8B)
# hi_cor_CD8B <- cor_CD8B[which(abs(cor_CD8B)>0.5)]
# 
# #screen the cytokines that have high correlation (>0.5) to both CD8A and CD8B
# hi_cor_CD8_both <- cor_CD8B[intersect(names(hi_cor_CD8A),names(hi_cor_CD8B))]+
#   cor_CD8A[intersect(names(hi_cor_CD8A),names(hi_cor_CD8B))]
# 
# par(cex=0.7,mar=c(5,2,2,0))
# barplot(hi_cor_CD8_both[order(hi_cor_CD8_both,decreasing = TRUE)],las=2,
#         main="Corr (CD8A, lig-rec) + corr (CD8B, lig-rec)")
# 
# cor(mydata_log2_tras_val[,"CD8A"],mydata_log2_tras_val[,"CD8B"])
# plot(mydata_log2_tras_val[,"CD8A"],mydata_log2_tras_val[,"CD8B"], 
#      xlab="CD8A",ylab="CD8B")
# 
# gene <- "CCL5"
# 
# plot(mydata_log2_tras_val[,"CD8A"],mydata_log2_tras_val[,"FASLG"],
#      xlab="CD8A",ylab="FASLG")
# abline(lm(mydata_log2_tras_val[,"FASLG"]~mydata_log2_tras_val[,"CD8A"]))
# cor(mydata_log2_tras_val[,"CD8A"],mydata_log2_tras_val[,"FASLG"])
#      
# plot(mydata_log2_tras_val[,"CD8A"],mydata_log2_tras_val[,gene],
#      xlab="CD8A",ylab=gene)
# abline(lm(mydata_log2_tras_val[,gene]~mydata_log2_tras_val[,"CD8A"]))
# cor(mydata_log2_tras_val[,"CD8A"],mydata_log2_tras_val[,gene])



