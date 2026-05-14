#***************
#*Comparing residual disease vs. complete pathological response (pCR) 

# GSE123845
library(GEOquery)

# gse <- getGEO("GSE123845")
# head(Meta(gse))
# gsedata <- gse[[1]]
# data_RNA <- exprs(gsedata)#no expression data
# data_clinical <- pData(gsedata)
# write.csv(data_clinical,"GSE123845_clinical_annotation.csv")

data_clinical <- read.csv("GSE123845_clinical_annotation.csv",header = TRUE, row.names = 1)

# eList <- getGEOSuppFiles("GSE123845")
# gunzip("GSE123845/GSE123845_exp_tpm_matrix.csv.gz")
data_RNA <- read.csv("GSE123845/GSE123845_exp_tpm_matrix.csv",row.names = 1,header=TRUE)

data_clinical$title==colnames(data_RNA)#TRUE for all 227 samples

table(data_clinical$pcr_status.ch1)#68 RD and 159 pCR
table(data_clinical$pcr_status.ch1, data_clinical$her2_status_diagnosis.ch1)#35 HER2 patients are RD and 48 are pCR

#select HER2 patients only
table(data_clinical$er_status_diagnosis.ch1)
HER2_id <- which(data_clinical$er_status_diagnosis.ch1==1)#102 samples are HER2+
data_RNA_HER2 <- data_RNA[,HER2_id]
data_clinical_HER2 <- data_clinical[HER2_id,]

#separate samples into pre and post therapy. Only proceed pre therapy data
grepl("^.+(_1)$",data_clinical_HER2$title)
data_RNA_HER2_pre <- data_RNA_HER2[,grepl("^.+(_1)$",data_clinical_HER2$title)]
data_clinical_HER2_pre <- data_clinical_HER2[grepl("^.+(_1)$",data_clinical_HER2$title),]


curr_gene <- "PLK1"
curr_gene %in% row.names(data_RNA_HER2_pre)
curr_gene_expr <- data_RNA_HER2_pre[which(row.names(data_RNA_HER2_pre)==curr_gene),]%>%as.numeric
hist(curr_gene_expr,breaks=100)
median_data = median(curr_gene_expr)
# idx_high <- which(curr_gene_expr>quantile(curr_gene_expr)[4])
# idx_low <- which(curr_gene_expr<=quantile(curr_gene_expr)[2])
idx_high <- which(curr_gene_expr>median_data)
idx_low <- which(curr_gene_expr<=median_data)
# idx_high <- which(curr_gene_expr>(-0.1))
# idx_low <- which(curr_gene_expr<=(-0.1))


#all samples with all 3 treatment types togetner
data_clinical_HER2_pre$pcr_status.ch1[which(data_clinical_HER2_pre$pcr_status.ch1==0)] <- "RD"
data_clinical_HER2_pre$pcr_status.ch1[which(data_clinical_HER2_pre$pcr_status.ch1==1)] <- "pCR"

id <- c(idx_high,idx_low)
gene_expr_plotting <- curr_gene_expr[id]
response_class <- data_clinical_HER2_pre$pcr_status.ch1[id]
gene_expr_group <- c(replicate(length(idx_high),"high"),replicate(length(idx_low),"low"))
table(response_class,gene_expr_group)

plotting_df <- data.frame(gene_expr_plotting,response_class)

library(hrbrthemes)
plotting_df %>%
  ggplot( aes(x=response_class, y=gene_expr_plotting, fill=response_class)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(curr_gene) +
  xlab("")+
  stat_compare_means(method="wilcox.test")













































