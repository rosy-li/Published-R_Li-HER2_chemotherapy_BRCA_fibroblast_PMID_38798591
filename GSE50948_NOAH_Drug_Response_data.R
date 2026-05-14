#***************
#*Comparing residual disease vs. complete pathological response (pCR) 

# GSE50948
library(GEOquery)

# gse <- getGEO("GSE50948")
# head(Meta(gse))
# gsedata <- gse[[1]]
# data_RNA <- exprs(gsedata)
# data_clinical <- pData(gsedata)
# write.csv(data_clinical,"GSE50948_clinical_annotation.csv")
# write.csv(data_RNA, "GSE50948_data_RNA_normalized.csv")

data_clinical <- read.csv("GSE50948_clinical_annotation.csv",header = TRUE, row.names = 1, fill=TRUE)
library(data.table)
data_RNA <- as.matrix(fread("GSE50948_data_RNA_normalized.csv"),rownames=1)%>%as.data.frame

#convert probe names to gene names
library(readr)
probe_list <- read_tsv("GPL570-55999.txt")
#* Source of the txt file: After downloading the platform probe to gene name conversion from NCBI ,
#* https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570
#* unzip. open the txt file AND remove the following rows at the beginning!!
#ID = Affymetrix Probe Set ID
#GB_ACC = GenBank Accession Number
#SPOT_ID = identifies controls
#Species Scientific Name = The genus and species of the organism represented by the probe set.
#Annotation Date = The date that the annotations for this probe array were last updated. It will generally be earlier than the date when the annotations were posted on the Affymetrix web site.
#Sequence Type = 
#Sequence Source = The database from which the sequence used to design this probe set was taken.
#Target Description = 
#Representative Public ID = The accession number of a representative sequence. Note that for consensus-based probe sets, the representative sequence is only one of several sequences (sequence sub-clusters) used to build the consensus sequence and it is not directly used to derive the probe sequences. The representative sequence is chosen during array design as a sequence that is best associated with the transcribed region being interrogated by the probe set. Refer to the "Sequence Source" field to determine the database used.
#Gene Title = Title of Gene represented by the probe set.
#Gene Symbol = A gene symbol, when one is available (from UniGene).
#ENTREZ_GENE_ID = Entrez Gene Database UID
#RefSeq Transcript ID = References to multiple sequences in RefSeq. The field contains the ID and Description for each entry, and there can be multiple entries per ProbeSet.
#Gene Ontology Biological Process = Gene Ontology Consortium Biological Process derived from LocusLink.  Each annotation consists of three parts: "Accession Number // Description // Evidence". The description corresponds directly to the GO ID. The evidence can be "direct", or "extended".
#Gene Ontology Cellular Component = Gene Ontology Consortium Cellular Component derived from LocusLink.  Each annotation consists of three parts: "Accession Number // Description // Evidence". The description corresponds directly to the GO ID. The evidence can be "direct", or "extended".
#Gene Ontology Molecular Function = Gene Ontology Consortium Molecular Function derived from LocusLink.  Each annotation consists of three parts: "Accession Number // Description // Evidence". The description corresponds directly to the GO ID. The evidence can be "direct", or "extended".

which(probe_list$ID!=rownames(data_RNA))#everything except 54629:54639 are the same
probe_list$ID[54628:54640]
rownames(data_RNA)[54628:54640]
#probes for genes are similar. Just 3 and 5 prime sequence differences

probe_list$`Gene Symbol`[1:500]#some gene symbols have the "//", remove everything after the first space

rownames(data_RNA) <- make.names(sub(" .*", "", probe_list$`Gene Symbol`),unique = TRUE)
"PLK1" %in% rownames(data_RNA)
"SERPINE1" %in% rownames(data_RNA)
"PAI1" %in% rownames(data_RNA)

#***
#*Run ESTIMATE--updated 03.19.24
#*
{
  {print("Run the following 2 lines for a new immunedeconvolution. Otherwise, open saved csv.")
    res_estimate<-immunedeconv::deconvolute(data_RNA, "estimate")
    write.csv(res_estimate,"ESTIMATE_GSE50948_NOAH_Immune_Stromal_Signature.csv",row.names = TRUE,col.names = TRUE)
  }
  ### read ESTIMATE results
  res_estimate<-read.csv("ESTIMATE_GSE50948_NOAH_Immune_Stromal_Signature.csv")
  rownames(res_estimate) <- res_estimate$cell_type
  res_estimate <- res_estimate[,-c(1,2)]
}


#***
#*select HER2 patients only
#*
HER2_id <- which(data_clinical$her2.ch1=="HER2+")
data_clinical_HER2 <- data_clinical[HER2_id,]
data_RNA_HER2 <- data_RNA[,HER2_id]
res_estimate_HER2 <- res_estimate[,HER2_id]


#Correlations between HER2 stroma scores and manually annotated stroma cell numbers
plot(res_estimate_HER2["stroma score",]%>%as.numeric, data_clinical_HER2$stroma_cells.....ch1)

#***
#*select HER2+ AND high stroma content patients only
#*MUST run the select HER2 patients only block first
#select top 30 patients or top 50%, 25% etc. with the highest ESTIMATE stroma score
hist(res_estimate_HER2["stroma score",]%>%as.numeric)
# idx_highstroma <- which(res_estimate_HER2["stroma score",]%>%as.numeric>mean(res_estimate_HER2["stroma score",]%>%as.numeric))
idx_highstroma <- sort(res_estimate_HER2["stroma score",]%>%as.numeric, index.return=TRUE, decreasing=TRUE)[[2]][1:20]

#OR select top 30 patients with high manually annotated stroma cell numbers
idx_highstroma <- which(data_clinical_HER2$stroma_cells.....ch1>20)

data_clinical_HER2 <- data_clinical_HER2[idx_highstroma,]
data_RNA_HER2 <- data_RNA_HER2[,idx_highstroma]
res_estimate_HER2 <- res_estimate_HER2[,idx_highstroma]


#all samples are pre-treatment
curr_gene <- "SERPINE1"
curr_gene %in% row.names(data_RNA_HER2)
curr_gene_expr <- data_RNA_HER2[which(row.names(data_RNA_HER2)==curr_gene),]%>%as.numeric
hist(curr_gene_expr,breaks=100)
median_data = median(curr_gene_expr)
# idx_high <- which(curr_gene_expr>quantile(curr_gene_expr)[4])
# idx_low <- which(curr_gene_expr<=quantile(curr_gene_expr)[2])
idx_high <- which(curr_gene_expr>median_data)
idx_low <- which(curr_gene_expr<=median_data)
# idx_high <- which(curr_gene_expr>(-0.1))
# idx_low <- which(curr_gene_expr<=(-0.1))


#all samples with all 3 treatment types togetner
id <- c(idx_high,idx_low)
gene_expr_plotting <- curr_gene_expr[id]
response_class <- data_clinical_HER2$pcr.ch1[id]
gene_expr_group <- c(replicate(length(idx_high),"high"),replicate(length(idx_low),"low"))
treatment_group <- data_clinical_HER2$treatment.ch1[id]
table(response_class,gene_expr_group)
table(gene_expr_group, treatment_group)

plotting_df <- data.frame(gene_expr_plotting,response_class)

treatment_group_selected <- "neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF) + Trastuzumab for 1 year"
#The other treatment group only had chemo but no HER2 targeted therapy
plotting_df <- plotting_df[which(treatment_group==treatment_group_selected),]

library(hrbrthemes)
plotting_df %>%
  ggplot( aes(x=response_class, y=gene_expr_plotting, fill=response_class)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width=0.75),aes(group=response_class))+
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(curr_gene) +
  xlab(treatment_group_selected)+
  stat_compare_means(method="wilcox.test")


#*************
#*Since the clinical annotation includes the fibroblast numbers, do a regression analysis between fibroblast and PLK1/SERPINE1 levels
#*And correlate gene expression to stroma score
#*
curr_gene <- "SERPINE1"
curr_gene %in% row.names(data_RNA_HER2)
curr_gene_expr <- data_RNA_HER2[which(row.names(data_RNA_HER2)==curr_gene),]%>%as.numeric
data_clinical_HER2$stroma_cells.....ch1[is.na(data_clinical_HER2$stroma_cells.....ch1)] <- 0
plot(curr_gene_expr,data_clinical_HER2$stroma_cells.....ch1)
median(data_clinical_HER2$stroma_cells.....ch1)

#Correlations using manually annotated stroma cell numbers
df_gene_stroma <- data.frame(
  gene_value = curr_gene_expr,
  stroma_value = data_clinical_HER2$stroma_cells.....ch1
)
#Correlations using ESTIMATE stroma scores updated 03.19.24
df_gene_stroma <- data.frame(
  gene_value = curr_gene_expr,
  stroma_value = res_estimate_HER2["estimate score",]%>%as.numeric
)
ggplot(df_gene_stroma, aes(x=gene_value, y=stroma_value)) +
  geom_point()+ 
  geom_smooth(method = "lm", fill = NA)+
  stat_cor()+
  xlab(curr_gene)+
  ylab("stroma score")

#grouping by stroma cell numbers. Plot gene expression in barplots, 3 groups
idx_high <- which(data_clinical_HER2$stroma_cells.....ch1>=30)
idx_med <-which(data_clinical_HER2$stroma_cells.....ch1<30 & data_clinical_HER2$stroma_cells.....ch1>10)
idx_low <- which(data_clinical_HER2$stroma_cells.....ch1<=10)

id <- c(idx_low,idx_med,idx_high)
library(forcats)
df_gene_stroma <- data.frame(
  gene_value = curr_gene_expr[id],
  stroma_group = c(replicate(length(idx_low),"low"),replicate(length(idx_med),"medium"),replicate(length(idx_high),"high"))
)%>%mutate(stroma_group = fct_relevel(stroma_group, 
                                 "low","medium","high")) #mutate forces the levels of therapy to follow a "pre" "post" order
df_gene_stroma %>%
  ggplot( aes(x=stroma_group, y=gene_value, fill=stroma_group)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle(curr_gene) +
  xlab("")+
  stat_compare_means(label = "p.signif", method = "t.test",
                              ref.group = "low")  


# #grouping by stroma cell numbers. Plot gene expression in barplots, 2 groups
# idx_high <- which(data_clinical_HER2$stroma_cells.....ch1>=20)
# idx_low <- which(data_clinical_HER2$stroma_cells.....ch1<20)
# 
# id <- c(idx_low,idx_high)
# library(forcats)
# df_gene_stroma <- data.frame(
#   gene_value = curr_gene_expr[id],
#   stroma_group = c(replicate(length(idx_low),"low"),replicate(length(idx_high),"high"))
# )%>%mutate(stroma_group = fct_relevel(stroma_group, 
#                                       "low","high")) #mutate forces the levels of therapy to follow a "pre" "post" order
# df_gene_stroma %>%
#   ggplot( aes(x=stroma_group, y=gene_value, fill=stroma_group)) +
#   geom_boxplot() +
#   theme_ipsum() +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=11)
#   ) +
#   ggtitle(curr_gene) +
#   xlab("")+
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "low")  







