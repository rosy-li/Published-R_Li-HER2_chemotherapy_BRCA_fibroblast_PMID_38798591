#* ***********************
#* data from https://www.cbioportal.org/study/summary?id=breast_ink4_msk_2021
#* Analysis of the MSK dataset (metastatic breast cancer)
#* Analysis of the association between patient survival and PLK/PAI1 expression
#* Rosy Li 3/5/2024
#* ***********************
#* 

# #unzip the downloaded dataset
# untar("breast_ink4_msk_2021.tar.gz")
# data_RNA <- read.table("breast_ink4_msk_2021/data_gene_panel_matrix.txt",header = TRUE, row.names = 1)
# "SERPINE1" %in% rownames(data_RNA)

library(cBioPortalData)
# library(AnVIL)

# cbio <- cBioPortal()
# studies <- getStudies(cbio, buildReport = TRUE)
# head(studies)

cbio <- cBioPortal()
cbio













