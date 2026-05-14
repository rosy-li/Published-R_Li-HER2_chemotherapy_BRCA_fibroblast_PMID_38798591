###Install packages and load libraries
{ 
  #install.packages("picante")
  #install.packages("sjPlot")
  #install.packages("xlsx")
  #install.packages("pls")
  
  #library(sjPlot)
  library(picante) # for doing matrix correlations
  library(readr) # for reading csv files
  library(matrixStats)
  library(factoextra) #visualize results from multidimensional analyses
  library(plyr)
  library(dplyr)
  library(gridExtra)
  library(ggplot2)
  library(gplots) # for heatmaps
  library(xlsx)
  library(pls)
  
  #install.packages("devtools")
  library(devtools)
  #install_github("vqv/ggbiplot")
  library(ggbiplot)
  
  # setwd("/Users/ioanniszervantonakis/Library/CloudStorage/Dropbox/Yanniscode/RCode/CMBEpaper2024")
  rm(list=ls()) #Clear environment
  options(stringsAsFactors = FALSE) #Do not convert strings (vector characters) to factors in data.frames and other variables
}

# setwd("/Users/ioanniszervantonakis/Library/CloudStorage/Dropbox/Yanniscode/RCode/CMBEpaper2024")

library(readxl)

#Data is protein expression FC from 0.1LapMonoculture to 0.1LapAR22CM conditions
FC_CM_M_0_1 <- read_excel("heatmapdata.xlsx")  %>% as.data.frame()
# dataforRosy <- read_excel("heatmapdata_for_Rosy.xlsx")%>%as.data.frame()

#Select list of proteins with >1.2 or <0.8 FC from 0.1LapMonoculture to 0.1LapAR22CM conditions in at least one FP cell line 
uniqueDE <- read_excel("heatmapdata.xlsx", sheet = 2, col_names = FALSE) %>% unique %>% as.matrix

DEGs_to_plot <- read_excel("heatmapdata_for_Rosy.xlsx", sheet = 2, col_names = FALSE) %>% unique %>% as.matrix
DEGs_to_plot <- gsub("'","",DEGs_to_plot)

#Select data for FP cell lines
FC_unique <- FC_CM_M_0_1[FC_CM_M_0_1$Protein %in% uniqueDE,c(1,3,4,6,8,9)]

#Log transform data for heatmap
hm_FC_CM_M_0_1 <- log2(FC_unique[,2:6])

#Shorten protein names (remove _R_V from end of protein names)
rownames(hm_FC_CM_M_0_1) <- substr(FC_unique[,1],1,nchar(FC_unique[,1])-4)


#Plot heatmap
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

ylim1=-1 # set limits
ylim1i=-.5
ylim2i=.5
ylim2=1 # set limits
dylim=0.01

col_breaks = c(seq(ylim1,ylim1i,length=34),  # for blue
               seq(ylim1i+dylim,ylim2i,length=33),           # for white
               seq(ylim2i+dylim,ylim2,length=34))             # for red

# pdf(file="heatmap_edited_by_RL_0329.pdf")
dev.new()
# rearrange the data matrix for heatmaps
hm_plotting <- t(as.matrix(hm_FC_CM_M_0_1[gsub('.{4}$', '',DEGs_to_plot),]))

hr_sp=hclust(dist(hm_plotting,method="euclidean"),method="average")
hc_sp=hclust(dist(t(hm_plotting),method="euclidean"),method="average") # transpose to compute distance


heatmap.2(hm_plotting,  
          scale="none", trace="none", 
          margins= c(10,7), 
          col=my_palette, 
          breaks=col_breaks, 
          cexRow=1, 
          # dendrogram = "none",
          Rowv=as.dendrogram(hr_sp),
          Colv=as.dendrogram(hc_sp),
          # Colv = FALSE,
          adjCol = c(NA,0.5),#adjust the alignment of the column names
          cexCol = 1,
          keysize =1,
          # lhei=c(2,6), lwid=c(2,10),#colorkey size
          key.par=list(cex=0.5,mar=c(3,3,3,3)))#adjust the size and margins of the color key
# this is for normalized ratios
dev.off()








