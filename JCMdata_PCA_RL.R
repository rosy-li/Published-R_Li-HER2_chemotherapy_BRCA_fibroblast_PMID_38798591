#*************
#*Script to analyze the JCM data sheet 1 biological replicates
#*Conditions: all cell lines except SUM
#*Rosy Li 3.12.2024
#*************
#*

#***
#*Load RPPA mean data
#*

library(dplyr)
library(readxl)
rawdata <- readr::read_csv("JCM Data Merged Biological Replicates_sheet1_RL.csv")  %>% 
  as.data.frame() %>% 
  filter(Timepoint == 48, 
         DrugConcentration != 0.3, 
         DrugConcentration != 1,
         CultureType %in% c("AR22CM","Monoculture"),
         CellLine != "SUM225H2BGFP")



table(rawdata$Timepoint)#all at 48
table(rawdata$DrugConcentration,rawdata$CultureType)

RPPA_means_data <- rawdata[,-c(1:7)]
row.names(RPPA_means_data) <- paste0(rawdata$CellLine,"_",rawdata$DrugConcentration,"_",rawdata$CultureType)
row.names(rawdata) <- row.names(RPPA_means_data)

#*** 
#*PCA for all points
#*

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('PCAtools')
# library("PCAtools")

# p <- pca(t(RPPA_means_data), metadata=rawdata[,c(1:3)])
# 
# biplot(p)
# 
# biplot(p,
#        # loadings parameters
#        showLoadings = TRUE,
#        lengthLoadingsArrowsFactor = 1.5,
#        sizeLoadingsNames = 4,
#        colLoadingsNames = 'red4',
#        # other parameters
#        lab = NULL,
#        colby = 'CellLine', #colkey = c('ER+'='royalblue', 'ER-'='red3'),
#        # hline = 0, vline = c(-25, 0, 25),
#        # vlineType = c('dotdash', 'solid', 'dashed'),
#        gridlines.major = FALSE, gridlines.minor = FALSE,
#        pointSize = 5,
#        legendPosition = 'left', legendLabSize = 14, legendIconSize = 8.0,
#        shape = 'DrugConcentration', #shapekey = c('Grade 1'=15, 'Grade 2'=17, 'Grade 3'=8),
#        drawConnectors = FALSE,
#        title = 'PCA bi-plot',
#        # subtitle = 'PC1 versus PC2',
#        # caption = '27 PCs ≈ 80%'
#        )
drugconcSelected=0
i <- which(rawdata$DrugConcentration==drugconcSelected)


PCA_data <- RPPA_means_data[i,]
pca_res <- prcomp(PCA_data, scale. = TRUE)

library(plotly)
library(ggfortify)
autoplot(pca_res)
p <- autoplot(pca_res, data = rawdata[i,], 
              colour = "CellLine", 
              shape = "CultureType",
              )+
  ggtitle(paste("Drug Concentration =",drugconcSelected))
ggplotly(p)



#***
#*Load RPPA single read data
library(dplyr)
library(readxl)
# install.packages("readxl")
rawdata <- read_excel("Raw RPPA Data (1).xlsx")
colnames(rawdata)[1:20]
table(rawdata$Timepoint)#all at 48
table(rawdata$DrugConcentration,rawdata$CultureType)

rawdata_subset <- rawdata%>%filter(Timepoint == 48, 
                                   DrugConcentration != 0.3, 
                                   DrugConcentration != 1,
                                   CultureType %in% c("AR22CM","Monoculture"),
                                   CellLine != "SUM225H2BGFP",
                                   CellLine != "UACC893H2BGFP")

table(rawdata_subset$CellLine,rawdata_subset$DrugConcentration,rawdata_subset$CultureType)

rawdata_subset$SampleDescription[which(rawdata_subset$CellLine=="HCC1569H2BGFP" & rawdata_subset$DrugConcentration=="0.1" & 
                                  rawdata_subset$CultureType=="Monoculture")]

RPPA_single_data <- rawdata_subset[,-c(1:6)]
# row.names(RPPA_single_data) <- paste0(rawdata$CellLine,"_",rawdata$DrugConcentration,"_",rawdata$CultureType)
row.names(RPPA_single_data) <- rawdata_subset$SampleDescription
# row.names(rawdata) <- row.names(RPPA_means_data)

#***
#*Select the drug or culture type grouping here
drugconcSelected=0
i <- which(rawdata_subset$DrugConcentration==drugconcSelected)

PCA_data <- RPPA_single_data[i,]
pca_res <- prcomp(PCA_data, scale. = TRUE)

library(plotly)
library(ggfortify)
autoplot(pca_res)
p <- autoplot(pca_res, data = rawdata_subset[i,], 
              colour = "CellLine", 
              shape = "CultureType",
)+
  ggtitle(paste("Drug Concentration =",drugconcSelected))
ggplotly(p)

#plot the color based on AR22 protection
protection <- c()
protection[which(rawdata_subset$CellLine%in%c("BT474H2BGFP",
                                              "EFM192H2BGFP",
                                              "HCC1569H2BGFP",
                                              "MDA-MB-361H2BGFP"))] <- "group1"

protection[which(rawdata_subset$CellLine%in%c("AU565H2BGFP",
                                              "HCC1419H2BGFP",
                                              "HCC1954H2BGFP",
                                              "HCC202H2BGFP",
                                              "UACC812H2BGFP"))] <- "group2"
rawdata_subset_protection <- data.frame(protection,rawdata_subset)


p <- c()

for (xidx in c(1:10)){
  for (yidx in c(1:10)){
    if (xidx > yidx){
      p[[10*(xidx-1)+yidx]] <- autoplot(pca_res,
                                        x=xidx,
                                        y=yidx,
                                        data = rawdata_subset_protection[i,],
                                        colour = "protection",
                                        shape = "CultureType",
                                        )+
        ggtitle(paste("Drug Concentration =",drugconcSelected))+
        guides(col="none",shape="none")
    }
    
  }
}

ggplotly(p[[11]])


library(gridExtra)
do.call("grid.arrange", c(p, ncol=10))

#single plot
p <- autoplot(pca_res,
              x=8,
              y=7,
              data = rawdata_subset_protection[i,],
              colour = "protection",
              shape = "CultureType",
)+
  ggtitle(paste("Drug Concentration =",drugconcSelected))
ggplotly(p)














#**********************
#*Plot by culture type
culturetypeSelected="Monoculture"
j <- which(rawdata_subset$CultureType==culturetypeSelected)
rawdata_subset$DrugConcentration <- as.character(rawdata_subset$DrugConcentration)

#***
#*Select the drug or culture type grouping here
PCA_data <- RPPA_single_data[j,]
pca_res <- prcomp(PCA_data, scale. = TRUE)

library(plotly)
library(ggfortify)
autoplot(pca_res)
p <- autoplot(pca_res, data = rawdata_subset[j,], 
              colour = "CellLine", 
              shape = "DrugConcentration",
)+
  ggtitle(paste("Culture Type =",culturetypeSelected))
ggplotly(p)

#plot the color based on AR22 protection
protection <- c()
protection[which(rawdata_subset$CellLine%in%c("BT474H2BGFP",
                                              "EFM192H2BGFP",
                                              "HCC1569H2BGFP",
                                              "MDA-MB-361H2BGFP"))] <- "group1"

protection[which(rawdata_subset$CellLine%in%c("AU565H2BGFP",
                                              "HCC1419H2BGFP",
                                              "HCC1954H2BGFP",
                                              "HCC202H2BGFP",
                                              "UACC812H2BGFP"))] <- "group2"
rawdata_subset_protection <- data.frame(protection,rawdata_subset)


p <- c()

for (xidx in c(1:10)){
  for (yidx in c(1:10)){
    if (xidx > yidx){
      p[[10*(xidx-1)+yidx]] <- autoplot(pca_res,
                                        x=xidx,
                                        y=yidx,
                                        data = rawdata_subset_protection[j,],
                                        colour = "protection",
                                        shape = "DrugConcentration",
      )+
        ggtitle(paste("Culture Type = ",culturetypeSelected))+
        guides(col="none",shape="none")
    }
    
  }
}

ggplotly(p[[11]])


library(gridExtra)
do.call("grid.arrange", c(p, ncol=10))

#single plot
p <- autoplot(pca_res,
              x=6,
              y=5,
              data = rawdata_subset_protection[j,],
              colour = "protection",
              shape = "DrugConcentration",
)+
  ggtitle(paste("Culture Type =",culturetypeSelected))
ggplotly(p)
































