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
  
  rm(list=ls()) #Clear environment
  options(stringsAsFactors = FALSE) #Do not convert strings (vector characters) to factors in data.frames and other variables
}

###Create VIP function
{
  ### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
  ### `pls' package.
  ### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $
  
  ### Copyright ? 2006,2007 Bj?rn-Helge Mevik
  ### This program is free software; you can redistribute it and/or modify
  ### it under the terms of the GNU General Public License version 2 as
  ### published by the Free Software Foundation.
  ###
  ### This program is distributed in the hope that it will be useful,
  ### but WITHOUT ANY WARRANTY; without even the implied warranty of
  ### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ### GNU General Public License for more details.
  
  ### A copy of the GPL text is available here:
  ### http://www.gnu.org/licenses/gpl-2.0.txt
  
  ### Contact info:
  ### Bj?rn-Helge Mevik
  ### bhx6@mevik.net
  ### R?dtvetvien 20
  ### N-0955 Oslo
  ### Norway
  
  ### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
  ### some variable selection methods when multicollinearity is present,
  ### Chemometrics and Intelligent Laboratory Systems 78, 103--112.
  
  ## VIP returns all VIP values for all variables and all number of components,
  ## as a ncomp x nvars matrix.
  VIP <- function(object) {
    if (object$method != "oscorespls")
      stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
    if (nrow(object$Yloadings) > 1)
      stop("Only implemented for single-response models")
    
    SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
    Wnorm2 <- colSums(object$loading.weights^2)
    SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
    sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
  }
  
  
  ## VIPjh returns the VIP of variable j with h components
  VIPjh <- function(object, j, h) {
    if (object$method != "oscorespls")
      stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
    if (nrow(object$Yloadings) > 1)
      stop("Only implemented for single-response models")
    
    b <- c(object$Yloadings)[1:h]
    T <- object$scores[,1:h, drop = FALSE]
    SS <- b^2 * colSums(T^2)
    W <- object$loading.weights[,1:h, drop = FALSE]
    Wnorm2 <- colSums(W^2)
    sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
  }
}

#Set working directory
setwd("/Users/ioanniszervantonakis/Library/CloudStorage/Dropbox/Yanniscode/RCode/CMBEpaper2024")

#Set cell line groups
FPCLs <- c("EFM192H2BGFP", "BT474H2BGFP", "HCC1569H2BGFP", "HCC202H2BGFP", "MDA-MB-361H2BGFP")
trainCLs <- c("EFM192H2BGFP", "HCC1419H2BGFP", "HCC202H2BGFP", "MDA-MB-361H2BGFP", "UACC812H2BGFP")
testCLs <- c("BT474H2BGFP", "HCC1569H2BGFP", "HCC1954H2BGFP", "SUM225H2BGFP")
  
#List of proteins that are inhibited by >20% in at least one AR22CM-protected cell line
uniquedownlap <- read_excel("C:/Users/jakec/Desktop/Working Directory 2024/JCM Data Merged Biological Replicates.xlsx", sheet = 2, col_names = FALSE) %>% unique %>% as.matrix

#Select  data for PLSR model
library(readxl)
rawdata <- read_excel("/Users/ioanniszervantonakis/Library/CloudStorage/Dropbox/Yanniscode/RCode/CMBEpaper2024/JCM Data Merged Biological Replicates.xlsx")  %>% as.data.frame() %>% filter(Timepoint == 48, DrugConcentration != 0.3, DrugConcentration != 1)

#To predict Day4Lap/Day4DMSO cell number
day4lapdmsoCN <- read_excel("C:/Users/jakec/Desktop/Working Directory 2024/Day4LapDay4DMSO.xlsx")  %>% as.data.frame() %>% filter(CultureType == "AR22CM" | CultureType == "Monoculture")
rawdata <- merge(rawdata, day4lapdmsoCN, by.x = c("CellLine", "CultureType"))

rawdata <- rbind.data.frame(rawdata[rawdata$CultureType == "Monoculture",], rawdata[rawdata$CultureType == "AR22CM",])
rawdata <- rawdata[rawdata$CultureType == "AR22CM",]
rawdata <- rawdata[rawdata$CultureType == "Monoculture",]
rawdata <- rawdata[rawdata$DrugConcentration == 0.1,]
rawdata <- rawdata[rawdata$DrugConcentration == 0,]
rawdata <- rawdata[rawdata$CellLine %in% FPCLs,]

#traindata <- rawdata[rawdata$CellLine %in% trainCLs,]
#testdata <- rawdata[rawdata$CellLine %in% testCLs,]

#To predict lap response with DMSO protein expression
#rawdata_lap <- rawdata[rawdata$DrugConcentration == 0.1,]
#rawdata$NormalizedCellGrowth_Endpoint <- rawdata_lap$NormalizedCellGrowth_Endpoint





#Set rownames of data (needed for scoreplot)
rownames(rawdata) <- paste(rawdata$CellLine, rawdata$CultureType, rawdata$DrugConcentration)
#rownames(rawdata) <- paste(rawdata$CellLine)

#Create string of all proteins (needed for PLSR function)
#paste(colnames(rawdata[,3:392]), collapse = "+")
#paste(rep(c("mean_"), length(uniquedownlap)), uniquedownlap, sep = "") %>% as.matrix %>% paste(collapse = "+")

#Model with all proteins
model <- plsr(Day4LapDay4DMSO~mean_x14_3_3_beta_R_V+mean_x14_3_3_epsilon_M_C+mean_x14_3_3_zeta_R_V+mean_x4E_BP1_R_V+mean_x4E_BP1_pS65_R_V+mean_x4E_BP1_pT37_T46_R_V+mean_x53BP1_R_V+mean_A_Raf_R_V+mean_A_Raf_pS299_R_C+mean_ACC1_R_C+mean_ACC_pS79_R_V+mean_ACSL1_R_V+mean_ACVRL1_R_C+
                mean_ADAR1_M_V+mean_Akt_R_V+mean_Akt1_R_V+mean_Akt1_pS473_R_V+mean_Akt2_R_V+mean_Akt2_pS474_R_C+mean_Akt_pS473_R_V+mean_Akt_pT308_R_V+mean_Ambra1_pS52_R_C+mean_AMPK_a2_pS345_R_V+mean_AMPKa_R_C+mean_AMPKa_pT172_R_C+mean_Annexin_I_M_V+mean_Annexin_VII_M_V+mean_AR_R_V+mean_ARID1A_R_C+
                mean_ASNS_R_V+mean_Atg3_R_V+mean_Atg4B_R_C+mean_Atg5_R_C+mean_Atg7_R_V+mean_ATM_R_V+mean_ATM_pS1981_R_V+mean_ATP5A_M_C+mean_ATRX_R_C+mean_ATR_pS428_R_C+mean_Aurora_B_R_V+mean_Axl_R_V+mean_b_Actin_R_C+mean_b_Catenin_R_V+mean_b_Catenin_pT41_S45_R_V+mean_B_Raf_R_C+mean_B_Raf_pS445_R_V+
                mean_B7_H3_R_C+mean_B7_H4_R_C+mean_Bad_pS112_R_V+mean_Bak_R_C+mean_BAP1_M_V+mean_Bax_R_V+mean_Bcl_xL_R_V+mean_BCL2A1_R_V+mean_Beclin_R_C+mean_Bid_R_C+mean_Bim_R_V+mean_BiP_GRP78_M_C+mean_BMK1_Erk5_pT218_Y220_R_V+mean_BRD4_R_V+mean_c_Abl_R_V+mean_c_Abl_pY412_R_C+mean_c_IAP2_R_C+mean_c_Jun_pS73_R_V+
                mean_c_Kit_R_V+mean_c_Met_pY1234_Y1235_R_V+mean_c_Myc_R_C+mean_C_Raf_R_C+mean_C_Raf_pS338_R_V+mean_CA9_R_C+mean_Caspase_3_cleaved_R_C+mean_Caspase_7_cleaved__R_C+mean_Caspase_8_M_Q+mean_Caveolin_1_R_V+mean_CD134_R_V+mean_CD171_M_V+mean_CD20_R_C+mean_CD26_R_V+mean_CD29_M_V+mean_CD31_M_V+mean_CD38_R_C+
                mean_CD4_R_V+mean_CD44_R_C+mean_CD45_M_V+mean_CD49b_M_V+mean_cdc25C_R_V+mean_cdc2_pY15_R_C+mean_Cdc42_R_C+mean_CDK1_pT14_R_C+mean_CDKN2A_R_C+mean_Chk1_M_C+mean_Chk1_pS296_R_V+mean_Chk1_pS345_R_C+mean_Chk2_M_V+mean_Chk2_pT68_R_C+mean_CIITA_R_C+mean_Claudin_7_R_V+mean_COG3_R_V+mean_Collagen_VI_R_V+mean_Complex_II_Subunit_M_V+mean_Connexin_43_R_C+mean_Coup_TFII_R_C+mean_Cox_IV_R_V+mean_Cox2_R_C+mean_Creb_R_C+mean_CSK_R_C+mean_CtIP_R_V+mean_Cyclin_B1_R_V+mean_Cyclin_D1_R_C+mean_Cyclin_D3_M_V+mean_Cyclophilin_F_M_V+mean_D_a_Tubulin_R_V+mean_DAPK1_pS308_M_C+mean_DAPK2_R_C+mean_DDB_1_R_V+mean_DDR1_R_V+
                mean_DJ1_R_V+mean_DM_K9_Histone_H3_M_V+mean_DNA_Ligase_IV_R_C+mean_DNA_POLG_R_V+mean_DNMT1_R_V+mean_DRP1_R_V+mean_DUSP4_R_V+mean_DUSP6_R_C+mean_E_Cadherin_R_V+mean_eEF2_R_C+mean_eEF2K_R_V+mean_EGFR_R_V+mean_EGFR_pY1173_R_V+mean_eIF4E_R_V+mean_eIF4E_pS209_R_V+mean_eIF4G_R_C+mean_Elk1_pS383_R_C+mean_EMA_M_C+mean_Enolase_2_R_V+mean_ENY2_M_C+mean_ER_a_R_V+mean_ER_a_pS118_R_V+mean_ERCC1_M_V+mean_ERCC5_R_C+mean_Erk5_R_V+mean_ERRalpha_R_V+mean_Ets_1_R_V+mean_FAK_R_C+mean_FAK_pY397_R_V+mean_FASN_R_V+mean_FGF_basic_R_C+mean_Fibronectin_R_V+mean_FN14_R_C+mean_FOXM1_R_V+mean_FOXO3_R_V+mean_G6PD_R_V+mean_Gab2_R_V+mean_GATA3_M_V+mean_GATA6_R_V+mean_GCLC_R_C+mean_GCN5L2_R_V+mean_Gli1_R_C+mean_Gli3_R_C+mean_Glutamate_D1_2_R_V+mean_Glutaminase_R_C+mean_Granzyme_B_R_V+mean_GRB7_R_V+mean_GSK_3a_b_M_V+mean_GSK_3a_b_pS21_S9_R_V+mean_GSK_3B_R_C+mean_Gys_R_V+mean_Gys_pS641_R_V+mean_H2AX_pS140_M_C+mean_HER2_M_V+mean_HER2_pY1248_R_C+mean_HER3_R_V+mean_HER3_pY1289_R_C+mean_Heregulin_R_V+mean_HES1_R_V+mean_Hexokinase_II_R_V+mean_Histone_H3_R_V+mean_HMHA1_R_V+mean_HSP27_M_C+mean_HSP27_pS82_R_V+mean_HSP60_R_V+mean_HSP70_R_C+mean_IDO_R_C+mean_IGF1R_pY1135_Y1136_R_V+mean_IGFBP2_R_V+mean_IGFBP3_M_V+
                mean_IGFRb_R_C+mean_IL_6_R_C+mean_IR_b_R_C+mean_IRF_1_R_C+mean_IRS1_R_V+mean_IRS2_R_C+mean_JAB1_M_C+mean_Jagged1_R_V+mean_Jak2_R_V+mean_LAD1_R_V+mean_Lasu1_R_V+mean_LC3A_B_R_C+mean_Lck_R_V+mean_LDHA_R_C+mean_LRP6_pS1490_R_V+mean_Mcl_1_R_V+mean_MCT4_R_V+mean_MDM2_pS166_R_V+mean_MEK1_R_V+mean_MEK1_p_S217_S221_R_V+mean_MEK2_R_V+mean_MelanA_R_C+mean_Melanoma_gp100_R_C+mean_MERIT40_R_C+mean_MERIT40_pS29_R_V+mean_Merlin_R_C+mean_MIF_R_C+mean_MIG6_M_V+mean_MITF_R_V+mean_Mitofusin_1_R_V+mean_Mitofusin_2_R_V+mean_MLH1_M_V+mean_MLKL_R_V+mean_MMP2_R_V+mean_Mnk1_R_V+mean_MR1_M_C+mean_MRAP_R_C+mean_MSH2_R_C+mean_MSH6_R_C+mean_MSI2_R_C+mean_MTCO1_M_V+mean_mTOR_R_V+mean_mTOR_pS2448_R_C+mean_MTSS1_M_C+mean_MYH11_R_C+mean_Myosin_IIa_R_C+mean_Myosin_IIa_pS1943_R_V+mean_Myt1_R_C+mean_N_Cadherin_R_V+mean_N_Ras_M_V+mean_NAPSIN_A_R_C+mean_NDRG1_pT346_R_V+mean_NF_kB_p65_pS536_R_C+mean_Notch1_R_V+mean_Notch1_cleaved_R_V+mean_Notch3_R_C+mean_NQO1_M_V+mean_NRF2_R_C+mean_Oct_4_R_C+mean_P_Cadherin_R_C+mean_p21_R_C+mean_p27_Kip1_R_V+mean_p27_pT157_R_C+mean_p27_pT198_R_V+mean_p38_a_M_V+mean_p38_MAPK_R_V+mean_p38_MAPK__pT180_Y182_R_V+mean_p44_42_MAPK_R_V+mean_p53_R_C+mean_p70_S6K1_R_V+mean_p70_S6K_pT389_R_V+mean_p90RSK_pT573_R_C+mean_PAI_1_M_V+mean_PAICS_R_C+mean_PAK1_R_V+mean_PAK4_R_V+mean_PAR_R_C+mean_PARG_R_C+mean_PARP_R_V+mean_Patched_R_C+mean_Paxillin_R_C+mean_PCNA_M_C+mean_PD_1_R_V+mean_Pdcd4_R_C+mean_PDH_M_V+mean_PDHA1_R_V+mean_PDHK1_R_C+mean_PDK1_R_V+mean_PDK1_pS241_R_V+mean_PEA_15_R_V+mean_PEA_15_pS116_R_V+mean_PHGDH_R_C+mean_PI3K_p110_a_R_C+mean_PI3K_p110_b_M_C+mean_PKA_a_R_V+mean_PKC_a_b_II_pT638_T641_R_V+mean_PKC_b_II_pS660_R_V+mean_PKC_delta_pS664_R_V+mean_PKCa_R_V+mean_PKM2_R_C+mean_PLK1_R_C+mean_PMS2_R_V+mean_Porin_M_V+mean_PR_R_V+mean_PRAS40_M_C+mean_PRAS40_pT246_R_V+mean_PREX1_R_V+mean_PTEN_R_V+mean_PTPN12_R_V+mean_Puma_R_C+mean_PYGB_R_V+mean_PYGM_M_C+mean_Pyk2_pY402_R_C+mean_Rab11_R_C+mean_Rab25_R_V+mean_Rad23A_R_C+mean_Rad50_R_V+mean_Rad51_R_C+mean_Raptor_R_V+mean_Rb_M_Q+mean_RBM15_R_V+mean_Rb_pS807_S811_R_V+mean_Rheb_M_C+mean_Rictor_R_C+
                mean_Rictor_pT1135_R_V+mean_RIP_R_C+mean_RIP3_R_C+mean_RPA32_R_V+mean_RPA32_pS4_S8_R_C+mean_RRM1_R_C+mean_RRM2_R_C+mean_RSK_R_C+mean_S100A4_R_V+mean_S6_M_V+mean_S6_pS235_S236_R_V+mean_S6_pS240_S244_R_V+mean_SCD_M_V+mean_SDHA_R_V+mean_SF2_M_V+mean_Shc_pY317_R_V+mean_SHP_2_pY542_R_C+mean_SHP2_R_V+mean_SLC1A5_R_C+mean_Slfn11_G_C+mean_Smac_M_Q+mean_Smad1_R_V+mean_Smad3_R_V+mean_Snail_M_Q+mean_SOD1_M_V+mean_SOD2_R_V+mean_Sox2_R_V+mean_Src_M_V+mean_Src_pY416_R_V+mean_Src_pY527_R_V+mean_Stat3_R_C+mean_Stat5a_R_V+mean_Stathmin_1_R_V+mean_STING_R_V+mean_Syk_M_V+mean_Tau_M_C+mean_TAZ_R_V+mean_TFAM_R_V+mean_TFRC_R_V+mean_TIGAR_R_V+mean_Transglutaminase_M_V+mean_TRAP1_M_V+mean_TRIM25_R_C+mean_TSC1_R_C+mean_TTF1_R_V+mean_Tuberin_R_V+mean_Tuberin_pT1462_R_V+mean_TUFM_R_V+mean_Twist_M_C+mean_Tyro3_R_V+mean_U_Histone_H2B_R_C+mean_UBAC1_R_V+mean_UBQLN4_M_C+mean_UGT1A_M_V+mean_ULK1_pS757_R_C+mean_UQCRC2_M_C+mean_UVRAG_R_C+mean_VASP_R_V+mean_VAV1_R_C+mean_VEGFR_2_R_V+mean_VHL_EPPK1_M_C+mean_Vinculin_M_V+mean_Wee1_R_C+mean_Wee1_pS642_R_C+mean_WIPI1_R_C+mean_WIPI2_R_C+mean_XBP_1_G_C+mean_XIAP_R_C+mean_XPA_M_V+mean_XPF_R_C+mean_XRCC1_R_C+mean_YAP_R_C+mean_YAP_pS127_R_V+
                mean_YB1_pS102_R_V+mean_ZAP_70_R_C, data=rawdata, scale=TRUE, validation="LOO", method = "oscorespls")

#Model with only proteins inhibited by lapatinib
model <- plsr(Day4LapDay4DMSO~mean_HER2_pY1248_R_C+mean_SHP_2_pY542_R_C+mean_S6_pS235_S236_R_V+mean_S6_pS240_S244_R_V+mean_Akt_pS473_R_V+mean_Src_pY416_R_V+mean_Cyclin_B1_R_V+mean_DUSP4_R_V+mean_PLK1_R_C+mean_CDK1_pT14_R_C+mean_Rb_pS807_S811_R_V+mean_Akt_pT308_R_V+mean_p70_S6K_pT389_R_V+mean_ACC1_R_C+mean_ARID1A_R_C+mean_TFRC_R_V+mean_RRM2_R_C+mean_EGFR_pY1173_R_V+mean_FOXM1_R_V+mean_Akt1_pS473_R_V+mean_Myt1_R_C+mean_FASN_R_V+mean_NDRG1_pT346_R_V+mean_cdc25C_R_V+mean_ASNS_R_V+mean_PRAS40_pT246_R_V+mean_PAK1_R_V+mean_Pyk2_pY402_R_C+mean_mTOR_pS2448_R_C+mean_ULK1_pS757_R_C+mean_Akt2_pS474_R_C+mean_FAK_pY397_R_V+mean_Hexokinase_II_R_V+mean_PKA_a_R_V+mean_DNMT1_R_V+mean_Chk1_M_C+mean_Atg5_R_C+mean_p90RSK_pT573_R_C+mean_x4E_BP1_pT37_T46_R_V+mean_MCT4_R_V+mean_x53BP1_R_V+mean_Sox2_R_V+mean_GSK_3a_b_pS21_S9_R_V+mean_Atg4B_R_C+mean_eEF2_R_C+mean_C_Raf_pS338_R_V+mean_Stat5a_R_V+mean_TRIM25_R_C+mean_Gys_R_V+mean_PKC_b_II_pS660_R_V+mean_AR_R_V+mean_RBM15_R_V+mean_EMA_M_C+mean_GSK_3B_R_C+mean_PTPN12_R_V+mean_A_Raf_R_V+mean_eIF4G_R_C+mean_Vinculin_M_V+mean_CtIP_R_V+mean_Wee1_R_C+mean_cdc2_pY15_R_C+mean_Bad_pS112_R_V+mean_Erk5_R_V+mean_UBAC1_R_V+mean_p27_pT157_R_C+mean_HER3_pY1289_R_C+mean_NQO1_M_V+mean_Gab2_R_V+mean_S6_M_V+mean_PKCa_R_V+mean_MDM2_pS166_R_V+mean_Mitofusin_1_R_V+mean_PKM2_R_C+mean_MLH1_M_V+mean_PARG_R_C+mean_Tau_M_C+mean_CSK_R_C+mean_Transglutaminase_M_V+mean_Rab25_R_V+mean_Chk1_pS296_R_V+mean_B_Raf_pS445_R_V+mean_ACC_pS79_R_V+mean_Cyclin_D3_M_V+mean_Notch1_R_V+mean_SCD_M_V+mean_x4E_BP1_pS65_R_V+mean_MSH6_R_C+mean_ATR_pS428_R_C+mean_VHL_EPPK1_M_C+mean_DRP1_R_V+mean_FN14_R_C+mean_PMS2_R_V+mean_XPA_M_V+mean_PDK1_R_V+mean_CD171_M_V+mean_ZAP_70_R_C+mean_Syk_M_V+mean_Shc_pY317_R_V, data=rawdata, scale=TRUE, validation="CV", method = "oscorespls")

#Model for train data
model <- plsr(NormalizedCellGrowth_Endpoint~mean_HER2_pY1248_R_C+mean_SHP_2_pY542_R_C+mean_S6_pS235_S236_R_V+mean_S6_pS240_S244_R_V+mean_Akt_pS473_R_V+mean_Src_pY416_R_V+mean_Cyclin_B1_R_V+mean_DUSP4_R_V+mean_PLK1_R_C+mean_CDK1_pT14_R_C+mean_Rb_pS807_S811_R_V+mean_Akt_pT308_R_V+mean_p70_S6K_pT389_R_V+mean_ACC1_R_C+mean_ARID1A_R_C+mean_TFRC_R_V+mean_RRM2_R_C+mean_EGFR_pY1173_R_V+mean_FOXM1_R_V+mean_Akt1_pS473_R_V+mean_Myt1_R_C+mean_FASN_R_V+mean_NDRG1_pT346_R_V+mean_cdc25C_R_V+mean_ASNS_R_V+mean_PRAS40_pT246_R_V+mean_PAK1_R_V+mean_Pyk2_pY402_R_C+mean_mTOR_pS2448_R_C+mean_ULK1_pS757_R_C+mean_Akt2_pS474_R_C+mean_FAK_pY397_R_V+mean_Hexokinase_II_R_V+mean_PKA_a_R_V+mean_DNMT1_R_V+mean_Chk1_M_C+mean_Atg5_R_C+mean_p90RSK_pT573_R_C+mean_x4E_BP1_pT37_T46_R_V+mean_MCT4_R_V+mean_x53BP1_R_V+mean_Sox2_R_V+mean_GSK_3a_b_pS21_S9_R_V+mean_Atg4B_R_C+mean_eEF2_R_C+mean_C_Raf_pS338_R_V+mean_Stat5a_R_V+mean_TRIM25_R_C+mean_Gys_R_V+mean_PKC_b_II_pS660_R_V+mean_AR_R_V+mean_RBM15_R_V+mean_EMA_M_C+mean_GSK_3B_R_C+mean_PTPN12_R_V+mean_A_Raf_R_V+mean_eIF4G_R_C+mean_Vinculin_M_V+mean_CtIP_R_V+mean_Wee1_R_C+mean_cdc2_pY15_R_C+mean_Bad_pS112_R_V+mean_Erk5_R_V+mean_UBAC1_R_V+mean_p27_pT157_R_C+mean_HER3_pY1289_R_C+mean_NQO1_M_V+mean_Gab2_R_V+mean_S6_M_V+mean_PKCa_R_V+mean_MDM2_pS166_R_V+mean_Mitofusin_1_R_V+mean_PKM2_R_C+mean_MLH1_M_V+mean_PARG_R_C+mean_Tau_M_C+mean_CSK_R_C+mean_Transglutaminase_M_V+mean_Rab25_R_V+mean_Chk1_pS296_R_V+mean_B_Raf_pS445_R_V+mean_ACC_pS79_R_V+mean_Cyclin_D3_M_V+mean_Notch1_R_V+mean_SCD_M_V+mean_x4E_BP1_pS65_R_V+mean_MSH6_R_C+mean_ATR_pS428_R_C+mean_VHL_EPPK1_M_C+mean_DRP1_R_V+mean_FN14_R_C+mean_PMS2_R_V+mean_XPA_M_V+mean_PDK1_R_V+mean_CD171_M_V+mean_ZAP_70_R_C+mean_Syk_M_V+mean_Shc_pY317_R_V, data=traindata, scale=TRUE, validation="LOO", method = "oscorespls")

#group C protein
model <- plsr(NormalizedCellGrowth_Endpoint~mean_AMPK_a2_pS345_R_V+mean_Atg7_R_V+mean_ATP5A_M_C+mean_Aurora_B_R_V+mean_BAP1_M_V+mean_c_Kit_R_V+mean_Caspase_7_cleaved__R_C+mean_cdc2_pY15_R_C+mean_CDK1_pT14_R_C+mean_Chk2_M_V+mean_COG3_R_V+mean_Cox2_R_C+mean_Cyclin_B1_R_V+mean_eIF4G_R_C+mean_Elk1_pS383_R_C+mean_Gli1_R_C+mean_Heregulin_R_V+mean_IDO_R_C+mean_IGFBP3_M_V+mean_Jagged1_R_V+mean_Lck_R_V+mean_MMP2_R_V+mean_Notch1_R_V+mean_NRF2_R_C+mean_PKC_delta_pS664_R_V+mean_PTEN_R_V+mean_Rad51_R_C+mean_S6_pS235_S236_R_V+mean_S6_pS240_S244_R_V+mean_Stat5a_R_V+mean_TAZ_R_V+mean_TFAM_R_V+mean_TSC1_R_C+mean_TTF1_R_V+mean_UVRAG_R_C+mean_VAV1_R_C, data=rawdata, scale=TRUE, validation="LOO", method = "oscorespls")




#Model summary
summary(model)

###Validation plots 
validationplot(model, val.type = "R2", legendpos = "topright", estimate = "all")
validationplot(model, val.type = "RMSEP", legendpos = "topright")
#validationplot(model, val.type = "MSEP", legendpos = "topright")
#R2(model, estimate = "all")

#VIP results
vip <- t(VIP(model))
vip2 <- model$projection

###Plot VIP scores
data_plot <- cbind.data.frame(rownames(vip), vip[,1], sign(vip2[,1]))
colnames(data_plot) <- c("protein", "VIP", "sign")
data_plot <- filter(data_plot, VIP > 1.3) #Select proteins to plot

ggplot(data_plot) + 
  geom_col(mapping = aes(x = reorder(protein, -VIP), y = VIP, fill = as.character(sign)), position = "dodge") +
  ggtitle("") + ylab("VIP Score") +
  labs(x = "", y = "", fill = "Sign") +
  theme(plot.title = element_text(hjust = 0.5, size = 13), axis.text = element_text(size=15, hjust = 1), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 20), legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  ylab("VIP Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#Plot scoreplot
scoreplot(model, labels = "names")


(cor(testdata$NormalizedCellGrowth_Endpoint, predict(model, testdata, ncomp = 1), use = "complete.obs"))^2

data_plot <- cbind.data.frame(`CellLine` = rawdata$CellLine, `Measured` = rawdata$NormalizedCellGrowth_Endpoint, `Predicted` = as.vector(predict(model, rawdata, ncomp = 3))) 
data_plot$Group[rawdata$CellLine %in% testCLs] <- traindata <- "Test Set"
data_plot$Group[rawdata$CellLine %in% trainCLs] <- traindata <- "Train Set"
data_plot <- filter(data_plot, CellLine != "AU565H2BGFP")

ggplot(data = data_plot) +
  geom_point(mapping = aes(y = Predicted, x = Measured, color = CellLine, shape = Group), size = 3) +
  # xlim(c(0,7)) +
  # ylim(c(0,7)) +
  labs(y = "Measured Viability", x = "Predicted Viability", color = "Cell Line") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text = element_text(size=20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  ggtitle("")

