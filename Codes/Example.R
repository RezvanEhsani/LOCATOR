
#' This is an example to clustering patients from keren et al. paper:
#'                https://www.cell.com/cell/pdf/S0092-8674(18)31100-0.pdf
#' Here we just use Macrophage cell type as cells of interest.

rm(list=ls())

# invoke source codes
source(".../Codes/ExtractFeatures.R")
source(".../Codes/Locator.R")

##
# libraries
library(data.table)
library(Rfast)
library(dplyr)
library(caramellar)
library(progressr)
library(parallel)
library(tidyverse)
library(pbmcapply)
library(Rtsne)
#
library(ConsensusClusterPlus)
library(pals)
library(pheatmap)
library(ggalluvial)
library(ggpubr)
library(ggrepel) 
library(RColorBrewer)
library(NLP)
#
library(ggplot2)
library(ggfortify)
library(survival)
library(survminer)


#######################################################
# Main Data must be included below columns:
#  - Image ID
#  - cell ID
#  - celltype names
#  - X position
#  - Y position
#  - Meta celltype (1 if cell is tumor and 0 otherwise)
#  - ExistingClasses (if there are any)

DataPath = ".../Data/"

KerenData <- fread(paste(DataPath,"KerenData.csv",sep=""))
KerenData <- data.frame(KerenData)
KerenData <- KerenData[(KerenData$MixingClass!="Cold"),]
#
Data = data.frame(KerenData)
#
Data = KerenData

# Set parameters and identify columns 
SampleID_col = 2
CellID_col = 1
CellType_col = 6
CellTypesOfInterest = 'Macrophages'
MarkersOfInterest_col = NULL
Zscore = TRUE
PositivityCutoff = 0.5
MetaCellType_col = 3
X_col = 4
Y_col = 5
r = 100 
minCount = 10
ExistingClass_col = 7
#
Nseed = 1234

#Number of Cell types in Data
NCellTypes <- 11
#

colnames(Data)[SampleID_col] <- "sample_id"
colnames(Data)[CellType_col] <- "CellType"
colnames(Data)[CellID_col] <- "cell_id"
colnames(Data)[MetaCellType_col] <- "MetaCellType"
colnames(Data)[X_col] <- "Pos_X"
colnames(Data)[Y_col] <- "Pos_Y"
#
SurvivalSampleID_col = 1
SurvivalTime_col = 2
SurvivalCensored_col = 3
#
KerenSurvival <- fread(paste(DataPath,"KerenSurvival.csv",sep=""))
SurvivalData <- data.frame(KerenSurvival)
#
colnames(SurvivalData)[SurvivalSampleID_col] <- "sample_id"
colnames(SurvivalData)[SurvivalTime_col] <- "SurvivalTime"
colnames(SurvivalData)[SurvivalCensored_col] <- "Censored"


##########################################################
#                                                        #
#               Feature extraction                       #
#                                                        #
##########################################################

start_time <- Sys.time()
KerenFeatures <- Extract_Features(Data  = Data)
end_time <- Sys.time()
print(end_time - start_time)


##########################################################
#                                                        #
#            Run main function "locator"                 #
#                                                        #
##########################################################

# 
start_time <- Sys.time()
OutPuts <- locator(FeatureData = KerenFeatures, 
                   MainData = Data,
                   SurvivalData = SurvivalData,
                   Cutoff = 'Mean')
end_time <- Sys.time()
print(end_time - start_time)

#
#files and plots

# csv file for clusterd cells information
print(head(OutPuts$ClusteredData))

# csv file for samples information
print(head(OutPuts$SampleData))

# Heatmap plot
print(OutPuts$Heatmap)

# Tsne plot
print(OutPuts$Tsne)

# Box plot
print(OutPuts$Boxplot)

# Bar plot
print(OutPuts$Barplot)

# Map plot
print(OutPuts$MapPlot)

# Alluvium plot
print(OutPuts$AlluviumPlot)

# Survival plot
print(OutPuts$SurvivalPlots)
