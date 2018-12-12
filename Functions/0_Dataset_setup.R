### 26 July 2018: Created to combine specimen list with metadata from museum databases

##############################################################
"########################  Data Set up  ######################"
##############################################################

## Load necessary libraries
library(dplyr)
library(stringr)
library(stringi)
library(data.table)

## Load in metadata
setwd("~/Documents/Ariel's Folder/Research/3 Australian Rodent Evo/Chp3_Crania_Morph/data/0_museum_data")

# Museum datasets
AM_all_metadata <- read.csv("Sydney_muridae_skulls.csv", header = T)
MV_all_metadata <- read.csv("Melbourne_muridae_skulls.csv", header = T)
QM_all_metadata <- read.csv("Brisbane_muridae_skulls.csv", header = T)

# Missing info from Melbourne (MV)
MV_extra_metadata <- read.csv("MelbourneUnknowns.csv", header = T)

## Load in specimen list
specimens <- read.csv("Data_forMorphosource.csv", header = T)


##############################################################
"########################  Data Joins  #######################"
##############################################################

## Set up museum metadata to have comparable column names / info
# Sydney / Australia Museum
AM_all_metadata$CatNum <- gsub("[.]", "", AM_all_metadata$RegNum)

# Melbourne / Museum Victoria
MV_all_metadata$CatNum <- paste(MV_all_metadata$Reg, MV_all_metadata$RegNum, sep = "")
MV_extra_metadata$CatNum <- paste(MV_extra_metadata$Reg, MV_extra_metadata$RegNum, sep = "")

## Join sex and age information to Morphosource file
AM_spec <- merge(specimens, AM_all_metadata, by = "CatNum")
AM_spec_sm <- AM_spec[ ,c(1:4, 16)]

MV_spec1 <- merge(specimens, MV_all_metadata, by = "CatNum")
MV_spec1_sm <- MV_spec1[ ,c(1:4, 9)]

MV_spec2 <- merge(specimens, MV_extra_metadata, by = "CatNum")
MV_spec2_sm <- MV_spec2[ ,c(1:4, 8)]

QM_spec <- merge(specimens, QM_all_metadata, by = "CatNum")
QM_spec_sm <- QM_spec[ ,c(1:4, 12)]

## Combine them back together

