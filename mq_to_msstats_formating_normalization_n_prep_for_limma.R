##################################################################################################################
### 
### Formating, normalization and generation of Limma-ready tabular file from MaxQuant output using MSstats #######
### Miguel Cosenza v1.0.1 
### 16.12.2019
##################################################################################################################

## This script will take you MaxQuant output (evidence.txt and proteinGroup.txt files) and an annotation file, 
# will run MSstats normalization and then generate a tsv file that can be used as Input into Limma 
# 

### Set conditions here ####

## Please provide the name of the MaxQuant output files to be processed

# Note: the files should be in the same R Project file from where the script would be executed. 

#evidence_file_location <- "evidence.txt"

#proteinGroups_file_location <- "proteinGroups.txt"

#annotation_file_location <- "annotation.csv"

# See how to create the annotation file in the README
# A more "automatic" way of creating this should be available 
# In a future version of this script.

### Please indicate if: 

## Do you want MSstats to remove proteins identified only with 1 peptide from the list?

# Default is FALSE, meaning: you want to keep even those proteins identified w 1 peptide. 

#remove_w1pep1 <- menu(c("Yes", "No"), 
#                      title= "Do you want to remove proteins identified only with 1 peptide from your analysis")

#remove_w1pep <- if(remove_w1pep1 == 1){TRUE} else {
#                  if(remove_w1pep1 == 2){
#                                    FALSE
#                                    }}

## What type of normalization to execute? 
# Default and recommended: "equalizeMedians"

# Options: "quantile", FALSE (no normalization) or "equalizeMedians"

#norm_type <- menu(c("quantile", "equalizeMedians", "No normalization"), 
#                 title= "What type  of normalization you want to execute? (equalizeMedians recommended)")

#NormalizationType <- if(norm_type == 1){"quantile"} else {
#                        if(norm_type == 2){"equalizeMedians"} else {
#                              if(norm_type == 3){
#                                    FALSE
#                              }
#                        }
#}

####################################### SCRIPT EXECUTION #######################

### Install required packages if necessary

packages <- c("dplyr", "here", "readr", "tidyr")

biopackgs <- c("MSstats")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
   install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   
   BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
   
}

### Load packages ####

library(MSstats)
library(dplyr)
library(here)
library(readr)
library(tidyr)

### Read and load input files ####

evidence <- read.table(file = here::here("evidence.txt"),
                       sep = "\t",
                       header = TRUE)

annotation <- read.csv(file = here::here("annotation.csv"),
                       header = TRUE)

proteingroups <- read.table(file = here::here("proteinGroups.txt"),
                            sep = "\t",
                            header = TRUE)

### Create label IDs for naming files according to conditions ####

# Label for the 'removePeptw1peptide' condition
#if(remove_w1pep == FALSE){
#      w1peptcond <- "protsw1pep"
#} else {
#      w1peptcond <- "noprotsw1pep"
#}


# Label for the normalization type condition

#if(NormalizationType == "equalizeMedians"){
#      normalization <- "equalizedMediansNorm"
#} else {
#      if(NormalizationType == "quantile"){
#            normalization <- "quantileNorm"
#      } else {
#            if(NormalizationType == FALSE){
#                  normalization <- "NoNormalization"
#            } else {
#                  stop("Error: You did not choose the right parameter for normalization. Maybe check for typos?",
#                       call. = FALSE)
#            }
#      }
#}

### Transform Max Quant format to MSstats format ####

if(dir.exists(here::here("MSstats_Output_data")) == FALSE){
   dir.create(here::here("MSstats_Output_data"))}

#if(file.exists(x = here::here("MSstats_Output_data/msts_data_w1pep.Rda")) == FALSE){ # check if Rda file with the data is already available. 
      
      msts_data_w1pep <- MaxQtoMSstatsFormat(evidence = evidence,
                                             annotation = annotation,
                                             useUniquePeptide=TRUE,
                                             proteinGroups = proteingroups,
                                             fewMeasurements="remove",
                                             removeProtein_with1Peptide	=TRUE)
      
      if(dir.exists(here::here("MSstats_Output_data/MSstats_formated_tables")) == FALSE){
         dir.create(here::here("MSstats_Output_data/MSstats_formated_tables"))}

      ### Correction of protein IDs ####      
      # Get the first listed proteinID when many Uniprot IDs are listed into one 
      # Protein Group (many identified peptides matching several proteins in the database)
      
      msts_formated_data <- dplyr::mutate(msts_data_w1pep,
                                          ProteinName = stringr::str_remove_all(ProteinName, ";.*$")) %>%
                            dplyr::mutate(ProteinName = stringr::str_trim(ProteinName))
      
      # File name definition
      file_name1 <- paste0("msstas_formated_data_bf_normalization.csv")
      
      write.csv(x = msts_formated_data, file = here::here(paste0("MSstats_Output_data/MSstats_formated_tables/",file_name1))) 
      
      ### Normalization ####
      
      normalized_data <- dataProcess(msts_formated_data,
                                     logTrans=2,
                                     normalization="equalizeMedians",
                                     nameStandards=NULL,
                                     address="",
                                     fillIncompleteRows=TRUE,
                                     featureSubset="highQuality",
                                     remove_uninformative_feature_outlier=TRUE,
                                     n_top_feature=3,
                                     summaryMethod="TMP",
                                     equalFeatureVar=TRUE,
                                     censoredInt="NA",
                                     cutoffCensored="minFeature",
                                     MBimpute=TRUE,
                                     remove50missing=FALSE,
                                     maxQuantileforCensored=0.999,
                                     clusters=NULL)
      
#}

### Getting tabular data into wide format as an input for Limma ####

tab_proccessed_data <- normalized_data$RunlevelData

tab_selection_msts_data <- dplyr::select(tab_proccessed_data,
                                         Protein, LogIntensities, Run = originalRUN, Group = GROUP_ORIGINAL,
                                         BioReplicate = SUBJECT_ORIGINAL) %>%
                           mutate(GROUP_REPLICATE_RUN = paste0(Group,"_",BioReplicate,"_",Run)) %>%
                           dplyr::select(Protein, LogIntensities, GROUP_REPLICATE_RUN)

tab_wide_msts_data <- tidyr::pivot_wider(data = tab_selection_msts_data,
                                  names_from = GROUP_REPLICATE_RUN, values_from = LogIntensities) # make data wide

#if(dir.exists(here::here("MSstats_Output_data")) == FALSE){
#   dir.create(here::here("MSstats_Output_data"))}

write_csv(x = tab_wide_msts_data, path = here::here("msstats_tabular_data_for_limma_input.csv"))


#message("File 'msstats_tabular_data_for_limma_input.csv' was stored into MSstats_Output_data")




