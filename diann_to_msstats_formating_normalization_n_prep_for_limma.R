##################################################################################################################
### 
### Formating, normalization and generation of Limma-ready tabular file from DIANN output using MSstats #######
### Miguel Cosenza v1.0
### 16.02.2021
##################################################################################################################

## This script will take your DIANN output and an annotation file, 
# will run MSstats normalization and then generate a tsv file that can be used as Input into Limma 
# 

### Set conditions here ####

### Set the parameters for dataProcess (summarization and normalization) function ----

logTrans = 2
normalization = "equalizeMedians"
nameStandards = NULL
address = ""
fillIncompleteRows = TRUE
featureSubset = "all"
remove_uninformative_feature_outlier = FALSE
n_top_feature = 3
summaryMethod = "TMP"
equalFeatureVar = TRUE
censoredInt = "NA"
cutoffCensored = "minFeature"
MBimpute = TRUE
remove50missing = FALSE
maxQuantileforCensored = 0.999
clusters = NULL

## Please provide the name of the DIANN output files to be processed

diann_file <- "DIANN_results.tsv" # this file should be located on the same RStudio project 
                                 # folder as this script.

# Note: the files should be in the same R Project file from where the script would be executed. 

# See how to create the annotation file in the README
# A more "automatic" way of creating this should be available 
# In a future version of this script.


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

diann_data <- suppressWarnings(
               read_tsv(file = here::here(diann_file)))

annotation <- read.csv(file = here::here("annotation_diann.csv"),
                         header = TRUE,
                         stringsAsFactors = FALSE)

annotation <- mutate(annotation,
                     Filename = str_remove(Filename, "ECMO"))

write.csv(annotation, file = here::here("samples/annotation_diann.csv"),
          row.names = FALSE, quote = FALSE)


### Transform DIANN format to MSstats format ####

if(dir.exists(here::here("MSstats_Output_data")) == FALSE){
   dir.create(here::here("MSstats_Output_data"))}

source("R/diann2msstats.R")

diann_msstats_formated <- DIANN_to_MSstats(diann_data, annotation)

      if(dir.exists(here::here("MSstats_Output_data/MSstats_formated_tables")) == FALSE){
         dir.create(here::here("MSstats_Output_data/MSstats_formated_tables"))}

      ### Correction of protein IDs ####      
      # Get the first listed proteinID when many Uniprot IDs are listed into one 
      # Protein Group (many identified peptides matching several proteins in the database)
      
      msts_formated_data <- dplyr::mutate(diann_msstats_formated,
                                          ProteinName = stringr::str_remove_all(ProteinName, ";.*$")) %>%
                            dplyr::mutate(ProteinName = stringr::str_trim(ProteinName))
      
      # File name definition
      file_name1 <- paste0("msstas_formated_diann_data_bf_normalization.csv")
      
      write.csv(x = msts_formated_data, file = here::here(paste0("MSstats_Output_data/MSstats_formated_tables/",file_name1))) 
      
      ### Normalization #####
      
      normalized_data <- dataProcess(msts_formated_data,
                                     logTrans=logTrans,
                                     normalization=normalization,
                                     nameStandards=nameStandards,
                                     address=address,
                                     fillIncompleteRows=fillIncompleteRows,
                                     featureSubset=featureSubset,
                                     remove_uninformative_feature_outlier=remove_uninformative_feature_outlier,
                                     n_top_feature=n_top_feature,
                                     summaryMethod=summaryMethod,
                                     equalFeatureVar=equalFeatureVar,
                                     censoredInt=censoredInt,
                                     cutoffCensored=cutoffCensored,
                                     MBimpute=MBimpute,
                                     remove50missing=remove50missing,
                                     maxQuantileforCensored=maxQuantileforCensored,
                                     clusters=maxQuantileforCensored)
      
      # File name definition
      file_name2 <- paste0("msstas_formated_data_after_normalization.csv")
      
      write.csv(x = msts_formated_data, file = here::here(paste0("MSstats_Output_data/MSstats_formated_tables/",file_name2))) 
      
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

write_csv(x = tab_wide_msts_data, path = here::here("msstats_tabular_data_for_limma_input.csv"))






