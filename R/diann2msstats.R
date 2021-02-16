### DIA-NN to MSstats formatting for LFQ normalization ----

## Install required packages if necessary ----

packages <- c("dplyr", "here", "readr", "tidyr", "diann", "tidyverse")

biopackgs <- c("MSstats")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
      install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(biopackgs, rownames(installed.packages()))) > 0){
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      
      BiocManager::install(setdiff(biopackgs, rownames(installed.packages())))
      
}


# Modified on 16.02.2021 ----
#install.packages("devtools")
library(tidyverse)
library(dplyr)

##Genes.MaxLFQ gives already normalized MQ-LFQ data for every protein group and is recommended
## Genes.MaxLFQ would be used to extract the quantitative information per feature from the DIANN output

DIANN_to_MSstats <- function(diann_data, annotation_file){

#remove filepath from File.Name
diann_data1 <- mutate(diann_data, File.Name = str_replace(diann_data[[1]], ".*\\\\", ""))
diann_data1 <- mutate(diann_data1, File.Name = str_replace(diann_data1[[1]], ".raw.mzml$", ""))

#select relevant columns -> choose which parameter to use for quantitation

for_msstats_prot <- diann_data1 %>% 
      group_by(Stripped.Sequence, 
               Protein.Group, 
               Precursor.Charge, 
               File.Name, 
               Genes.MaxLFQ) %>% 
      summarize()

for_msstats_prot1 <- mutate(for_msstats_prot,
                            File.Name = str_remove(File.Name, ".raw$"))

#merge with MSstats annotation file
annotation <- annotation_file
colnames(annotation)[colnames(annotation) == "Filename"] <- "File.Name"
for_msstats_prot2 <- left_join(for_msstats_prot1, annotation, by = "File.Name")

#change column names to MSstats format
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "BioReplicate"] <- "BioReplicate"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Stripped.Sequence"] <- "PeptideSequence"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Protein.Group"] <- "ProteinName"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Precursor.Charge"] <- "PrecursorCharge"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "File.Name"] <- "Run"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Precursor.Quantity"] <- "Intensity"
colnames(for_msstats_prot2)[colnames(for_msstats_prot2) == "Genes.MaxLFQ"] <- "Intensity"

#add missing columns
for_msstats_prot2$IsotopeLabelType <- "L"
for_msstats_prot2$ProductCharge <- NA
for_msstats_prot2$FragmentIon <- NA

#reorder to MSstats format
for_msstats_prot3 <- for_msstats_prot2[, c("ProteinName", "PeptideSequence", "PrecursorCharge", 
                                           "FragmentIon", "ProductCharge", "IsotopeLabelType", 
                                           "Condition", "BioReplicate", "Run", "Intensity")]

return(for_msstats_prot3)
}
