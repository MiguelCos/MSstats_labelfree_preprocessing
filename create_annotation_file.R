<<<<<<< HEAD
##################################################################################################################
### 
### 'Interactive' script to create the annotation.csv file required as input into MSstats                  #######
### Miguel Cosenza v1.0 13.12.2019
###
##################################################################################################################

# Small README note ---- 

# Please execute this script to create your annotation file for MSstats and corroborate that your sample names/codes 
# correspond with your experimental conditions 

# ----

sample_init <- readline("Please indicate the two-letter identifier for your samples (i.e. 'MC') ")

n_samples <- readline("How many samples do you have? ")

min_code_num <- as.numeric(readline("Please indicate the number code of your first sample (i.e. 21) "))

max_code_num <- as.numeric(readline("Please indicate the number code of your last sample (i.e. 38) "))

bioreplicates <- menu(c("Yes", "No"), 
                      title= "Is every sample representative of an unique Biological replicate")

n_conditions <- as.numeric(readline("Please indicate the number of conditions in your experiment \n (p.e. treated vs. control would be two conditions)"))

equal_n_samples <- menu(c("Yes", "No"), 
                      title= "Do you have the same number of samples for every experimental condition?")


if(equal_n_samples == 2){
      condition <- vector(mode = "character")
      n_per_condition <- vector(mode = "numeric")
      
      condition_col <- vector(mode = "character")
      for(i in 1:n_conditions){
            condition[i] <- readline(paste0("Type the name/code of your condition_",as.character(i),"\n "))
            n_per_condition[i] <- as.numeric(readline(paste0("How many samples are in your condition_",as.character(i),"\n ")))
      }
      condition_data <- data.frame(condition = as.character(condition),
                                   n_cond = n_per_condition)
     
      for(i in 1:n_conditions){
            condition_col <- c(condition_col,as.character(rep(condition_data[i,1],condition_data[i,2])))
      }
} else {if(equal_n_samples == 1){
   condition <- vector(mode = "character")
   n_per_condition <- as.numeric(readline("Please indicate the number of samples per condition "))
   
   condition_col <- vector(mode = "character")
   for(i in 1:n_conditions){
      condition[i] <- readline(paste0("Type the name/code of your condition_",as.character(i),"\n "))
   }
   
   condition_data <- data.frame(condition = as.character(condition),
                                n_cond = rep(n_per_condition, n_conditions))
   
      for(i in 1:n_conditions){
      condition_col <- c(condition_col,as.character(rep(condition_data[i,1],condition_data[i,2])))
   }
   
}}

annotation_tab <- data.frame(Raw.file = paste0(sample_init, 
                                               seq(from = min_code_num, 
                                                   to = max_code_num), 
                                                   sep = ""),
                             Condition = condition_col,
                             BioReplicate = seq(from = 1,
                                                to = n_samples),
                             IsotopeLabelType = rep("L", n_samples))

write.table(annotation_tab,
            file = "annotation.csv",
            row.names = FALSE,
            sep = ",")

message("annotation.csv file created. Please corroborate that your sample names correspond to your experimental conditions")
=======
##################################################################################################################
### 
### 'Interactive' script to create the annotation.csv file required as input into MSstats                  #######
### Miguel Cosenza v1.0 13.12.2019
###
##################################################################################################################

# Small README note ---- 

# Please execute this script to create your annotation file for MSstats and corroborate that your sample names/codes 
# correspond with your experimental conditions 

# ----

sample_init <- readline("Please indicate the two-letter identifier for your samples (i.e. 'MC') ")

n_samples <- readline("How many samples do you have? ")

min_code_num <- as.numeric(readline("Please indicate the number code of your first sample (i.e. 21) "))

max_code_num <- as.numeric(readline("Please indicate the number code of your last sample (i.e. 38) "))

bioreplicates <- menu(c("Yes", "No"), 
                      title= "Is every sample representative of an unique Biological replicate")

n_conditions <- as.numeric(readline("Please indicate the number of conditions in your experiment \n (p.e. treated vs. control would be two conditions)"))

equal_n_samples <- menu(c("Yes", "No"), 
                      title= "Do you have the same number of samples for every experimental condition?")


if(equal_n_samples == 2){
      condition <- vector(mode = "character")
      n_per_condition <- vector(mode = "numeric")
      
      condition_col <- vector(mode = "character")
      for(i in 1:n_conditions){
            condition[i] <- readline(paste0("Type the name/code of your condition_",as.character(i),"\n "))
            n_per_condition[i] <- as.numeric(readline(paste0("How many samples are in your condition_",as.character(i),"\n ")))
      }
      condition_data <- data.frame(condition = as.character(condition),
                                   n_cond = n_per_condition)
     
      for(i in 1:n_conditions){
            condition_col <- c(condition_col,as.character(rep(condition_data[i,1],condition_data[i,2])))
      }
} else {if(equal_n_samples == 1){
   condition <- vector(mode = "character")
   n_per_condition <- as.numeric(readline("Please indicate the number of samples per condition "))
   
   condition_col <- vector(mode = "character")
   for(i in 1:n_conditions){
      condition[i] <- readline(paste0("Type the name/code of your condition_",as.character(i),"\n "))
   }
   
   condition_data <- data.frame(condition = as.character(condition),
                                n_cond = rep(n_per_condition, n_conditions))
   
      for(i in 1:n_conditions){
      condition_col <- c(condition_col,as.character(rep(condition_data[i,1],condition_data[i,2])))
   }
   
}}

annotation_tab <- data.frame(Raw.file = paste0(sample_init, 
                                               seq(from = min_code_num, 
                                                   to = max_code_num), 
                                                   sep = ""),
                             Condition = condition_col,
                             BioReplicate = seq(from = 1,
                                                to = n_samples),
                             IsotopeLabelType = rep("L", n_samples))

write.table(annotation_tab,
            file = "annotation.csv",
            row.names = FALSE,
            sep = ",")

message("annotation.csv file created. Please corroborate that your sample names correspond to your experimental conditions")
>>>>>>> b8a34b150f5a48f649f42e3afa5f996e8b0a6cde
