# MSstats label-free preprocessing

This repo contains a script and a Rmd file for the pre-processing and normalization of MaxQuant output files through the `MSstats` R package. 

## Instructions for using the script and to create a reproducible report of your normalization step  
1. Download/clone the contents of this repo into your local computer. This should create a R project folder with the script to run the preprocessing. 

2. Add the next three MaxQuant output files into this folder 

1. The R script would take three files, two of which can be taken from the txt output folder from your MaxQuant analysis:

- evidence.txt
- proteinGroups.txt
- annotation.csv (not included in the MaxQuant txt folder).

**NOTE**: These files should be in the same folder as the R script, and this folder should be an initiated RStudio project (There should be a .Rproj file in the same folder). 

2. 



### Creating the annotation file (provisional)  

At some point I will add something in the processing script (or create a new script) for trying to automotize this, but for now, the annotation.csv file can be created manually as described now.







