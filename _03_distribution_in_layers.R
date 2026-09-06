################################################################################
### 03 COUNT TABLE LAYER DISTRIBUTION ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate


## Given the count table with the sequences and their abundance in CD4 and CD8 
# samples, we divided the data based on in how many samples the sequence is 
#  present (layer) and how many of those were CD8 
# (categories from CD8_0 = 0 CD8, the sequence is only in CD4 samples)

# loading libraries
library(tidyverse)

# layers = clonal abundance of each sequence
  # second layer contains sequences that are present in exactly two samples
  # third layer contains seq present in three samples etc.
  # n-th layer contains seq present in all samples where n is number of samples
  # also used term "depth" in more computational aspect of graphs
  # but both terms are interchangable

# groups CD8_0 to CD8_n
  # these group represent the exact distribution within a layer considering
  # CD4 or CD8 samples
  # for example, in second layer, we have three scenarios happening
    # CD4 + CD4
    # CD4 + CD8
    # CD8 + CD8 
  # each of these scenarios can be distribed by the number of CD8 samples
    # the CD8 variable was chosen without loss of generality (WLOG)
  # example
    # CD4 + CD4 = CD8_0
    # CD4 + CD8 = CD8_1
    # CD8 + CD8 = CD8_2

# function returning the number of samples in each category in particular layer
# input:
  # Table = count table generated from 02 step
  # layer = clonal abundance; in how many samples is the seq present
  # n_cd4 = number of CD4 samples in the dataset
  # n_cd8 = number of CD8 samples in the dataset
# output: combinations of groups template
generate_combinations <- function(table,layer, n_cd4, n_cd8) {
  combinations <- c()
  # filter out only the sequences that are present in exactly the number of samples as the current layer
  datatable <- dplyr::filter(table, table$Occurrence == layer)
  # initial set-up
  combinations <- c("depth" = layer,"SeqCount" = nrow(datatable))
  for (col in 0:n_cd8){
    # each column correspond to the exact number of CD8 samples 
    # that the sequence is present in 
    col_name <- paste0("CD8_", col)
    combinations[[col_name]] <- NA
  }
  for (i in 0:layer) {
    # cd4_count = how many CD4 samples are present in the group
    # cd8_count = how many CD8 samples are present in the group
    if (layer >= i ){ # 
      cd4_count <- layer - i # number of CD4 samples in that group in that layer
      cd8_count <- i # # number of CD8 samples in that group in that layer
    }
    else {
      cd8_count <- n_cd8
      cd4_count <- layer - n_cd8
    }
    # Check if combination is possible given sample sizes
    if (cd4_count <= n_cd4 && cd8_count <= n_cd8) {
      # starting from 3rd column as first two are depth and SeqCount
      # fills in number of sequences in group with cd4_count CD4 samples and 
      # cd8_count CD8 samples
      combinations[cd8_count+3] <- length(which(
        datatable$Occurr_CD4 == cd4_count & datatable$Occurr_CD8 == cd8_count))  # Initialize count to 0
    }
  }
  return(combinations)
}


# prepare a table that will be a template for the plot
# stores abundance of sequences in each layer
# input:
  # countTable = count table generated from 02 step
  # n_cd4 = number of CD4 samples in the dataset
  # n_cd8 = number of CD8 samples in the dataset
  # finalTable = at input empty table
# output: table with each group CD8_0 to CD8_n and number of seq in each
  # layer for that group
table_prep <- function(countTable, n_cd4,n_cd8,finalTable){
  max_layer = n_cd4+n_cd8
  for (layer in 2:max_layer){
    print(layer)
    new_row = generate_combinations(countTable,layer,n_cd4,n_cd8)
    finalTable <- rbind(finalTable,new_row)
  }
  finalTable = data.frame(finalTable)
  return(finalTable)
}
