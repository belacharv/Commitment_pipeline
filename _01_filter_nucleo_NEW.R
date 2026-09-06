################################################################################
### 01 FILTER NUCLEO - COMMITMENT PIPELINE ###
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

### The goal is to have table with non-unique sequences from where we'll create count
# table in the next step

## MERGE THE NUCLEO SEQ WITH SIMILAR AA SEQ
# In our project, we focus on the determining the commitment of a sequence to 
# one of the T cell lineages. For that we focused on non-unique sequences (seq),
# so we can observe in which samples the sequence is present. In other words,
# we study the commitment of each sequence to the cell lineage based on its 
# presence in CD4 and CD8 samples.
# Hence it is in our favor to have an overlap of sequence across our samples as big as possible.
# For that purpose, we focused on amino acid seq instead of nucleo acid seq. 
# For two nucleo sequences to be merged, they need to have similar:
  # - amino acid sequence
  # - group (CD4 or CD8)
  # - mouse (animal tag)
  # - Sample.bio.name (name tag combining group and mouse)

## Loading libraries
library(tidyverse)
library(this.path)
library(data.table)
library(NLP)
library(hash)
library(this.path)


## LOAD DATA
# input: 
  # name = name of the file
  # directory = adress to the directory
load_data <- function(name, directory){
  setwd(directory)
  dt <- fread(name) # inicial input table
  # add Umi proportion, which describes UMI proportion of this specific seq from
  # the whole sample
  dt[, Umi.proportion  := Umi.count  / sum(Umi.count),  by = Sample.bio.name]
  return(dt)
}

# Recount proportions per sample after filtering
recountProportions <- function(dt) {
  dt[, Umi.proportion  := Umi.count  / sum(Umi.count),  by = Sample.bio.name]
  dt
}

# Full pipeline: merge nucleotide duplicates -> filter singletons -> recount
# input:
  # input_data = datatable from the current dataset after customization 
    # (for example after renaming columns so it fits the pipeline)
# output = table with merged nucleo seq
generate_result <- function(input_data, dt) {
  dt <- data.table(input_data)
  
  # Merge rows with same aa seq + same sample (summing counts)
  merged <- dt[, .(
    Umi.count      = sum(Umi.count),
    Umi.proportion = sum(Umi.proportion)
  ), by = .(CDR3.amino.acid.sequence, Sample.bio.name, group, Mouse)]
  
  # Compute group size (rows per aa sequence) as its own column
  merged[, n_rows := .N, by = CDR3.amino.acid.sequence]
  
  # Filter: keep only sequences with >1 UMI or appearing in >1 row
  merged <- merged[Umi.count > 1 | n_rows > 1]
  merged[, n_rows := NULL]
  
  # Recount proportions cleanly per sample
  recountProportions(merged)
}