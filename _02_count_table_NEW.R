################################################################################
### 02 COUNT TABLE - COMMITMENT PIPELINE #####
################################################################################
## Author: Bela Charvatova
## Laboratory of Adaptive Immunity, IMG CAS
## Project Cell Fate

### The goal is to provide a count table where 
# rows are non-unique sequences (present in more than one sample)
# in columns are samples
# if a sequence x is present in sample y:
#count_table[x,y] = 1 (or there is number of UMI - check out function
#umiCountSample when umi_value = TRUE)
# if not
#count_table[x,y] = 0 

# This count table is essential for our further equations and estimations as we
# are interested only in those aa seq, that are present across at least two 
# samples to apply them to our probabilistic model

# loading libraries
library(xlsx)
library(tidyverse)
library(tibble)
library(data.table)


## creating a table with metadata info that can be added to the end of the 
## count table
## last 6 rows
## input : 
  # input_data_01 = datatable from current dataset edited in 01 step
  # sample_names = list containing Sample.bio.name, unique identifier for each sample
## output: a datatable with
  # in rows - information about id, sample name, mouse, tissue, tetramer, cell lineage
  # columns each sample
metadata_table <- function(input_data_01,sample_names){
  df_unique <- input_data_01 %>%
    distinct(Sample.bio.name, .keep_all = TRUE) # collect unique samples with metadata info
  df <- data.frame(ID = c(1:length(sample_names)),df_unique)

  # transposition of the metadata table so it can be easily attached to the big
  # count table by the last 6 rows
  dfT <- t(df)
  dfT <- data.frame(dfT)
  dfT <- rownames_to_column(dfT,"aa_seq")
  colnames(dfT)[2:ncol(dfT)] <- sample_names
  return(dfT)
}

# creating template table filled with zeros where columns are samples and
# rows are sequences
# input 
  # count_table = singular column with aa seq from the original data
  # sample_names = list of sample names
# output - template table filled with 0
create_countTable_layout <- function(count_table,sample_names){
  sorted_aa <- data.frame(aa_seq = sort(count_table[[1]])) # sorting sequences
  count_Table1 <- unique(sorted_aa)
  for (name in sample_names) {
    # the template for a count table with only zeros
    # in the next step, we put there actual umi numbers or binary values
    count_Table1 <- mutate(count_Table1, !!name := 0)
  }
  return(count_Table1)
}

# fills the template table with 0 or 1 or with UMI values
# input:
  # sample_name = list of sample names
  # countTable = template of a count table generated with create_countTable_layout
  # umiValue = true/false variable which determinates if the table will be binary
    # or contain umi values 
  # input_data_01 = datatable from current dataset edited in 01 step
# output - count table where in rows are aa seq, in columns each sample
  # filled with either binary or UMI values
umiCountSample <- function(sample_name, countTable, umiValue, input_data_01) {
  sample_dt <- as.data.table(input_data_01)
  countTable <- as.data.table(countTable)
  # select only sample names, umi count and aa seq
  sample_table <- sample_dt[`Sample.bio.name` %in% sample_name, 
                            .(`CDR3.amino.acid.sequence`, `Umi.count`)]
  # if umiValue is true -> generate table with UMI values
  if (umiValue) {
    # for each aa seq fill in UMI in each sample
    # the absent seq in sample will remain as 0
    countTable[sample_table, (sample_name) := i.Umi.count, 
               on = .(aa_seq = CDR3.amino.acid.sequence)]
  # umiValue is false -> generate binary val table
    } else {
      # for each aa seq fill in presence (1) in each sample
      # the absent seq in sample will remain as 0
    countTable[sample_table, (sample_name) := 1L, 
               on = .(aa_seq = CDR3.amino.acid.sequence)]
  }
  return(countTable)
}



# preparing the count table across all the samples
#input:
# metadata - table prepared with function metadata_table
# count_table - changing table from template with zeros, slowly filled with 
# values column by column
# input:
  # metadata = table generated with metadata_table function
  # count_table = only singular column with aa seq
  # sample_names = list of samples
  # input_data_01 = datatable from current dataset edited in 01 step
# output - one count table with UMI vals, other with binary vals
generating_count_table <- function(metadata,count_table,sample_names,input_data_01){
  sample_names <- as.character(sample_names)
  count_Table1 <- create_countTable_layout(count_table,sample_names)
  count_Table2 = data.frame(count_Table1)
  for (sample in sample_names){
    # count table with umi values
    count_Table1 = umiCountSample(sample,count_Table1,TRUE,input_data_01)
    # count table with binary values
    count_Table2 = umiCountSample(sample,count_Table2,FALSE,input_data_01)
  }
  # make the numbers numeric (so we can calculate occurrence later)
  count_Table1 %>% mutate(across(c(2:ncol(count_Table1)),as.numeric))
  count_Table2 %>% mutate(across(c(2:ncol(count_Table2)),as.numeric))
  # output is both of the tables, presented as list
  return(list(count_Table1,count_Table2))
}


# final function that will calculate the occurrence of each seq across all the 
# samples and then gets rid of the sequences, that are present in just one sample
# input is a list with count table for spleen samples and metadata table
# output is the final table
occurrence_and_order <- function(countTableSpleen){
  # extract count table and metadata from the list
  countTableSpleen.a <- as.data.frame(countTableSpleen[[1]])
  metadata.spleen.a  <- as.data.frame(countTableSpleen[[2]])
  countTableSpleen <- list(countTableSpleen.a,metadata.spleen.a)
  # determinate number of samples (15 for alpha, 31 for beta)
  No_of_samples = ncol(countTableSpleen.a)-1
  countTableSpleen.a <- countTableSpleen.a %>% 
    # count occurrence of each seq across all the samples
    mutate(Occurrence = rowSums(countTableSpleen.a[,2:ncol(countTableSpleen.a)])) %>%
    dplyr::filter(Occurrence > 1) %>% # only seqs that are Occurring in multiple samples
    dplyr::filter(!grepl("\\*", aa_seq)) # get rid of the unproductive seqs
  # add occurrence column to the metadata and merge it later with count table
  metadata.spleen.a <- metadata.spleen.a %>% mutate(Occurrence = NA)
  countTableSpleen <- rbind(countTableSpleen.a,metadata.spleen.a)
  
  # determine which samples are CD4 and which CD8
  group <- as.character(metadata.spleen.a[4,])  
  cd4_samples <- which(group == "CD4")
  cd8_samples <- which(group == "CD8")
  
  countTableSpleen <- countTableSpleen %>% 
    # add columns of occurrence in CD4 and CD8 samples
    mutate(Occurr_CD4 = 0) %>%
    mutate(Occurr_CD8 = 0)
  
  sample_cols <- 2:(No_of_samples+1)
  sample_data <- cbind(0,countTableSpleen[, sample_cols]) # sample data table
  sample_data <- as.data.frame(lapply(sample_data, as.numeric))
  # count occurrence in the CD4 samples or CD8 samples
  countTableSpleen$Occurr_CD4 <- rowSums(sample_data[, cd4_samples] > 0)
  countTableSpleen$Occurr_CD8 <- rowSums(sample_data[, cd8_samples] > 0)
  
  # order the columns properly  
  cols <- colnames(countTableSpleen)
  non_sample_cols <- c("aa_seq", "Occurrence", "Occurr_CD4", "Occurr_CD8")
  sample_cols <- setdiff(cols, non_sample_cols)
  sample_groups <- group[1:length(sample_cols)] 
  # reorder the columns so the CD4 samples are together and CD8 together
  ordered_samples <- sample_cols[order(factor(sample_groups, levels = c("CD4", "CD8", NA)))]
  # set up the order 
  final_col_order <- c("aa_seq", ordered_samples, non_sample_cols[non_sample_cols != "aa_seq"])
  countTableSpleen <- countTableSpleen[, final_col_order]
  return(countTableSpleen) # final version of the count table
}

