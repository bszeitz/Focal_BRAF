#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# USAGE: Rscript Parse_PepQuery_outputs.R <folder with PepQuery results> <filename for export>
# for example: Rscript Parse_PepQuery_outputs.R E:/MM_Pilot/foc_braf/pepquery PepQuery_results

library(seqinr)


########################################################################
# Function definitions

# Function to select the reference / modified peptide with the highest score for a given spectrum returned by PepQuery
select_best_match <- function(df){
  library(dplyr)
  df$Row_name <- seq(1,nrow(df))
  df <- df %>%
    group_by(spectrum_title) %>%
    arrange(desc(score)) %>%
    slice(1) %>% # selects the first row of each group (spectrum_title) after the data has been arranged based on the score in descending order
    ungroup() %>% # removes the grouping
    select(-Row_name) # removes the Row_name column as it's no longer needed
  
  return(df)
}
# Function to create an empty df formatted according to PepQuery output
create_empty_df <- function(mgf_filename,known_or_novel){
 return(data.frame(file = mgf_filename,
                   input_peptide_type  =  known_or_novel,
                   spectrum_title = NA, 
                   input_peptide = NA, 
                   validation_info = NA, 
                   ref_peptide = NA, 
                   ref_modification = NA, 
                   modref_peptide = NA, 
                   modref_modification = NA,
                   input_modification = NA, 
                   input_n = NA, 
                   input_charge = NA, 
                   input_exp_mass = NA, 
                   input_tol_ppm = NA, 
                   input_tol_da = NA,
                   input_isotope_error = NA, 
                   input_pep_mass = NA, 
                   input_mz = NA, 
                   input_score = NA, 
                   input_n_db = NA,
                   input_total_db = NA, 
                   input_n_random = NA, 
                   input_total_random = NA, 
                   input_pvalue = NA, 
                   input_rank = NA, 
                   input_n_ptm = NA,
                   input_confident = NA,
                   input_ref_delta_score = NA, 
                   input_mod_delta_score = NA, 
                   ref_exp_mass = NA,
                   ref_pep_mass = NA,
                   ref_tol_ppm = NA, 
                   ref_tol_da = NA, 
                   ref_isotope_error = NA, 
                   ref_score = NA,
                   modref_charge = NA,
                   modref_exp_mass = NA, 
                   modref_pep_mass = NA, 
                   modref_tol_ppm = NA,
                   modref_tol_da = NA,
                   modref_isotope_error = NA,
                   modref_score = NA
 ))
}



#foldername <- folders[2]
summarize_PepQ_res <- function(foldername){
  paramfile <- read.delim(paste0(foldername,"/parameter.txt"), header = F)
  paramfile <- as.matrix(apply(paramfile, 1, function(Y){gsub("/mnt/e/","E:/",Y)}))
  
  commands_splitted <- strsplit(paramfile[2,1], split=" -", fixed=T)[[1]]
  mgf_file <- gsub("ms ","",commands_splitted[startsWith(commands_splitted,"ms ")])
  mgf_filename <- strsplit(mgf_file, split="/", fixed=T)[[1]]
  mgf_filename <- gsub(".mgf", "",mgf_filename[grep(".mgf",mgf_filename, fixed=T)])
  known_or_novel <- ifelse(gsub("s ","",commands_splitted[startsWith(commands_splitted,"s ")]) =="2", "reference","non-reference")

  ###########################################################
  if (!file.exists(paste0(foldername,"/psm_rank.txt"))){
    print(paste0(foldername, " is done - ", Sys.time()))
    print(grep(foldername, folders, fixed = T))
    return(create_empty_df(mgf_filename,known_or_novel))
  }
  
  if (class(try(read.delim(paste0(foldername,"/psm_rank.txt")), silent = T))=='try-error') {
    print(paste0(foldername, " is done - ", Sys.time()))
    print(grep(foldername, folders, fixed = T))
    return(create_empty_df(mgf_filename,known_or_novel))
  }
  
  
  psm_rank <- read.delim(paste0(foldername,"/psm_rank.txt"))
  colnames(psm_rank)[-grep("spectrum_title",colnames(psm_rank))] <- paste("input",colnames(psm_rank)[-grep("spectrum_title",colnames(psm_rank))],sep="_")
  
  psm_rank$validation_info <- sapply(seq(1,nrow(psm_rank)), function(i){
    if (psm_rank$input_n_db[i] !=0){
      y <- "better_ref_pep"
    } else if (psm_rank$input_confident[i] =="Yes"){
      y <- "confident"
    } else if (psm_rank$input_n_ptm[i] ==-1){
      y <- "high_p_value"
    } else if (psm_rank$input_n_ptm[i] >0){
      y <- "better_ref_pep_with_mod"
    } else {
      y <- "ERROR?"
    }
    return(y)
  })
  
  
  
  ###########################################################
  # Load ref.pept.search results and paralelly, select the best matching peptide
  detail <- select_best_match(read.delim(paste0(foldername,"/detail.txt")))
  colnames(detail)[2:ncol(detail)] <- paste("ref",colnames(detail)[2:ncol(detail)],sep="_")
  all_results <- merge(psm_rank, detail, by="spectrum_title", all.x=T)
  
  
  ###########################################################
  # Load mod.ref.pept.search results and paralelly, select the best matching peptide
  if (file.exists(paste0(foldername,"/ptm.txt"))){
    ptm <- select_best_match(read.delim(paste0(foldername,"/ptm.txt")))
  } else {
    ptm <- data.frame(spectrum_title=NA,
                      peptide=NA, charge=NA, exp_mass=NA, pep_mass=NA, tol_ppm=NA, tol_da=NA, isotope_error=NA, modification=NA, score=NA)
  }
  colnames(ptm)[2:ncol(ptm)] <- paste("modref",colnames(ptm)[2:ncol(ptm)],sep="_")
  all_results <- merge(all_results, ptm, by="spectrum_title", all.x=T)
  
  all_results[is.na(all_results)] <- ""
  
  first_columns <- c("spectrum_title", "charge","input_peptide", "input_peptide_modification","validation_info", 
                     "ref_peptide", "ref_modification","modref_peptide", "modref_modification")

  all_results_reordered <- cbind(mgf_filename, known_or_novel,all_results[,colnames(all_results) %in% first_columns], 
                                 all_results[,!colnames(all_results) %in% first_columns])
  colnames(all_results_reordered)[1:2] <- c("file","input_peptide_type")
  all_results_reordered[is.na(all_results_reordered)] <- ""
  print(paste0(foldername, " is done - ", Sys.time()))
  print(grep(foldername, folders, fixed = T))
  return(all_results_reordered)
}

setwd(args[1])
folders <- gsub("parameter.txt","",list.files(pattern="parameter.txt", recursive = T))


result_list <- do.call(rbind,lapply(folders, summarize_PepQ_res))

write.table(result_list, 
            paste0(args[2],".txt"),
            sep="\t", row.names = F, quote = F, col.names = T)
