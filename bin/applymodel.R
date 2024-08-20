#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-05-15
# Applies ranger model to read list

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(ranger))

#args <- data.frame(reads = "/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/frameAnnotation/ppiatrans.feather",
#           model = "/hpc/hub_oudenaarden/mvanins/RPF/analysis/mlrtest/latest.tuned.model.RData",stringsAsFactors = FALSE)

option_list = list(
  make_option(c("-r", "--reads"), type = "character", default = NULL, help = "feather file containing read information dataframe", metavar="character"),
  make_option(c("-m", "--model"), type = "character", default = NULL, help = "RData file containing tuned.model", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)

if(is.null(args$reads) | is.null(args$model)){
  print_help(opt_parser)
  stop("Please supply reads and model", call.= FALSE)
}


# Function to process data in chunks
process_in_chunks <- function(file, model, chunk_size = 3e6) {
  load(model)
  if (!exists("tuned.model")) {
    stop("tuned.model not found in model file")
  }
  
  if (grepl("feather", file)) {
    dt <- as.data.table(read_feather(file))
  } else {
    dt <- fread(file)
  }
  
  total_rows <- nrow(dt)
  chunks <- ceiling(total_rows / chunk_size)
  predictions <- list()
  
  for (i in 1:chunks) {
    start <- ((i - 1) * chunk_size) + 1
    end <- min(i * chunk_size, total_rows)
    
    chunk <- dt[start:end, ]
    
    predReads <- chunk %>%
      select(length, 
             base5, base5m1, base5m2, base5m3, base5m4, base5m5, base5m6, base5m7, base5m8, 
             base5p1, base5p2, base5p3, base5p4, base5p5, base5p6, base5p7, base5p8, 
             base3, base3m1, base3m2, base3m3, base3m4, base3m5, base3m6, base3m7, base3m8, 
             base3p1, base3p2, base3p3, base3p4, base3p5, base3p6, base3p7, base3p8) %>%
      mutate_if(is.character, as.factor)
    
    pred <- predict(tuned.model$learner.model, predReads, num.threads = 10)
    
    chunk <- chunk %>%
      mutate(offset = as.numeric(levels(pred$predictions))[pred$predictions]) %>%
      mutate(psite = cut5 - (offset + 3))
    
    predictions[[i]] <- chunk
    rm(predReads, chunk)
    gc()
  }
  
  result <- rbindlist(predictions)
  result
}

predictedReads <- process_in_chunks(args$reads, args$model)

if (grepl("feather", args$reads)) {
  write_feather(predictedReads, paste0(tools::file_path_sans_ext(basename(args$reads)), ".predicted.feather"))
} else {
  fwrite(predictedReads, paste0(tools::file_path_sans_ext(basename(args$reads)), ".predicted.csv"))
}
