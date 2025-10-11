.libPaths("/home/data/cqh/R/x86_64-pc-linux-gnu-library/4.4")
wkdir <- "~/project/crc/wgs_process/cnv_calling/"
setwd(wkdir)

library(sequenza)
library(CNAqc)
library(tidyverse)
library(data.table)
library(parallel)

seqz_file_dir <- "/public/data/cqh_project/crc/wgs_process/res_seqz/"
sample_list <- readLines("/public/data/cqh_project/crc/data/tumor_list")


mclapply(sample_list, mc.cores = 9, function(sample_name){
	try({
	  res <- Sequenza_CNAqc(
		sample_id = sample_name, 
		seqz_file = paste0(seqz_file_dir, "/", sample_name, ".seqz.gz"),
		mutations = NULL,
		sex = "M",
		cellularity = c(0.05, 1),
		ploidy = c(1.0, 8.0),
		reference = "GRCh37",
		normalization.method = "median",
		window = 1e5,
		gamma = 40,
		kmin = 10,
		min.reads = 30,
		min.reads.baf = 30,
		min.reads.normal = 30,
		max.mut.types = 1,
		delta_cellularity = 0.01,
		delta_ploidy = 0.05,
		verbose = T)  
	})
	paste0(Sys.time(), ": ", sample_name, "parsing is done !") %>% cat()
})






