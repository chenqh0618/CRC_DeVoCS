pjdir <- "~/project/crc/combined/"
setwd(pjdir)

.libPaths("/home/data/cqh/R/x86_64-pc-linux-gnu-library/4.4")
library(data.table)
library(tidyverse)
library(parallel)



rewrite_clonefinder_fasta <- function(res_clonefinder_dir){
  input_snv_path <- paste0(res_clonefinder_dir, "/input_snv")
  keywd <- list.files(res_clonefinder_dir, pattern = ".meg$", full.names = T) %>% 
    basename() %>% str_remove(pattern = ".meg")
  raw_fa_path <- paste0(res_clonefinder_dir, "/", keywd, ".fa")
  
  
  site_info_raw <- fread(input_snv_path)
  site_info <- select(site_info_raw, ref=Wild, alt=Var) %>% 
    mutate(ref_len = nchar(ref), alt_len = nchar(alt),
           seq_len = if_else(ref_len>alt_len, ref_len, alt_len)) %>% 
    mutate(ref_new = paste0(ref, str_dup("-", seq_len - ref_len)),
           alt_new = paste0(alt, str_dup("-", seq_len - alt_len))) %>% 
    select(ref = ref_new, alt = alt_new)
  
  clone_gt <- ape::read.FASTA(raw_fa_path)
  clone_gt <- as.matrix(clone_gt)
  clone_snv <- clone_gt
  clone_snv <- apply(clone_snv, 1, function(x){
    x <- as.vector(x) %>% as.character()
    x[x=="88"] = 0
    x[x=="18"] = 1
    return(x)
  })
  write.table(clone_snv, paste0(res_clonefinder_dir, "/", keywd, ".gt"), 
              quote = F,
              row.names = F,
              sep = "\t")
  
  
  temp <- apply(clone_gt, 1, function(x){
    x <- as.vector(x) %>% as.character()
    x[x=="88"] = site_info$ref[x=="88"]
    x[x=="18"] = site_info$alt[x=="18"]
    return(x)
  })
  res <- lapply(colnames(temp), function(seq_name){
    c(paste0(">", seq_name),
      paste0(temp[,seq_name], collapse = ""))
  }) %>% 
    do.call(what = c)
  
  length(res)
  writeLines(res, paste0(res_clonefinder_dir, "/final_", keywd, ".fa"))
}

run_clonefinder <- function(case_name,
                            distance_thinned_length = 1e5,
                            run_idx,
                            raw_result_dir){

  keywd <- paste0(case_name,"_distanceThined", distance_thinned_length, "_run", run_idx)

  res_dir <- paste0(raw_result_dir, "/res_", keywd)
  dir.create(res_dir)
  
  snv_rds_path <- paste0("~/project/crc/snv_calling/", case_name,"_snv.rds")
  snv_annotated_path <- paste0("~/project/crc/snv_calling_", 
                               case_name,
                               "/res_maf_0.05/res_annotate/annovar_.hg19_multianno.txt")
  cnv_rds_path <- paste0("~/project/crc/res_wgs_cnv_seqz/", 
                         case_name,".rds")
  
  
  # input data
  snv_data <- read_rds(snv_rds_path) %>%
    mutate(chrom = paste0("chr", chrom),
           snv_id = paste0(chrom, "_", pos))
  
  cnv_data <- read_rds(cnv_rds_path)
  
  seg_data <- lapply(1:nrow(cnv_data), function(i){
    cnv_data$sequenza[[i]]$segments %>%
      mutate(sample = cnv_data$sample[[i]]) %>%
      na.omit()
  }) %>%
    do.call(what = rbind)
  
  purity_data <- select(cnv_data, sample, purity, ploidy)
  
  data_high_impact <- fread(snv_annotated_path) %>%
    dplyr::filter(ExonicFunc.refGene %in% c("nonsynonymous SNV",
                                            "stopgain",
                                            "stoploss",
                                            "startloss")) %>%
    mutate(snv_id = paste0("chr",Chr, "_", Start))
  
  
  # combined datas
  seg_data_temp <- dplyr::select(seg_data,
                                 sample, chrom = chr, from, to, major = Major, minor)
  
  data_combined_raw <- left_join(snv_data, seg_data_temp,
                                 by = c("sample", "chrom"),
                                 relationship = "many-to-many")  %>%
    dplyr::filter(pos >= from, pos <= to)
  
  ## filter no CNV info SNV
  
  snv_noCNA_list <- filter(data_combined_raw, 
                           major == 1, 
                           minor == 1) %>% 
    count(snv_id) %>% 
    dplyr::filter(n == length(unique(snv_data$sample))) %>% 
    pull(snv_id)
  
  snv_noCNA_IMPACT_list <- intersect(snv_noCNA_list,
                                     data_high_impact$snv_id)
  
  length(snv_noCNA_list)

  data_cnv_neu <- dplyr::filter(data_combined_raw,
                             snv_id %in% snv_noCNA_list)  
  ## downsampling snv by distance
  snv_keep_list <- select(data_cnv_neu, snv_id, chrom, pos) %>%
    distinct() %>%
    mutate(label = floor(pos/distance_thinned_length)) %>%
    slice_sample(n=1, by = c(chrom, label)) %>%
    pull(snv_id)
  snv_keep_list <- c(snv_keep_list, snv_noCNA_IMPACT_list) %>% unique()
  
  data_selected <- dplyr::filter(data_cnv_neu,
                                 snv_id %in% snv_keep_list) %>%
    left_join(purity_data, by = c("sample")) 

  ## generated result input
  data_res <- select(data_selected, 
                     SNVID = snv_id, sample, 
                     Wild = ref, Var = alt, 
                     ref = rd, alt = ad) %>% 
    pivot_longer(cols = c(ref, alt), 
                 names_to = "type", 
                 values_to = "dp") %>% 
    pivot_wider(id_cols = c(SNVID, Wild, Var),
                names_from = c(sample, type), names_sep = ":",
                values_from = dp)
  
  fwrite(data_res, paste0(res_dir, "/", keywd, ".tab"),
         row.names = F, quote = F, sep = "\t")
  
  
  # run clonefinder model
  paste0("bash ~/software/script_tools/process_clonefinder.sh ",
         res_dir, "/", keywd, ".tab ",
         " 20 ",
         " 3 ",
         " 0.03 ",
         " &> ", res_dir, "/process_", keywd, ".log") %>% 
    system()
  
  
  rewrite_clonefinder_fasta(paste0(res_dir, "/res_clonefinder"))
  
  
  # run iqtree
  res_iqtree_dir <- paste0(res_dir, "/res_iqtree")
  dir.create(res_iqtree_dir)
  
  paste0("source ~/.bashrc; ",
         "source /home/data/cqh/miniforge3/etc/profile.d/conda.sh; ",
         "conda activate iqtree;",
         "iqtree2 ",
         " -s ", res_dir, "/res_clonefinder/final_", keywd, ".fa ",
         " -m MFP --mset beast2 ",
         " -o hg19 ",
         " -T 8 ",
         " -B 1000 ",
         " --prefix ", res_dir, "/res_iqtree/", keywd) %>%
    system()
  
  
  
  
  # plot 
  library(treedataverse)
  region_clone_ccf <- paste0(res_dir, "/res_clonefinder/", keywd, ".txt") %>% 
    fread() %>% 
    pivot_longer(cols = starts_with("Clo"),
                 names_to = "clone",
                 values_to = "ccf") %>%
    rename(sample = Tumor) %>% 
    mutate(tumor = case_when(grepl("LM|HCC", sample) ~ "metastasis",
                             .default = "primary"),
           region = case_when(sample %in% c("C0", "C1", "C2") ~ "PI",
                              sample %in% c("C10", "C12", "C14", "CB") ~ "PO",
                              grepl("LM4|HCC4", sample) ~ "LM4",
                              grepl("LM3|HCC3", sample) ~ "LM3",
                              grepl("LM2|HCC2", sample) ~ "LM2",
                              grepl("LM1|HCC1", sample) ~ "LM1",
                              .default = str_sub(sample, 1, 2)))  %>% 
    reframe(ccf = mean(ccf, na.rm = T), .by = c(clone, tumor, region)) %>% 
    pivot_wider(names_from = region, values_from = ccf, id_cols = clone) %>% 
    column_to_rownames("clone")
  
  
  clone_iqtree <- paste0(res_dir, "/res_iqtree/", keywd, ".treefile") %>% 
    read.iqtree()
  clone_idx=fortify(clone_iqtree) %>%
    subset(isTip) %>%
    arrange(y) %>%
    pull("label")
  
  p1 <- ggtree(clone_iqtree) +
    geom_tiplab(hjust = -0.02, align = T) +
    geom_nodepoint(aes(subset = UFboot > 60),
                   color = "grey") +
    geom_nodepoint(aes(subset = UFboot > 80),
                   color = "black")
  
  offset_value = p1$data$x %>% max()*0.1
  plot <- gheatmap(p1, region_clone_ccf,
                   offset = offset_value,
                   width = 0.6,
                   colnames_angle = 45,
                   legend_title = "Clone freq",
                   high="#D32F2F",
                   low = "#FFF7F3",
                   font.size = 4) +
    ggtitle(label = keywd)
  
  ggsave(plot = plot, paste0(res_dir, "/res_mltree_region.pdf"))
}





test_tab <- expand.grid(case_name = c("case1", "case2"),
                        run_idx = 1:3,
                        stringsAsFactors = F)


mclapply(1:nrow(test_tab), mc.cores = 10, function(i){
  case_name = test_tab$case_name[i]
  run_idx = test_tab$run_idx[i]
  distance_thinned_length = 1e5
  
  
  paste0(case_name, "_run", run_idx, " is runing! \n") %>% cat()
  
  
  run_clonefinder(case_name = case_name,
                  distance_thinned_length = distance_thinned_length,
                  run_idx = run_idx,
                  raw_result_dir = "~/project/crc/combined/res_wgs_deconvolution/res_clonefinder3")
})

