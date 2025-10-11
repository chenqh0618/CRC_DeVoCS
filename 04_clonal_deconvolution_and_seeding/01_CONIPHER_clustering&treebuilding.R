pjdir <- "~/project/crc/"
setwd(pjdir)

.libPaths("/home/data/cqh/R/x86_64-pc-linux-gnu-library/4.4")
library(data.table)
library(tidyverse)
library(parallel)
library(CONIPHER)
source("~/software/script_tools/R_color.R")


snv_rds_path <- "~/project/crc/snv_calling/case1_snv.rds"
snv_annotated_path <- "~/project/crc/snv_calling_case1/res_annotate/annovar_.hg19_multianno.txt"
cnv_rds_path <- "~/project/crc/res_wgs_cnv_seqz/case1.rds"


plot_clone_tree <- function(res_dir,
                            case_name,
                            clone_color = NULL,
                            plot_clone_size = T,
                            log_scale_size = F,
                            scale_size_factor = 20,
                            plot_ccf = F){
  raw_dir <- setwd(res_dir)
  library(cloneMap)
  library(igraph)
  
  trees_file_path <- "allTrees.txt"
  clones_file_path <- "clusterInfo.txt"
  
  raw_tree_content <- readLines(trees_file_path) 
  idx_start = grep(raw_tree_content, pattern = "tree 1$") + 1
  idx_end  = grep(raw_tree_content, pattern = "tree 2$") - 1
  tab_tree <- fread(text = raw_tree_content[idx_start:idx_end]) %>% 
    mutate(across(everything(), as.character))
  tab_clone <- fread(clones_file_path) %>% 
    mutate(clusterID = as.character(clusterID))
  
  clone_list <- as.vector(tab_tree) %>% 
    unlist() %>%  unique() #%>% c(tab_clone$clusterID) %>% unique() 
  
  if(is.null(clone_color)){
    library(qualpalr)
    get_colors <- function(color_num = color_num){
      pal <- qualpal(n = color_num, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))
      
      return(rownames(pal$HSL))
    }
    
    clone_color <- get_colors(color_num = length(clone_list))
    names(clone_color) <-  as.character(clone_list)
    clone_color <- clone_color[order(as.numeric(names(clone_color)))]
  }
  

  
  ## plot clone tree
  data_plot_tree <- graph_from_data_frame(d = tab_tree, directed = T)
  V(data_plot_tree)$color = clone_color[as.character(V(data_plot_tree)$name)]
  if(plot_clone_size){
    clone_size <- filter(tab_clone, clusterID %in% clone_list) %>% 
      distinct(clones = clusterID, size = nMuts) 
    clone_size <- clone_size$size %>% 
      `names<-`(clone_size$clones)
    if(log_scale_size){
      clone_size = log2(clone_size)
      clone_size <- clone_size/max(clone_size)*scale_size_factor
      }
    V(data_plot_tree)$size = clone_size[as.character(V(data_plot_tree)$name)]
    
  }
  pdf("clone_tree.pdf", width = 10, height = 10)
  plot(data_plot_tree, layout = layout_as_tree)
  graphics.off()
  
  
  ## plot clone CCF
  if(plot_ccf){
    tab_clone <- select(tab_clone, 
                        clones = clusterID,
                        sample = SAMPLE,
                        ccf = meanCCF) %>% 
      dplyr::filter(clones %in% clone_list) %>% 
      mutate(ccf = ccf/100)
    sample_list <- unique(tab_clone$sample) %>% {.[!is.na(.)]}
    
    paste0("clone_ccf.pdf") %>% pdf(width = 12, height = 12)
    if(case_name == "case1"){
      par(mfrow = c(4, 4))
    } else if(case_name == "case2"){
      par(mfrow = c(3, 3))
    }
    for(sample_name in sample_list){
      tab_ccf <- dplyr::filter(tab_clone, 
                               sample  %in% c(sample_name)) %>% 
        select(clones, CCF = ccf) %>% 
        dplyr::filter(CCF > 0) 
      
      clone_map_eg <- cloneMap(tree.mat = tab_tree, 
                               CCF.data = tab_ccf,
                               plot.data = F, 
                               output.Clone.map.obj = T)
      cloneMap(clone_map = clone_map_eg,
               clone.cols =  clone_color)  
      
      title(main = sample_name)
      
    }
    graphics.off()
  }

  ## finished 
  paste0("plots are generated! \n") %>% cat()
  setwd(raw_dir)
}



write_conipher_input <- function(case_name, 
                                 sample_method = "distance_thinned",  
                                 distance_thinned_length = 1e6,
                                 cluster_sampling_num = 3000,
                                 raw_result_dir,
                                 run_idx = 1){
  require(tidyverse)
  require(parallel)
  require(data.table)
  require(CONIPHER)
  # set.seed(run_idx)
  
  if(sample_method == "clusters"){
    temp_dir <- paste0(raw_result_dir, "/",case_name,"_clusterSampling", cluster_sampling_num, "_run", run_idx)
  } else if(sample_method == "distance_thinned"){
    temp_dir <- paste0(raw_result_dir, "/",case_name,"_distanceThined", distance_thinned_length, "_run", run_idx)
  } else {
    print("sample methods ony distance_thinned and clusters !")
  }
  
  dir.create(temp_dir)
  snv_rds_path <- paste0("~/project/crc/combined/snv_calling/", case_name, "_rescued_snv.rds")
  snv_annotated_path <- paste0("~/project/crc/combined_v4/snv_calling_", case_name, "/res_maf_0.05/res_annotate/annovar_.hg19_multianno.txt")
  cnv_rds_path <- paste0("~/project/crc/combined/res/res_wgs_cnv_seqz_v2/", case_name, ".rds")
  
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
  
  
  ## combined datas
  seg_data_temp <- dplyr::select(seg_data,
                                 sample, chrom = chr, from, to, major = Major, minor)
  
  data_combined_raw <- left_join(snv_data, seg_data_temp,
                                 by = c("sample", "chrom"),
                                 relationship = "many-to-many")  %>%
    dplyr::filter(pos >= from, pos <= to)
  
  ## filter no CNV info SNV
  list_snv_with_cnv <- count(data_combined_raw, snv_id) %>%
    dplyr::filter(n == length(unique(snv_data$sample))) %>%
    pull(snv_id)
  length(list_snv_with_cnv)
  list_snv_high_impact <- intersect(list_snv_with_cnv, 
                                    data_high_impact$snv_id)
  data_with_cnv <- filter(data_combined_raw, 
                          snv_id %in% list_snv_with_cnv)

  
  ## downsampling to fit 
  if(sample_method == "distance_thinned"){
    ## thined by distance
    snv_thined_list <- select(data_with_cnv, snv_id, chrom, pos) %>%
      distinct() %>%
      mutate(label = floor(pos/distance_thinned_length)) %>%
      slice_sample(n = 1, by = c(chrom, label)) %>%
      pull(snv_id)
    
    snv_keep_list <- union(snv_thined_list, list_snv_high_impact) 
    
    data_selected <- dplyr::filter(data_combined_raw, snv_id %in% snv_keep_list) %>%
      left_join(purity_data, by = c("sample"))
    
  } else if (sample_method == "clusters") {
    ## thined by clusters
    sample_tab <- data.frame(sample = unique(data_with_cnv$sample)) %>% 
      mutate(idx = 2^(1:nrow(.)))
    
    snv_tab <- filter(data_with_cnv, vaf >= 0.05) %>% 
      select(snv_id, sample, vaf) %>% 
      left_join(sample_tab) %>% 
      reframe(idx = sum(idx), .by = c(snv_id))
    snv_clusters <- count(snv_tab, idx) %>% 
      filter(n >= 5) %>% 
      mutate(prop = n/sum(n),
             sample_n = round(cluster_sampling_num * prop)) %>% 
      filter(sample_n > 0)
    sum(snv_clusters$sample_n)
    snv_select <- left_join(snv_tab, snv_clusters, by = c("idx")) %>% 
      rename(cluster_idx = idx) %>% 
      na.omit()  %>% 
      mutate(snv_idx = runif(nrow(.), 1, nrow(.))) %>% 
      mutate(rank = rank(snv_idx), .by = c(cluster_idx)) %>% 
      arrange(cluster_idx, snv_idx) %>% 
      filter(rank <= sample_n)
    
    snv_keep_list <- union(snv_select$snv_id, list_snv_high_impact) 
    
    data_selected <- dplyr::filter(data_with_cnv, 
                                   snv_id %in% snv_keep_list) %>% 
      left_join(purity_data, by = c("sample"))
    
  } else {
    print("ony distance_thinned and clusters !")
  }
  
  
  
  ## generated result input
  data_res <- select(data_selected,
                     CASE_ID = sample, #fix
                     SAMPLE = sample,
                     CHR = chrom,
                     POS = pos,
                     REF = ref,
                     ALT = alt,
                     REF_COUNT = rd,
                     VAR_COUNT = ad,
                     DEPTH = rd, #fix
                     COPY_NUMBER_A = major,
                     COPY_NUMBER_B = minor,
                     ACF = purity,
                     PLOIDY = ploidy) %>%
    mutate(CASE_ID = case_name,
           SAMPLE = str_replace(SAMPLE, "HCC", "LM") %>%
             str_replace("_", "."),
           CHR = str_remove(CHR, "chr"),
           ACF = round(ACF, 3),
           PLOIDY = round(PLOIDY, 3),
           # COPY_NUMBER_A = 1,
           # COPY_NUMBER_B  = 1,
           # PLOIDY = 2,
           DEPTH = REF_COUNT + VAR_COUNT) %>%
    mutate(SAMPLE = paste0(CASE_ID, "_", SAMPLE))
  
  fwrite(data_res, paste0(temp_dir, "/clustering_input.tab"),
         row.names = F, quote = F, sep = "\t")
  
  
}




run_conipher_tree <- function(case_name, 
                              distance_thinned = 0,
                              raw_result_dir,
                              run_idx = 1,
                              max_thread = 4){
  require(tidyverse)
  require(parallel)
  require(data.table)
  require(CONIPHER)
  
  temp_dir <- paste0(raw_result_dir, "/", case_name, "_thinned", distance_thinned, "_run", run_idx) 
  # run clustering 
  conipher_clustering(case_id = case_name,
                      prefix = case_name,
                      out_dir = paste0(temp_dir, "/Clustering/"),
                      input_tsv_loc = paste0(temp_dir, "/clustering_input.tab"),
                      nProcs = max_thread,
                      PyClone_path = "/home/data/cqh/miniforge3/envs/conipher/bin/PyClone"
  )
  
  # build tree
  conipher_treebuilding(input_tsv_loc = paste0(temp_dir, "/Clustering/", case_name, ".SCoutput.CLEAN.tsv"),
                        out_dir = paste0(temp_dir, "/res_trees/"),
                        prefix = case_name,
                        ccf_buffer = 10,
                        pval_cutoff = 0.01,
                        use_boot = TRUE,
                        merge_clusters = TRUE,
                        correct_cpn_clusters = TRUE,
                        adjust_noisy_clusters = TRUE,
                        adjust_noisy_clusters_prop = 0.05,
                        min_ccf = 0.05,
                        min_cluster_size = 20,
                        multi_trees = T,
                        nProcs = max_thread
  )
  
  # plot clone CCF  
  try({
	plot_clone_tree(res_dir = paste0(temp_dir, "/res_trees/"), 
                    case_name = case_name)
   })
  

}



# RUN PROCESS ------------------------------------------------------------------
raw_result_dir <- "res/res_conipher/"; dir.create(raw_result_dir)

test_tab <- expand.grid(case = c("case1", "case2"),
                        run = 1:3,
                        stringsAsFactors = F) %>%
  arrange(desc(distance))


# GET INPUT
mclapply(1:nrow(test_tab), mc.cores = 2, function(i){
  run_idx <- test_tab$run[i]
  case_name <- test_tab$case[i]
  distance_thinned <- 1e6

  paste0(case_name, "_", distance_thinned, "thined_run", run_idx, " starts processing! \n") %>% cat()
  write_conipher_input(raw_result_dir = raw_result_dir,
                       case_name = case_name,
                       distance_thinned = distance_thinned,
                        run_idx = run_idx)

})

# PROCESS CONIPHER
mclapply(1:nrow(test_tab), mc.cores = 6, function(i){
  Sys.sleep(i*30-30)
  run_idx <- test_tab$run[i]
  case_name <- test_tab$case[i]
  distance_thinned <- test_tab$distance[i]

  paste0(case_name, "_", distance_thinned, "thined_run", run_idx, " starts processing! \n") %>% cat()
  run_conipher_tree(raw_result_dir = raw_result_dir,
                    case_name = case_name,
                    distance_thinned = distance_thinned,
                    run_idx = run_idx,
                    max_thread = 36)

})


