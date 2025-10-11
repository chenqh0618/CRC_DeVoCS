.libPaths("/home/data/cqh/R/x86_64-pc-linux-gnu-library/4.4")
wkdir <- "~/project/crc/combined//"
setwd(wkdir)

library(tidyverse)
library(data.table)
library(parallel)
library(OptM)
source("~/software/script_tools/treemix_plotting_funcs.R")


get_treemix_input_CCF <- function(data_path = data_path,
                              run_name = "case1_test", 
                              object_type = "sample",
                              res_dir = "res_treemix/res_case1", 
                              if_with_normal = T,
                              write_input_file = F){
  require(tidyverse)
  require(data.table)
  if (grepl("pyclone", data_path)){
    method = "pycloneCCF"
  } else if(grepl("conipher", data_path)){
    method = "conipherCCF"
  }
  run_name <- paste0(run_name, "_", method)
  

	data <- fread(data_path)

	data_input <- reframe(data, 
						  sample = str_remove(SAMPLE, "case._"),
						  vaf = CCF_OBS  * ccf_scale_factor, 
						  ad = round(600 * vaf),
						  rd = 600 - ad,
						  chrom = paste0("chr",CHR),
						  pos = POS,
						  site = paste0(CHR, "_", pos))  

  
  temp_ad <- pivot_wider(data_input,
                         id_cols = c(site), 
                         names_from = sample, 
                         values_from = ad, 
                         values_fill = 0) %>% 
    pivot_longer(cols = -site, values_to = "ad", names_to = "sample")
  temp_rd <- pivot_wider(data_input,
                         id_cols = c(site), 
                         names_from = sample, 
                         values_from = rd, 
                         values_fill = 600) %>% 
    pivot_longer(cols = -site, values_to = "rd", names_to = "sample")
  
  data_input <- left_join(temp_ad, temp_rd, by = c("sample", "site"))
  
  ## region or sample 
  if(object_type == "region"){
    frq_res <- mutate(data_input, 
                      region = case_when(sample %in% c("C0", "C1", "C2") ~ "PI",
                                         sample %in% c("C10", "C12", "C14", "CB") ~ "PO",
                                         sample %in% c("HCC1", "LM1") ~ "LM1",
                                         sample %in% c("HCC2", "LM2") ~ "LM2",
                                         sample %in% c("HCC3", "LM3") ~ "LM3",
                                         grepl("HCC4", sample) ~ "LM4",
                                         grepl("LM4", sample) ~ "LM4",
                                         grepl("LU", sample) ~ "LU",
                                         grepl("LD", sample) ~ "LD",
                                         grepl("RU", sample) ~ "RU",
                                         grepl("RD", sample) ~ "RD",
                                         # grepl("LM", sample) ~ "LM",
                                         .default = sample) %>% 
                        factor(levels = c( "PO", "PI", 
                                           "LU", "LD", "RU", "RD",
                                           "LM1", "LM2", "LM3", "LM4",
                                           "LM"))) %>% 
      group_by(site, region) %>% 
      reframe(
        ad = sum(ad), rd = sum(rd), #raw data
      ) %>% 
      rename(sample = region)
  } else if(object_type == "sample"){
    frq_res <- select(data_input, site, sample, ad, rd)
  }
  
  
  ## if with normal sample 
  if(if_with_normal){
    frq_normal <- reframe(frq_res, 
                          ad = 0, 
                          rd = 600,
                          .by = c(site)) %>% 
      mutate(sample = "NORMAL")
    frq_res <- bind_rows(frq_res, frq_normal)
  }
  
  frq_res <- mutate(frq_res, idx = paste0(ad, ",", rd)) %>% 
    select(site, sample, idx) %>% 
    arrange(site, desc(sample)) %>% 
    pivot_wider(names_from = sample, 
                values_from = idx,
                values_fill = "0,600") %>% 
    column_to_rownames("site")  
  
  if(write_input_file){
    run_dir <- paste0(res_dir, "/", run_name)
    dir.create(run_dir, recursive = T)
    rawdir <- setwd(run_dir)
    
    writeLines(colnames(frq_res), "pop_list") 
    
    fwrite(frq_res, 
           file = "input.allel.frq",
           sep = "\t", quote = F, row.names = F)
    paste0("gzip -f input.allel.frq") %>% system()
    
    setwd(rawdir)
  }
  return(frq_res)
}

run_treemix <- function(run_dir = NULL, 
                        max_mig = 10,
                        max_rep = 5,
                        k_num = 50,
                        num_cores = 30,
                        se = F,
                        global = F,
                        noss = T,
                        bootstrap = T,
                        root_name = "NORMAL",
                        input_file = "input.allel.frq.gz",
                        pop_file = "pop_list"){
  require(tidyverse)
  require(parallel)
  require(OptM)
  source("~/software/script_tools/treemix_plotting_funcs.R")
  
  if(dir.exists(paste0(run_dir, "/res"))){return(NULL)}
  if(dir.exists(paste0(run_dir, "/process"))){file.remove(paste0(run_dir, "/process/.*"))}
  
  setwd(run_dir)
  process_dir <- paste0(run_dir, "/process"); dir.create(process_dir)
  file.copy(input_file, paste0(process_dir, "/", input_file))
  file.copy(pop_file, paste0(process_dir, "/", pop_file))
  
  setwd(process_dir)
  mig_tab <- expand.grid(mig = 0:max_mig, rep = 1:max_rep)
  
  mclapply(1:nrow(mig_tab), mc.cores = num_cores, function(i){
    mig_num = mig_tab$mig[i]
    rep_num = mig_tab$rep[i]
    paste0("source ~/.bashrc; ",
           "source /home/data/cqh/miniforge3/etc/profile.d/conda.sh; ", 
           "conda activate treemix;",
           "treemix -i ", input_file,
           " -m ", mig_num,
           if_else(se, " -se ", ""),
           if_else(global, " -global ", ""),
           if_else(noss, " -noss ", ""),
           if_else(bootstrap, " -bootstrap ", ""),
           " -k ", k_num,
           " -root ", root_name,
           " -seed ", 2^rep_num-2,
           " -o mig", mig_num, "_rep", rep_num,
           " &> process_mig", mig_num, "_rep", rep_num, ".log") %>% 
      system()
  })
  
  
  res_dir <- paste0(run_dir, "/res"); dir.create(res_dir)
  # best mig_num  ----
  try({
    graphics.off()
    pdf(paste0(res_dir, "/optM_liner_model.pdf"))
    linear=optM(folder = ".",
                method = "linear")
    plot_optM(linear, method = "linear")
    dev.off()
  })
  
  try({
    graphics.off()
    pdf(paste0(res_dir, "/optM_Evanno_model.pdf"))
    evanno=optM(folder = ".",
                method = "Evanno")
    par(mfrow=c(2,1))
    plot_optM(evanno, method = "Evanno")
    dev.off()
  })
  
  if(file.exists(paste0(process_dir,"/Rplots.pdf"))){
    file.copy(from = paste0(process_dir,"/Rplots.pdf"), 
              to = paste0(res_dir, "/optM_Evanno_model2.pdf"),
              overwrite = T)
    }
  
  
  
  # phylogeny net ---- 
  lapply(0:max_mig, function(mig_num){
    lapply(1:max_rep, function(rep_num){
      stem_name <- paste0("mig", mig_num, "_rep", rep_num)
      
      try({
        graphics.off()
        paste0(res_dir, "/res_",stem_name, "_tree.pdf") %>% 
          pdf(width = 16, height = 16)
        plot_tree(stem = stem_name, ,
                  disp = 0.0001,
                  plus = 0.0008,
                  cex = 2, lwd = 4, font = 1,
                  scale = T,
                  arrow = 0.25,
                  ybar = 0.15,
                  plotmig = T, mbar = T, plotnames = T,
                  xmin = 0)
        dev.off()
      })
      
      try({
        graphics.off()
        paste0(res_dir, "/res_",stem_name, "_SE.pdf") %>% 
          pdf(width = 16, height = 16)
        plot_resid(stem = stem_name, "pop_list")
        dev.off()
      })
      
    })
  })
}

treemix_plot <- function(res_dir, stem_name){
  rawdir <- setwd(res_dir)
  require(corrplot)
  # stem_name <- paste0("mig", mig_num, "_rep", rep_num)

  plot_tree <- plot_tree(stem = stem_name,
                         disp = 0.0001,
                          plus = 0.0008,
                          cex = 1.6, lwd = 4, font = 1,
                          scale = F,
                          arrow = 0.25,
                         ybar = 0.15,
                         plotmig = T,
                         mbar = T,
                         plotnames = T,
                          xmin = 0)

  temp <- plot_resid(stem = stem_name, "pop_list")
  rownames(temp) <- colnames(temp)
  temp <- as.matrix(temp)
  se_value <- sd(temp) / (sqrt(nrow(temp))-1)
  temp <- temp/se_value


  corrplot::corrplot(corr = temp,
                     method = "color",
                       type = "lower",
                       is.corr = F,
                       col = c(rev(Blues), "white", Reds),
                       outline = "black",
                       tl.pos = "ld",
                       tl.cex = 1.2,
                       tl.col = "black",
                       tl.offset = 0.7,
                       tl.srt = 0,
                       cl.pos = "r",
                       cl.length = 5,
                       cl.cex = 1,
                       cl.ratio = 0.2,
                       number.digits = 4,
                       addshade = "negative" )

    pheatmap::pheatmap(temp, clustering_method = "ward.D2", color = c(rev(Blues), "white", Reds))

}



# get treemix input -------
test_tab <- expand.grid(case = c("case1", "case2"),
                        run = 1:3,
                        stringsAsFactors = F)
mclapply(1:nrow(test_tab), mc.cores = 6, function(i){
  run_idx = test_tab$run[i]
  case_name = test_tab$case[i]

  data_path <- paste0(res_conipher_path, "/", case_name, "/run", run_idx,"/Clustering/", case_name, ".SCoutput.CLEAN.tsv")

  get_treemix_input_directCCF(data_path = data_path,
                    run_name = paste0(case_name, "_run", run_idx),
                    object_type = object_type,
                    res_dir = res_dir,
                    write_input_file = T)
  })


# run treemix
dir_list <- list.files(res_dir, 
                       full.names = T,
                       pattern = "root")

lapply(dir_list, function(dir_name){
  root_name <- str_extract(dir_name, "(?<=root).*?(?=_)")
  run_treemix(run_dir = dir_name, 
              max_mig = 10,
              max_rep = 5,
              k_num = 50,
              num_cores = 30,
              se = T,
              global = T,
              noss = T,
              bootstrap = T,
              root_name = root_name,
              input_file = "input.allel.frq.gz",
              pop_file = "pop_list")
})

