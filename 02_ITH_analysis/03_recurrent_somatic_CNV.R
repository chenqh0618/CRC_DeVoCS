#!/usr/bin/env Rscript

library(GenomicRanges)
library(data.table)
library(tidyverse)
library(parallel)


# subject data  -----------
all_lesions <- fread("data/database/BRCA.gistic/all_lesions.conf_99.txt",
                     header = TRUE, sep = "\t")

peak_tab <- str_split(all_lesions$`Peak Limits`, 
                      pattern = ":|-|\\(", 
                      simplify = T) %>% 
  {.[,1:3]} %>% 
  `colnames<-`(c("chr", "start", "end")) %>%
  as.data.frame() %>% 
  mutate(name = all_lesions$Descriptor,
         p_value = all_lesions$`Residual q values after removing segments shared with higher peaks`,
         start = as.numeric(start),
         end = as.numeric(end),
         type = case_when(grepl("Deletion", all_lesions$`Unique Name`) ~ "del",
                          grepl("Amplification", all_lesions$`Unique Name`) ~ "amp")) %>% 
  distinct() 


amp_peak = filter(peak_tab, type == "amp") %>% 
  mutate(cnv_id = paste0("amp_", 1:nrow(.)))
del_peak = filter(peak_tab, type == "del") %>% 
  mutate(cnv_id = paste0("del_", 1:nrow(.)))

peak_tab <- rbind(amp_peak, del_peak)

amp_gr <- GRanges(seqnames = amp_peak$chr, ranges = IRanges(start = amp_peak$start, end = amp_peak$end))
del_gr <- GRanges(seqnames = del_peak$chr, ranges = IRanges(start = del_peak$start, end = del_peak$end))



# query data ------------
get_hit_seg <- function(data, amp_gr, del_gr, case_name = NULL){
  hit_cnv <- lapply(data$sample, function(sample_temp){
    temp <- filter(data, sample == sample_temp)  
    temp <- temp$cnaqc[[1]]$cna %>% 
      mutate(type = dplyr::case_when(
        minor == 0 & Major == 1 ~ "del",
        minor == 0 & Major == 2 ~ "del",
        minor == 0 & (Major > 2 | Major == 0) ~ "del",
        minor > 0 & ((Major + minor) >= 3) ~ "amp", 
        .default = "none"
      )) 
    amp_tab <- filter(temp, type == "amp")
    del_tab <- filter(temp, type == "del")
    
    
    
    user_amp_gr <- GRanges(seqnames = amp_tab$chr,
                           ranges = IRanges(start = amp_tab$from, end = amp_tab$to))
    
    user_del_gr <- GRanges(seqnames = del_tab$chr,
                           ranges = IRanges(start = del_tab$from, end = del_tab$to))
    
    
    amp_overlaps <- findOverlaps(query = user_amp_gr,
                                 subject = amp_gr,
                                 type = "any", 
                                 select = "all")
    
    del_overlaps <- findOverlaps(query = user_del_gr,
                                 subject = del_gr,
                                 type = "any", 
                                 select = "all")
    
    res <- data.frame(
      cnv_idx = c(amp_overlaps@from, del_overlaps@from),
      hit_idx = c(amp_overlaps@to, del_overlaps@to),
      type = c(rep("amp", length(amp_overlaps@from)),
               rep("del", length(del_overlaps@from))),
      sample = sample_temp)  
    
    return(res)
  }) %>% 
    do.call(what = rbind)
  
  
  
  res_recurrent_cnv <- distinct(hit_cnv, hit_idx, type, sample) %>% 
    mutate(cnv_id = paste0(type, "_", hit_idx),
           case = case_name)  %>% 
    filter(length(cnv_id) > 0.2* length(unique(hit_cnv$sample)),
           .by = c(cnv_id)) 
  
  return(res_recurrent_cnv)
}



data_case1 <- list.files("~/project/crc/wgs_process/cnv_calling/case1/", full.names = T, pattern = "rds") %>%
	read_rds() %>% bind_rows()
data_case2 <- list.files("~/project/crc/wgs_process/cnv_calling/case2/", full.names = T, pattern = "rds") %>%
	read_rds() %>% bind_rows()

res_recurrent_cnv_case1 <- get_hit_seg(data_case1, 
                                       amp_gr = amp_gr,
                                       del_gr = del_gr,
                                       case_name = "case1" )

res_recurrent_cnv_case2 <- get_hit_seg(data_case2, 
                                       amp_gr = amp_gr,
                                       del_gr = del_gr,
                                       case_name = "case2" )



## cbind reults ----
res_recurrent_cnv <- rbind(res_recurrent_cnv_case1,
                           res_recurrent_cnv_case2) %>% 
  filter(length(unique(case)) == 2, .by = c(cnv_id)) %>%
  left_join(peak_tab, by = c("cnv_id", "type")) %>% 
  filter(p_value < 0.01) %>% 
  mutate(sample = factor(sample, level = c("PI_0", "PI_1", "PI_2", "PO_12", "PO_14", "PO_10", "PO_b",
                                           "LM1_1", "LM2_1", "LM3_1", 
                                           "LM4_0", "LM4_1", "LM4_10", "LM4_3", "LM4_4",
                                           "PLU_1", "PLU_2", "PLD_1", "PLD_4",
                                           "PRU_2", "PRU_3", "PRD_1", "PRD_2",
                                           "PLM1_1")))

x_intercepts <- seq(1.5, 24 -0.5, by = 1)
y_intercepts <- seq(1.5, length(unique(res_recurrent_cnv$cnv_id))-0.5, by = 1)

plot_cnv <- ggplot(res_recurrent_cnv, aes(x = sample, y = name)) +
  geom_point(aes(color = type), shape = 15, size = 6) +
  geom_vline(xintercept = x_intercepts, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = y_intercepts, color = "grey", linetype = "dashed") +
  scale_color_manual(values = c(del = "#77b3d4", amp = "#d45e00" )) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(angle = 45))

plot_cnv

