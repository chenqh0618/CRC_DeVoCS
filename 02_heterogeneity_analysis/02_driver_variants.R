#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(parallel)
library(clusterProfiler)
library(maftools)
library(fanyi)
library(patchwork)

source("~/software/script_tools/R_color.R")

# SNV ----
vc_cols = c(Missense_Mutation = "#76afda",
            Nonsense_Mutation = "#7CC767",
            Translation_Start_Site = "#BD6263",
            Multi_Hit = "black")
fabcolors <- list(sample_type = c(primary = "darkgreen", metastasis = "darkred"))

driver_gene_list <- readLines("/public/public_data/biotrainee/tumor_driver_gene/Alex_nature_2024_CRC_driver_gene_list.txt")

# case1 data input ----
data <- fread("/public/data/cqh_project/crc/wgs_process/hc_somatic_merge/annovar_merged_hc_somatic_CASE1.hg19_multianno.txt")

data$ExonicFunc.refGene %>% unique()
data <- filter(data,
               Otherinfo10 == "PASS",
               ExonicFunc.refGene %in% c("nonsynonymous SNV", 
										 "frameshif",
                                         "stopgain", 
                                         "stoploss", 
                                         "startloss"),
)  %>% 
  mutate(label = paste0(Chr, "__", Start))

data <- lapply(1:nrow(data), function(i){
  x <- str_split(data$Gene.refGene[i], ";") %>% unlist()
  if(any(x %in% driver_gene_list)){return(data[i,])}
}) %>% 
  do.call(what = rbind)


case1_gene_tab <- dplyr::select(data, 
                                chrom = Chr, pos = Start, 
                                snv_id = label,
                                gene = Gene.refGene,
                                gene_func = Func.refGene,
                                snv_func = ExonicFunc.refGene
) %>% 
  mutate(case = "case1")

vaf_case1 <- read_rds("snv_calling/case1_rescued_snv.rds") %>% 
  mutate(snv_id = paste0(chrom, "__", pos)) %>% 
  filter(snv_id %in% case1_gene_tab$snv_id) %>% 
  left_join(case1_gene_tab, by = c("snv_id", "chrom", "pos"))




# case2 data input ----
data <- fread("/public/data/cqh_project/crc/wgs_process/hc_somatic_merge/annovar_merged_hc_somatic_CASE2.hg19_multianno.txt")
head(data); dim(data)
data <- filter(data, 
               Otherinfo10 == "PASS",
               ExonicFunc.refGene %in% c("nonsynonymous SNV", 
										 "frameshif",
                                         "stopgain",
                                         "stoploss", 
                                         "startloss")
) %>% 
  mutate(label = paste0(Chr, "__", Start))

data <- lapply(1:nrow(data), function(i){
  x <- str_split(data$Gene.refGene[i], ";") %>% unlist()
  if(any(x %in% driver_gene_list)){return(data[i,] )}
}) %>% 
  do.call(what = rbind)

case2_gene_tab <- dplyr::select(data,
                                chrom = Chr, pos = Start, 
                                snv_id = label,
                                gene = Gene.refGene,
                                gene_func = Func.refGene,
                                snv_func = ExonicFunc.refGene
) %>% 
  mutate(case = "case2",
         chrom = as.character(chrom))

vaf_case2 <- read_rds("snv_calling/case2_rescued_snv.rds") %>% 
  mutate(snv_id = paste0(chrom, "__", pos))%>% 
  filter(snv_id %in% case2_gene_tab$snv_id) %>% 
  left_join(case2_gene_tab, by = c("snv_id", "chrom", "pos"))


# combine two samples
data <- rbind(vaf_case1, vaf_case2) %>% 
  filter(sum(vaf > 0.05) > 2, .by = c(case, snv_id)) %>% 
  mutate(sample = factor(sample, level = c("PI_0", "PI_1", "PI_2", "PO_12", "PO_14", "PO_10", "PO_b",
                                           "LM1_1", "LM2_1", "LM3_1", 
                                           "LM4_0", "LM4_1", "LM4_10", "LM4_3", "LM4_4",
                                           "PLU_1", "PLU_2", "PLD_1", "PLD_4",
                                           "PRU_2", "PRU_3", "PRD_1", "PRD_2",
                                           "PLM1_1")))  

x_intercepts <- seq(1.5, 24 -0.5, by = 1)
y_intercepts <- seq(0.5, 21 -0.5, by = 1)
y_intercepts <- y_intercepts[!(y_intercepts %in% c(2.5, 7.5))]

data2 <- filter(data,  vaf > 0.5)


plot_snv <- ggplot(data = data2 , aes(x = sample, y = paste0(gene, "_", snv_id), color = snv_func)) +
  geom_point(aes(size = vaf )) +
  # geom_point(aes( alpha = vaf, size = (vaf > 0.05)*4),   shape = 15) +
  scale_color_manual(values = c("#66c2a4", "#e789c3" )) +
  geom_vline(xintercept = x_intercepts, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = y_intercepts, color = "grey", linetype = "dashed") +
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0, "in"))
plot_snv



plot_snv/plot_cnv +
  plot_layout(ncol = 1, heights = c(5,3))

