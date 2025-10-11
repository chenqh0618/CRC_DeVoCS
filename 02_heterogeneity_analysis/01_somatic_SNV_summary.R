
library(tidyverse)
library(data.table)
library(webr)
library(vcfR)
library(ggpattern)
maf_thred <- 0.05


# case1 ----
## input data ----
vcf_file <- "/public/data/cqh_project/crc/wgs_process/hc_somatic_merge/merged_hc_somatic_CASE1.vcf"
vcf_data <- read.vcfR(vcf_file, verbose = F)
colnames(vcf_data@gt) ; dim(vcf_data@gt)

site <- extract.gt(vcf_data, "AD") %>% rownames()
ad <- extract.gt(vcf_data, "AD") %>% 
  apply(2, as.numeric) %>% 
  as.data.frame() %>% 
  mutate(site = site) %>% 
  pivot_longer(cols = -site, names_to = "sample", values_to = "ad")
rd <- extract.gt(vcf_data, "RD") %>% 
  apply(2, as.numeric) %>% 
  as.data.frame() %>% 
  mutate(site = site) %>% 
  pivot_longer(cols = -site, names_to = "sample", values_to = "rd")

data_raw <- left_join(ad, rd, by = c("site", "sample")) %>% 
  mutate(ad = replace_na(ad, 0),
         rd = replace_na(rd, 50),
         dp = ad + rd,
         vaf = ad/dp,
         vaf = if_else(is.na(vaf), 0, vaf)) %>% 
  mutate(region = case_when(sample %in% c("PI_0", "PI_1", "PI_2") ~ "PI",
                            sample %in% c("PO_12", "PO_14", "PO_10", "PO_b") ~ "PO",
                            sample == "LM1_1" ~ "LM1",
                            sample == "LM2_1" ~ "LM2",
                            sample == "LM3_1" ~ "LM3",
                            grepl("LM4", sample) ~ "LM4",
                            .default = sample),
         tumor = if_else(grepl("LM", sample), "metastasis", "primary"),
         gt = if_else(vaf >= maf_thred , 1, 0)) %>% 
  filter(!sample %in% c("NORMAL")) %>% 
  rename(snv_id = site)
				   
										   
										   
										   
snv_type_tab <- group_by(data_raw, snv_id, tumor) %>% 
  reframe(gt = max(gt)) %>% 
  pivot_wider(names_from = tumor, values_from = gt) %>% 
  mutate(idx = case_when(primary > 0 & metastasis > 0 ~ "shared",
                         primary > 0 & metastasis == 0 ~ "primary",
                         primary == 0 & metastasis > 0 ~ "metastasis")) %>% 
  select(snv_id, snv_type = idx) 
data_raw <- left_join(data_raw, snv_type_tab)


## SNV groups --------------
res_region <- lapply(unique(data_raw$snv_type), function(idx){
  temp <- filter(data_raw, snv_type == idx, gt > 0) %>% 
    group_by(snv_id, region) %>% 
    reframe(gt_region = max(gt)) 
  max_region_num = length(unique(temp$region))
  
  res <- group_by(temp, snv_id) %>% 
    reframe(gt_region_num = sum(gt_region)) %>% 
    mutate(region_idx = case_when(gt_region_num == max_region_num ~ "Fixed",
                                  gt_region_num == 1 ~ "Specific",
                                  .default = "Partial")) %>% 
    select(snv_id, region_idx) %>% 
    mutate(snv_type = idx)
  return(res)
}) %>% 
  do.call(what = rbind)

res_sample <- lapply(unique(data_raw$snv_type), function(idx){
  temp <- filter(data_raw, snv_type == idx, gt > 0) %>% 
    group_by(snv_id, sample) %>% 
    reframe(gt_sample = max(gt)) 
  max_sample_num = length(unique(temp$sample))
  
  res <- group_by(temp, snv_id) %>% 
    reframe(gt_sample_num = sum(gt_sample)) %>% 
    mutate(sample_idx = case_when(gt_sample_num == max_sample_num ~ "Fixed",
                                  gt_sample_num == 1 ~ "Specific",
                                  .default = "Partial")) %>% 
    select(snv_id, sample_idx) %>% 
    mutate(snv_type = idx)
  return(res)
}) %>% 
  do.call(what = rbind)
res <- left_join(res_region, res_sample, by = c("snv_id", "snv_type"))

case1_snv_pieplot = PieDonut(res, 
                             mapping =aes(x = snv_type,  y = region_idx),
                             addPieLabel = F,
                             addDonutLabel = F,
                             showRatioDonut = TRUE,
                             showRatioPie = F,
                             ratioByGroup = F,
                             showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.0001),
                             labelposition = getOption("PieDonut.labelposition", 1),
                             labelpositionThreshold = 0.001,
                             # explode = 1,
                             # explodePie = TRUE,
                             # explodeDonut=F,
                             showPieName = F,
                             showDonutName = FALSE,
                             title  = "Region_based summary",
                             titlesize = 7,
                             r0=0.6,
                             start=3*pi/2, 
                             pieAlpha = 0.8,
                             donutAlpha = 1)  



## upset plot ---------------
library(UpSetR)
temp <- select(data_raw, snv_id, gt, region) %>% 
  group_by(region, snv_id) %>% 
  reframe(gt = max(gt)) %>% 
  pivot_wider(names_from = region, 
              values_from = gt) %>% 
  column_to_rownames("snv_id")

list_intersection = list("PI", "PO", "LM1", "LM2", "LM3", "LM4",
                         c("PI", "PO"),
                         list("LM1", "LM2", "LM3", "LM4"),
                         list("LM1", "LM2", "LM3", "LM4", "PI"),
                         list("LM1", "LM2", "LM3", "LM4", "PO"),
                         list("LM1",  "PI","PO" ),
                         list("LM2",  "PI","PO" ),
                         list("LM3",  "PI","PO" ),
                         list("LM4",  "PI","PO" ),
                         list("LM1", "LM2", "LM3", "LM4", "PI","PO"))

case1_upset_plot1 <- upset(temp, 
                           nsets = 7, 
                           # nintersects = NA,
                           intersections = list_intersection,
                           keep.order = T,
                           sets = c("PI", "PO", "LM1", "LM2", "LM3", "LM4"),
                           group.by = "degree",
                           order.by = c("freq", "degree"),
                           # scale.sets = "log10",
                           # scale.intersections = "log10",
                           point.size = 4,
                           text.scale = 2
)


temp <- select(data_raw, sample, snv_id, gt) %>% 
  filter(!is.na(gt)) %>% 
  pivot_wider(names_from = sample, values_from = gt) %>% 
  column_to_rownames("snv_id")

upset(temp, nsets = 15,  
      sets = colnames(temp),
      keep.order = T,
      order.by = c("freq"),
      decreasing = c(T))  




## sample SNV component ----

case1_sample_snv_type_tab <- filter(data_raw, vaf > maf_thred) %>% 
  count(sample, snv_type) %>% 
  mutate(snv_type = factor(snv_type, levels = c("shared", "primary", "metastasis")),
         sample = factor(sample, level = c("PI_0", "PI_1", "PI_2", "PO_12", "PO_14", "PO_10", "PO_b",
                                           "LM1_1", "LM2_1", "LM3_1", 
                                           "LM4_0", "LM4_1", "LM4_10", "LM4_3", "LM4_4",
                                           "PLU_1", "PLU_2", "PLD_1", "PLD_4",
                                           "PRU_2", "PRU_3", "PRD_1", "PRD_2",
                                           "PLM1_1")))

case1_sample_private_snv_tab <- filter(data_raw, vaf > maf_thred) %>% 
  filter(length(sample) == 1, .by = c(snv_id)) %>% 
  count(sample) 


case1_snv_component_plot <- ggplot(case1_sample_snv_type_tab) +
  geom_bar(aes(x = sample, y = n, fill = snv_type),
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = tumor_color_list) +
  geom_bar_pattern(data = case1_sample_private_snv_tab,
                   mapping = aes(x = sample, y = n),
                   alpha = 0.1,
                   pattern = "circle",
                   fill = "white",
                   pattern_spacing = 0.02,
                   stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "bottom")


# case2 ----
## data input ----
vcf_file <- "/public/data/cqh_project/crc/wgs_process/hc_somatic_merge/merged_hc_somatic_CASE2.vcf"
vcf_data <- read.vcfR(vcf_file, verbose = F)
colnames(vcf_data@gt) ; dim(vcf_data@gt)

site <- extract.gt(vcf_data, "AD") %>% rownames()
ad <- extract.gt(vcf_data, "AD") %>% 
  apply(2, as.numeric) %>% 
  as.data.frame() %>% 
  mutate(site = site) %>% 
  pivot_longer(cols = -site, names_to = "sample", values_to = "ad")
rd <- extract.gt(vcf_data, "RD") %>% 
  apply(2, as.numeric) %>% 
  as.data.frame() %>% 
  mutate(site = site) %>% 
  pivot_longer(cols = -site, names_to = "sample", values_to = "rd")

data_raw <- left_join(ad, rd, by = c("site", "sample")) %>% 
  mutate(ad = replace_na(ad, 0),
         rd = replace_na(rd, 50),
         dp = ad + rd,
         vaf = ad/dp,
         vaf = if_else(is.na(vaf), 0, vaf)) %>% 
  mutate(region = str_sub(sample, 1, 3),
         tumor = if_else(grepl("LM", sample), "metastasis", "primary"),
         gt = if_else(vaf >= maf_thred , 1, 0),
         sample = factor(sample, 
                         levels = c("PLU_1", "PLU_2", "PLD_1", "PLD_4",
									"PRU_2", "PRU_3", "PRD_1", "PRD_2",
                                    "LM1_1"))) %>% 
  filter(!sample %in% c("NORMAL")) %>% 
  rename(snv_id = site)

snv_type_tab <- group_by(data_raw, snv_id, tumor) %>% 
  reframe(gt = max(gt)) %>% 
  pivot_wider(names_from = tumor, values_from = gt) %>% 
  mutate(idx = case_when(primary > 0 & metastasis > 0 ~ "shared",
                         primary > 0 & metastasis == 0 ~ "primary",
                         primary == 0 & metastasis > 0 ~ "metastasis")) %>% 
  select(snv_id, snv_type = idx) 
data_raw <- left_join(data_raw, snv_type_tab)




## SNV groups --------------
res_region <- lapply(unique(data_raw$snv_type), function(idx){
  temp <- filter(data_raw, snv_type == idx, gt > 0) %>% 
    group_by(snv_id, region) %>% 
    reframe(gt_region = max(gt)) 
  max_region_num = length(unique(temp$region))
  
  res <- group_by(temp, snv_id) %>% 
    reframe(gt_region_num = sum(gt_region)) %>% 
    mutate(region_idx = case_when(gt_region_num == max_region_num ~ "Fixed",
                                  gt_region_num == 1 ~ "Specific",
                                  .default = "Partial")) %>% 
    select(snv_id, region_idx) %>% 
    mutate(snv_type = idx)
  return(res)
}) %>% 
  do.call(what = rbind)

res_sample <- lapply(unique(data_raw$snv_type), function(idx){
  temp <- filter(data_raw, snv_type == idx, gt > 0) %>% 
    group_by(snv_id, sample) %>% 
    reframe(gt_sample = max(gt)) 
  max_sample_num = length(unique(temp$sample))
  
  res <- group_by(temp, snv_id) %>% 
    reframe(gt_sample_num = sum(gt_sample)) %>% 
    mutate(sample_idx = case_when(gt_sample_num == max_sample_num ~ "Fixed",
                                  gt_sample_num == 1 ~ "Specific",
                                  .default = "Partial")) %>% 
    select(snv_id, sample_idx) %>% 
    mutate(snv_type = idx)
  return(res)
}) %>% 
  do.call(what = rbind)
res <- left_join(res_region, res_sample, by = c("snv_id", "snv_type"))

PieDonut(res, 
         mapping =aes(x = snv_type,  y = region_idx),
         addPieLabel = F,
         addDonutLabel = F,
         showRatioDonut = TRUE,
         showRatioPie = F,
         ratioByGroup = F,
         showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.001),
         labelposition = getOption("PieDonut.labelposition", 1),
         labelpositionThreshold = 0.001,
         # explode = 1,
         # explodePie = TRUE,
         # explodeDonut=F,
         showPieName = F,
         showDonutName = FALSE,
         title  = "Region_based summary",
         titlesize = 7,
         r0=0.6,
         start=3*pi/2, 
         pieAlpha = 0.8,
         donutAlpha = 1,)  



## upset plot ---------------
temp <- select(data_raw, snv_id, gt, region) %>% 
  group_by(region, snv_id) %>% 
  reframe(gt = max(gt)) %>% 
  pivot_wider(names_from = region, 
              values_from = gt) %>% 
  column_to_rownames("snv_id")

list_intersection = list("PLU", "PLD", "PRU", "PRD", "LM1") %>% 
  c(lapply(c("PLU", "PLD", "PRD", "PRU"), function(x) {c(x, "LM1")}),
    list(c("PLD", "PRD", "PRU"),
         c("PLU", "PRD", "PRU"),
         c("PLU", "PLD", "PRU"),
         c("PLU", "PLD", "PRD"),
         c("PLU", "PLD", "PRD", "PRU"),
         c("PLU", "PLD", "PRD", "PRU", "LM1")))

case2_upset_plot <- upset(temp, 
                          nsets = 7, 
                          order.by = c("freq",  "degree"),
                          intersections = list_intersection,
                          keep.order = T,
                          sets = c("PLU", "PLD", "PRU", "PRD", "LM1"),
                          group.by = "degree",
                          point.size = 4,
                          text.scale = 2)



temp <- select(data_raw, sample, snv_id, gt) %>% 
  filter(!is.na(gt)) %>% 
  pivot_wider(names_from = sample, values_from = gt) %>% 
  column_to_rownames("snv_id")

upset(temp, nsets = 15,  
      sets = colnames(temp),
      keep.order = T,
      order.by = c("freq"),
      decreasing = c(T))




## sample SNV component ----
case2_sample_snv_type_tab <- filter(data_raw, vaf > maf_thred) %>% 
  count(sample, snv_type) %>% 
  mutate(snv_type = factor(snv_type, levels = c("shared", "primary", "metastasis")),
         sample = factor(sample, level = c("PI_0", "PI_1", "PI_2", "PO_12", "PO_14", "PO_10", "PO_b",
                                           "LM1_1", "LM2_1", "LM3_1", 
                                           "LM4_0", "LM4_1", "LM4_10", "LM4_3", "LM4_4",
                                           "PLU_1", "PLU_2", "PLD_1", "PLD_4",
                                           "PRU_2", "PRU_3", "PRD_1", "PRD_2",
                                           "PLM1_1"))) 

case2_sample_private_snv_tab <- filter(data_raw, vaf > maf_thred) %>% 
  filter(length(sample) == 1, .by = c(snv_id)) %>% 
  count(sample)


case2_snv_component_plot <- ggplot(case2_sample_snv_type_tab) +
  geom_bar(aes(x = sample, y = n, fill = snv_type),
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = tumor_color_list) +
  geom_bar_pattern(data = case2_sample_private_snv_tab,
                   mapping = aes(x = sample, y = n),
                   alpha = 0.1,
                   pattern = "circle",
                   fill = "white",
                   pattern_spacing = 0.02,
                   stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "bottom")



# combine 

case1_upset_plot 
case2_upset_plot 


case1_snv_component_plot + 
  case2_snv_component_plot +
  plot_layout(ncol =  2, widths = c(15,9))


sample_snv_type_tab <- rbind(case1_sample_snv_type_tab, case2_sample_snv_type_tab)
sample_private_snv_tab <- rbind(case1_sample_private_snv_tab, case2_sample_private_snv_tab)  

ggplot(sample_snv_type_tab) +
  geom_bar(aes(x = sample, y = n, fill = snv_type),
           width = 0.6, 
           stat = "identity",
           position = "stack") +
  scale_fill_manual(values = tumor_color_list) +
  geom_bar_pattern(data = sample_private_snv_tab,
                   mapping = aes(x = sample, y = n),
                   alpha = 0.1,
                   width = 0.6,
                   pattern = "circle",
                   fill = "white",
                   pattern_spacing = 0.02,
                   stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45),
        axis.text = element_text(size = 14),
        legend.position = "bottom")
