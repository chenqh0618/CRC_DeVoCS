library(tidyverse)
library(parallel)
library(data.table)
library(ggpubr)

get_fst <- function(data){
  require(parallel)
  require(tidyverse)
  ## 3 columns needed: snv_id, sample, vaf
  ## Nei's FST
  calc_fst <- function(p1, p2) {
    q1 <- 1 - p1
    q2 <- 1 - p2
    p_bar <- (p1 + p2) / 2
    q_bar <- 1 - p_bar
    hs <- (2 * p1 * q1 + 2 * p2 * q2) / 2
    ht <- 2 * p_bar * q_bar
    if(ht == 0) return(NA)  # avoid 0
    fst <- (ht - hs) / ht
    return(fst)
  }
  
  ## sample pairwise 
  pairs <- t(combn(unique(data$sample), 2))
  
  ## data impute
  vaf_matrix <- dplyr::select(data, snv_id, sample, vaf) %>% 
    pivot_wider(names_from = sample, values_from = vaf)
  
  ## calculation
  res_fst <- mclapply(1:nrow(pairs), mc.cores = 12,function(i){
    pop1 <- pairs[i,1]
    pop2 <- pairs[i,2]
    vaf1 <- vaf_matrix[, pop1] %>% t
    vaf2 <- vaf_matrix[, pop2] %>% t
    fst_per_locus <- mapply(calc_fst, vaf1, vaf2)
    fst_overall <- mean(fst_per_locus, na.rm = TRUE)
    res <- data.frame(pop1 = pop1,
                      pop2 = pop2,
                      fst = fst_overall)
  }) %>% 
    do.call(what = rbind)
  
  ## result
  return(res_fst)
}


# WGS case1 ----
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
				   
		

res <- get_fst(data_raw)


dist_mat <- dplyr::mutate(res, pop1_new = pop2, pop2_new = pop1) %>% 
  dplyr::select(pop1 = pop1_new, pop2 = pop2_new, fst) %>% 
  rbind(res) %>% 
  pivot_wider(names_from = pop1, 
              values_from = fst) %>% 
  column_to_rownames("pop2") 
dist_mat <- dist_mat[colnames(dist_mat), colnames(dist_mat)]
dist_mat <- as.dist(dist_mat, upper = T, diag = T) 


tree_res <- hclust(dist_mat, method = "complete")  
plot(tree_res)

sample_idx <- tree_res$labels[tree_res$order]
dist_mat <- as.matrix(dist_mat)[sample_idx, sample_idx]
diag(dist_mat) <- NA

pheatmap::pheatmap(dist_mat,
                   cluster_rows = F, 
                   cluster_cols = F, 
                   color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                             "RdBu")))(100),
                   border_color = "white", angle_col = 45)



# WGS case2 ----
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

res <- get_fst(data_raw)

dist_mat <- dplyr::mutate(res, pop1_new = pop2, pop2_new = pop1) %>% 
  dplyr::select(pop1 = pop1_new, pop2 = pop2_new, fst) %>% 
  rbind(res) %>% 
  pivot_wider(names_from = pop1, 
              values_from = fst) %>% 
  column_to_rownames("pop2") 
dist_mat <- dist_mat[colnames(dist_mat), colnames(dist_mat)]
dist_mat <- as.dist(dist_mat, upper = T, diag = T) 

tree_res <- hclust(dist_mat, method = "complete")  
plot(tree_res)

sample_idx <- tree_res$labels[tree_res$order]
dist_mat <- as.matrix(dist_mat)[sample_idx, sample_idx]
diag(dist_mat) <- NA

pheatmap::pheatmap(dist_mat,
                   cluster_rows = F, 
                   cluster_cols = F, 
                   color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                             "RdBu")))(100),
                   border_color = "white", angle_col = 45)

# AMP case1 ----
data_raw <- read_rds("/public/data/cqh_project/crc/amp_data/filtered_CASE1.rds")
res_fst <- dplyr::select(data_raw,
                     snv_id = site, 
                     vaf = maf, 
                     sample) %>% 
  get_fst()


temp <- dplyr::select(data_raw, sample, region) %>% 
  distinct()
res_fst <- left_join(res_fst, temp, by = c("pop1" = "sample")) %>% 
  rename(region1 = region) %>% 
  left_join(temp, by = c("pop2" = "sample")) %>% 
  rename(region2 = region)

## inter tumor ---- 

inter_tumor_fst <- filter(res_fst, region1 != region2) %>% 
  mutate(tumor1 = if_else(region1 %in% c("PO", "PI"), 0, 1),
         tumor2 = if_else(region2 %in% c("PO", "PI"), 0, 1),
         type = tumor1 + tumor2,
         type = case_when(type == 0 ~ "Primary",
                          type == 1 ~ "Primary-Metastasis",
                          type == 2 ~ "Metastasis") %>% 
           factor(levels = c("Primary",  "Primary-Metastasis", "Metastasis" )))

plot1 <- ggplot(inter_tumor_fst, aes(x = type, y = fst, fill = type)) +
  geom_violin() +
  geom_signif(comparisons = list(c("Primary","Primary-Metastasis"),
                                 c("Metastasis","Primary-Metastasis"),
                                 c("Primary","Metastasis")), 
              test = t.test,   
              y_position = c(0.09, 0.118, 0.125),
              map_signif_level = T,
              size= 0.5, textsize = 5,
              color="black") +
  scale_fill_manual(values = c("darkgreen", "steelblue", "darkred")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom")



## primary tumor ----
primary_fst <- filter(res_fst, 
               region1 %in% c("PO", "PI"),
               region2 %in% c("PO", "PI")) %>% 
  mutate(type = paste0(region1, "_to_", region2))

plot2 <- ggplot(primary_fst, aes(x = type, y = fst, fill = type)) +
  geom_violin() +
  geom_signif(comparisons = list(c("PI_to_PO","PI_to_PI"),
                                 c("PI_to_PO","PO_to_PO"),
                                 c("PO_to_PO","PI_to_PI")), 
              test = t.test, 
              map_signif_level = T,
              size= 0.5, textsize = 5,
              y_position = c(0.04,0.04,0.042),
              color="black") +
  theme_bw() +
  scale_fill_manual(values = c("#044827", "#009E73", "#9DCF97" )) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom") 
  



## metastasis tumor ----
metastasis_fst <- filter(res_fst, 
               !region1 %in% c("PO", "PI"),
               !region2 %in% c("PO", "PI")) %>% 
  mutate(type = if_else(region1 == region2, region1, "Intertumor"))


plot3 <- ggplot(metastasis_fst, aes(x = type, y = fst, fill = type)) +
  geom_violin()  + 
  geom_signif(comparisons = list(c("LM4","Intertumor"),
                                 c("LM1","LM2"),
                                 c("LM1","LM3"),
                                 c("LM1","LM4"),
                                 c("LM3","LM4")), 
              test = t.test,  ##计算方法
              y_position = c(0.12, 0.015, 0.03, 0.066, 0.06), 
              map_signif_level = T,
              size= 0.5, textsize = 5, 
              color="black") +
  theme_bw() +
  scale_fill_manual(values = c("grey", "#F7C2C1", "#FC7E47", "#E82EB8", "#991E23")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom") 



plot2 + plot1 + plot3 +
  plot_layout(ncol = 3, widths = c(3,3,5))

# AMP case2 ----
data_raw <- read_rds("public/data/cqh_project/crc/amp_data/filtered_CASE2.rds")

res_fst <- dplyr::select(data_raw,
                         snv_id = site, 
                         vaf = maf, 
                         sample) %>% 
  get_fst()


temp <- dplyr::select(data_raw, sample, region) %>% 
  distinct()
res_fst <- left_join(res_fst, temp, by = c("pop1" = "sample")) %>% 
  rename(region1 = region) %>% 
  left_join(temp, by = c("pop2" = "sample")) %>% 
  rename(region2 = region)


filter(res_fst, region1 == "LM", region2 == "LM" )

## inter tumor 
inter_tumor_fst <- filter(res_fst) %>% 
  mutate(tumor1 = if_else(region1 %in% c("LM", "normal"), 1, 0),
         tumor2 = if_else(region2 %in% c("LM", "normal"), 1, 0),
         type = tumor1 + tumor2,
         type = case_when(type == 0 ~ "Primary",
                          type == 1 ~ "Primary-Metastasis",
                          type == 2 ~ "Metastasis") %>% 
           factor(levels = c("Primary", "Primary-Metastasis" , "Metastasis")))

plot1 <- ggplot(inter_tumor_fst, aes(x = type, y = fst, fill = type)) +
  geom_violin() +
  geom_signif(comparisons = list(c("Primary","Primary-Metastasis"),
                                 c("Metastasis","Primary-Metastasis"),
                                 c("Primary","Metastasis")), 
              test = t.test,  ##计算方法
              y_position = c(0.18,  0.19, 0.20), 
              
              size=0.5,color="black") +
  scale_fill_manual(values = c("darkgreen", "steelblue", "darkred")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom")



## primary tumor
metastasis_fst <- filter(res_fst, 
                         !region1 %in% c("LM", "normal"),
                         !region2 %in% c("LM", "normal")) %>% 
  mutate(type = if_else(region1 == region2, region1, "Inter-region") %>% 
           factor(levels = c("Inter-region", "LU", "RU", "LD", "RD" )))


plot2 <- ggplot(metastasis_fst, aes(x = type, y = fst , fill = type)) +
  geom_violin()  + 
  geom_signif(comparisons = list(c("LU","Inter-region"),
                                 c("RU","Inter-region"),
                                 c("LD","Inter-region"),
                                 c("RD","Inter-region")), 
              test = t.test,  ##计算方法
              y_position = c(0.06, 0.068, 0.075, 0.082), 
              size=0.5,color="black") +
  scale_fill_manual(values = c("grey","#B7DBAD", "#91AB28", "#34B388", "#044827")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom")


plot2 + plot1 +
  plot_layout(ncol = 2, widths = c(3, 5))


