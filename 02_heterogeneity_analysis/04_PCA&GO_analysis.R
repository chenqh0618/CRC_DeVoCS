pjdir <- "~/project/crc/combined/"
library(tidyverse)
library(data.table)
library(parallel)
library(vcfR)
library(clusterProfiler)
library(ggrepel)
library(ggExtra)

#case1 -----------------------------
## data vaf 
vcf_file <- "/public/data/cqh_project/crc/wgs_process/hc_somatic_merge/annotated_CASE1.vcf"
data_vcf <- read.vcfR(vcf_file)

data_vaf <- extract_gt_tidy(data_vcf, format_fields = c("FREQ"))
data_vaf <- mutate(data_vaf, 
                   vaf = str_remove(gt_FREQ, "%") %>% 
                     as.numeric() %>% 
                     {./100} ,
                   vaf = if_else(is.na(vaf) , 0, vaf)) %>% 
  select(snv_id = Key, vaf, sample = Indiv) %>% 
  pivot_wider(names_from = snv_id, values_from = vaf, values_fill = 0) %>% 
  column_to_rownames("sample")
dim(data_vaf)

## SNV annoatation
data_site <- data_vcf@fix %>% as.data.frame() %>%
  separate(col = INFO, 
           sep = ";", 
           into = c("ac", "an", "dp", "gpv", "sf",
                    "snv_type", "spv", "ss", "ssc",
                    "ann", "lof", "nmd")) %>% 
  mutate(across(everything(.), function(x){str_remove(x, ".*=")}),
         chrom = getCHROM(data_vcf),
         pos = getPOS(data_vcf),
         snv_id = 1:nrow(.))



data_anno <- mclapply(1:nrow(data_site), mc.cores = 24, function(i){
  ann <- data_site$ann[i]
  temp_ann <- ann %>% 
    str_split(pattern = ",", simplify = T) %>% 
    unlist() %>% 
    str_split("\\|", simplify = T) %>% 
    as.data.frame() %>% 
    mutate(snv_id = data_site$snv_id[i],
           chrom = data_site$chrom[i],
           pos = data_site$pos[i])
  return(temp_ann)
}) %>% do.call(what = bind_rows) %>% 
  `colnames<-`(c("allele", "annotation", "impact",
                 "gene", "gene_id", 
                 "feature", "feature_id", "transcript_biotype", "rank",
                 "hgvs_c", "hgvs_p", 
                 "cdnaPos_2_cdnaLength", "cdsPos_2_cdsLength", "aaPos_2_aaLength",
                 "distance", "info",
                 "snv_id", "chrom", "pos"))
count(data_anno, impact)

## pca analysis
res_pca <- prcomp(x = data_vaf)
head(res_pca$x)
summary(res_pca)

temp <- res_pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(region = if_else(grepl("LM", sample), "Metastatsis", "Primary"))



ggplot(temp, aes(x = PC1, y =PC2)) +
  geom_point(aes(color = region), size = 3) +
  geom_text_repel(aes(label=sample), size  = 5) +
  labs(x = "PC1 (40.95%)",
       y = "PC2 (10.37%)") +
  theme(axis.title = element_text(size = 12))

## outlier genes in PC2
idx <- res_pca$rotation[,2]
length(idx)
idx_outlider <- idx[which(idx > mean(idx)+3*sd(idx) | idx < mean(idx)-3*sd(idx))]
length(idx_outlider)
outlier_snv <- names(idx_outlider)

temp <- filter(data_anno,
               snv_id %in% outlier_snv,
               impact %in% c("HIGH", "MODERATE", "LOW"))
outlier_gene1 <- unique(temp$gene)



# case2 ---------------------------
## data vaf
vcf_file <- "/public/data/cqh_project/crc/wgs_process/hc_somatic_merge/annotated_CASE2.vcf"

data_vaf <- extract_gt_tidy(vcf_file, format_fields = c("FREQ"))
data_vaf <- mutate(data_vaf, 
                   vaf = str_remove(gt_FREQ, "%") %>% 
                     as.numeric() %>% 
                     {./100} ,
                   vaf = if_else(is.na(vaf) , 0, vaf)) %>% 
  select(snv_id = Key, vaf, sample = Indiv) %>% 
  pivot_wider(names_from = snv_id, values_from = vaf, values_fill = 0) %>% 
  column_to_rownames("sample")
dim(data_vaf)

## data annotation
data_site <- data_vcf@fix %>% as.data.frame() %>%
  separate(col = INFO, 
           sep = ";", 
           into = c("ac", "an", "dp", "gpv", "sf",
                    "snv_type", "spv", "ss", "ssc",
                    "ann", "lof", "nmd")) %>% 
  mutate(across(everything(.), function(x){str_remove(x, ".*=")}),
         chrom = getCHROM(data_vcf),
         pos = getPOS(data_vcf),
         snv_id = 1:nrow(.))



data_anno <- mclapply(1:nrow(data_site), mc.cores = 24, function(i){
  ann <- data_site$ann[i]
  temp_ann <- ann %>% 
    str_split(pattern = ",", simplify = T) %>% 
    unlist() %>% 
    str_split("\\|", simplify = T) %>% 
    as.data.frame() %>% 
    mutate(snv_id = data_site$snv_id[i],
           chrom = data_site$chrom[i],
           pos = data_site$pos[i])
  return(temp_ann)
}) %>% do.call(what = bind_rows) %>% 
  `colnames<-`(c("allele", "annotation", "impact",
                 "gene", "gene_id", 
                 "feature", "feature_id", "transcript_biotype", "rank",
                 "hgvs_c", "hgvs_p", 
                 "cdnaPos_2_cdnaLength", "cdsPos_2_cdsLength", "aaPos_2_aaLength",
                 "distance", "info",
                 "snv_id", "chrom", "pos"))



## pca analysis
res_pca <- prcomp(x = data_vaf)
head(res_pca$x)
summary(res_pca)

temp <- res_pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(region = if_else(grepl("LM", sample), "Metastatsis", "Primary"))

ggplot(temp, aes(x = PC1, y =PC2)) +
  geom_point(aes(color = region), size = 3) +
  geom_text_repel(aes(label=sample), size  = 5) +
  labs(x = "PC1 (52.47%)",
       y = "PC2 (16.52%)") +
  theme(axis.title = element_text(size = 12))

## outlier genes in PC2
idx <- res_pca$rotation[,2]
length(idx)
idx_outlider <- idx[which(idx > mean(idx)+3*sd(idx) | idx < mean(idx)-3*sd(idx))]
length(idx_outlider)
outlier_snv <- names(idx_outlider)

temp <- filter(data_anno,
               snv_id %in% outlier_snv,
               impact %in% c("HIGH", "MODERATE", "LOW"))
outlier_gene2 <- unique(temp$gene)





# intersection between two cases
intersect_gene_list <- intersect(outlier_gene1, outlier_gene2)

res_enrich <- enrichGO(gene = intersect_gene_list,
                       OrgDb = "org.Hs.eg.db",
                       ont = "all",
                       pAdjustMethod = "BH",
                       minGSSize = 10,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       keyType="SYMBOL")

ggplot(res_enrich,
       aes(x = GeneRatio, y = Description, color = -log10(pvalue)))+
  geom_point(aes(size =Count))+
  theme_bw()+
  scale_y_discrete(labels = function(y) str_wrap(y,width=50))+
  labs(size = "Counts", x = "GeneRatio", y = "GO terms", title = "GO Enrichment") +
  scale_color_gradient(low="blue",high ="red")+
  theme(axis.text = element_text(size=18,color="black"),
        axis.text.x = element_text(angle = 90),
        axis.title=element_text(size=16),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title=element_text(size=16),
        legend.text = element_text(size = 16),
        title=element_text(size=20),
        strip.text = element_text(size = 14))+
  facet_grid(ONTOLOGY~., scales ="free",space="free")+
  guides(color = guide_colorbar (reverse = TRUE))



gene_list <- res_enrich@result %>% 
  filter(ID %in% c("GO:0007156", "GO:0098742")) %>% 
  pull(geneID) %>% 
  str_split(pattern = "/") %>% 
  unlist %>% 
  unique()

driver_gene_list <- readLines("/public/public_data/biotrainee/tumor_driver_gene/Alex_nature_2024_CRC_driver_gene_list.txt")
intersect(gene_list, driver_gene_list)

# "PCDH15"  "DSCAML1" "KIRREL3" "IGSF21"  "CNTN6"   "CNTN4"   "ROBO2"   "ROBO1"   "CADM2"  
# "TENM3"   "CDH18"   "CDH12"   "CDH9"    "PCDHGA3" "PCDHGA2" "PCDHGA1" "SDK1"   
# "DAB1"    "TGFB2"   "NRXN1"   "GRID2"   "PTPRD"  
