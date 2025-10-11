library(tidyverse)
library(data.table)
library(parallel)
library(vcfR)
library(treedataverse)

res_dir <- "/public/data/cqh_project/crc/sample_phylogeny/wgs_MLtree/"

## input data
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
         vaf = if_else(is.na(vaf), 0, vaf)) 
				   
				
seq <- data.frame(site = snv_id,
                  ref = getREF(vcf_data),
                  alt = getALT(vcf_data))
res <- left_join(data_raw, seq, by = "site")
head(res)

## get sequences
vaf_thd = 0.05
res_seq_tumor <- mutate(res, seq = if_else(vaf > vaf_thd, alt, ref))  %>% 
	select(site, sample , seq) 

res_seq_normal <- mutate(seq, sample = "NORMAL") %>% 
	select(site, sample, seq = ref)

res_seq <- rbind(res_seq_normal, res_seq_tumor) %>% 
	arrange(site, sample) %>% 
	group_by(sample) %>% 
	reframe(seq = paste0(seq, collapse = "")) %>% 
	mutate(sample = paste0(">", sample)) %>% 
	t() %>% as.vector()  
length(res_seq)
writeLines(res_seq, con = paste0(res_dir, "/wgs_case1_cutoff0.05.fa"))


## run iqtree with modelfinder
paste0("iqtree2 "
		" -s ", res_dir, "/wgs_case1_cutoff0.05.fa",
		 " -o NORMAL -m MFP -T 12 -B 1000 ", 
		 " --prefix ", res_dir, "/wgs_case1_cutoff0.05",  
		 " --redo ") %>% 
	system()


## plot result
tree_file_path <- paste0(res_dir, "/wgs_case1_cutoff0.05.treefile")
data_iqtree <- read.iqtree(tree_file_path) %>% 
as_tibble() %>% 
mutate(region = case_when(label %in% c("PI_0", "PI_1", "PI_2") ~ "PI",
                          label %in% c("PO_12", "PO_14", "PO_10", "PO_b") ~ "PO",
                          label == "LM1_1" ~ "LM1",
                          label == "LM2_1" ~ "LM2",
                          label == "LM3_1" ~ "LM3",
                          grepl("LM4", label) ~ "LM4",
                          .default = label),
         tumor = if_else(grepl("LM", label), "metastasis", "primary")) %>%  
as.treedata()

ggtree(data_iqtree,
	 layout = "rectangular", ) +
	geom_tiplab(hjust = -0.02, align = F, size = 6) +
	geom_tippoint(aes(color = region), size = 3) +
	geom_nodepoint(aes(subset = UBoot > 60), color = "grey") +
	geom_nodepoint(aes(subset = UBoot > 90), color = "black") +
	theme(legend.position = c(0.1,0.8),
		  legend.key.size = 18,
		  text=element_text(size= 18 )) +
	theme(axis.title.x = element_text(hjust=0.3)) +
	ggtitle(label = paste0("vaf threshold: ",vaf_thd)) +
	scale_color_manual(values =  c("#A9D179", "brown", "#004529",  "pink",
								   "black", "steelblue", "#abddff")) +
	theme_tree2() 

ggsave(filename = paste0(res_dir, "wgs_case1_cutoff0.05_MLtree.pdf"),
	 width = 7, height = 8)



