pjdir <- "/public/data/cqh_project/crc/amp_case1/"
setwd(pjdir)
library(tidyverse)
library(data.table)
library(parallel)
library(treedataverse)
library(poppr)
library(ggpmisc)


source("~/software/script_tools/R_color.R")
source("~/software/script_tools/function_lib.R")

nor_min_max=function(x){
  y=na.omit(x)
  return((x - min(y))/(max(y) - min(y)))
}

vaf_thred = 0.05
# GET DATA ---------------------------------------------------------------------
data <- read_rds("data/processed/s1_DepthAndNA_filtered_data.rds")

sample_info <- fread("data/raw/sample_info.csv") %>% 
  mutate(z = if_else(region  %in% c("N", "W"), z *0.03, z)) %>% 
  select(sample ,  x, y, z) %>% 
  filter(sample != "", !is.na(sample))

data_snv <- left_join(data, sample_info, by = "sample") %>% 
  mutate(maf = if_else(maf <= 0.01, 0, maf),
         region2 = if_else(region %in% c("PI", "PO"), "CRC", region))



data_snv <- group_by(data_snv, region2) %>% 
  mutate(x = nor_min_max(x),
          y = nor_min_max(y),
          z = nor_min_max(z)) %>% 
  filter(region %in% c("PI", "PO", "LM4")) %>% 
  ungroup()

group_by(data_snv, region) %>% 
  reframe(x = range(x),
          y = range(y),
          z = range(z))

## SPATIAL DISTANCE ----
d3_spatial_dist <- select(data_snv, sample, x, y, z) %>% 
  distinct() %>% 
  column_to_rownames("sample") %>% 
  dist(method = "euclidean") %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("from") %>% 
  pivot_longer(cols = -from, names_to = "to", values_to = "d3_dist")


dim(d3_spatial_dist)

temp <- select(data_snv, sample, region, z) %>% distinct()
spatial_dist <- left_join(d3_spatial_dist,
            select(temp, sample, region_from = region, z_from = z),
            by = c("from" = "sample")) %>% 
  left_join(select(temp, sample, region_to = region, z_to = z),
            by = c("to" = "sample"))



# GENETIC DISTANCE  -----
### fa distance
#edwards.dist, nei.dist, rogers.dist, reynolds.dist, prevosti.dist
gene_data <- read.FASTA("data/backup//raw_cutoff0.05.fa")
gene_data <- DNAbin2genind(gene_data)


nei_dist <- nei.dist(gene_data) %>% 
  as.matrix() %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  pivot_longer(cols = -from, values_to = "nei_dist",
               names_to = "to") %>% 
  filter(nei_dist != 0)
head(nei_dist); dim(nei_dist)

### maf distance
af <- select(data, site, sample, maf) %>% 
  pivot_wider(id_cols = sample, names_from = site, values_from = maf, values_fill = 0) %>% 
  column_to_rownames("sample")

maf_dist <- dist(af, method = "euclidean") %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("from") %>% 
  pivot_longer(cols = -from, names_to = "to", values_to = "maf_dist") %>% 
  filter(from !=  to)

genetic_dist <- full_join(nei_dist, maf_dist, by = c("from", "to"), )
head(genetic_dist); dim(genetic_dist)


ggplot(genetic_dist, aes(x = nei_dist, y = maf_dist))+
  geom_point()

## merge data ----
data_raw <- left_join(spatial_dist, genetic_dist, by = c("from", "to"))
dim(data_raw); head(data_raw)


# CORRELATION IN 3D SPACE ------------------------------------------------------
## Prmary
data_plot <- filter(data_raw,
                    region_from == region_to, 
                    region_to %in% c("PI", "PO"))

my.formula <- y ~ x
ggplot(data = data_plot, 
       aes( x = d3_dist, y = maf_dist )) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  geom_jitter(aes( color = region_to), size = 2, alpha = 0.8) + 
  geom_smooth(aes( color = region_to), method = "lm", linewidth = 2 )+ 
  geom_smooth(method = "lm", level = 0.95, formula = my.formula,
              color="black", linewidth = 2)+ 
  
  ggpubr::stat_cor(aes( color = region_to), size = 10)  +
  ggpubr::stat_cor( size = 10)  +
  theme_classic() +
  theme(axis.text = element_text(size =16 ),
        title = element_text(size = 24)) +
  scale_color_manual(values = c(Greens[c(3,8)])) +
  ggtitle("samples in CRC of case1")


## LM4

data_plot <- filter(data_raw,
                    region_from == region_to, 
                    region_to %in% c("LM4"))

my.formula <- y ~ x
ggplot(data = data_plot, 
       aes( x = d3_dist, y = maf_dist)) +
  geom_jitter(size = 2, alpha = 0.6, color = "darkred") + 
  geom_smooth(method = "lm",  color="black", linewidth = 2)+ 
  ggpubr::stat_cor(size = 10)  +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               parse = TRUE) +
  theme_classic() +
  theme(axis.text = element_text(size =16 ),
        title = element_text(size = 24)) +
  scale_color_manual(values = c(Darkreds[9])) +
  ggtitle("samples in LM4 of case1")



