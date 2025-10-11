library(data.table)
library(tidyverse)
library(parallel)
source("script/function/calculate_moran_and_geary.R")


case_name <- "case1"

# SNV ----
snv_data_path <- paste0("data/amp_", case_name, ".rd/")

data_raw <- read_rds(snv_data_path) %>% 
  dplyr::filter(!region %in% c("LM1","LM2", "LM3")) %>% 
  mutate(region = case_when(region %in% c("PI", "PO", "PRU", "PRD", "PLU", "PLD") ~  "primary",
                            region %in% c("LM4", "LM") ~ "metastasis",
                            .default = "others"))

data <- select(data_raw, sample, region, x, y, z, starts_with("site")) %>% 
  pivot_longer(cols = starts_with("site"), 
               names_to = "site", 
               values_to = "vaf")
site_list <- unique(data$site)

test_tab <- expand.grid(site = site_list,
                        region = c("metastasis", "primary"),
                        k = 4:12,
                        stringsAsFactors = F)

res_moran_raw <- mclapply(1:nrow(test_tab), mc.cores = 36, function(i){
  temp <- dplyr::filter(data, 
                 site == test_tab$site[i],
                 region == test_tab$region[i])
  
  calculate_spatial_stats_3d(temp, 
                             feature_col = "vaf", 
                             k = test_tab$k[i]) %>% 
    mutate(k = test_tab$k[i],
           region = test_tab$region[i],
           site = test_tab$site[i],
           type = case_when(p.value <= 0.01 & statistic > expected ~ "concentrated",
                            p.value <= 0.01 & statistic < expected ~ "spread",
                            .default = "nosig"))
}) %>% 
  bind_rows()


vaf_thrd = 0.01
sample_num_thread = 5
p_value_thread <- 0.01
k_nums_thread <- 3

sample_num <- dplyr::filter(data, vaf >= vaf_thrd) %>% 
  reframe(sample_num = length(vaf),
          mean_vaf = mean(vaf),
          sd_vaf = sd(vaf),
          r50_vaf = quantile(vaf, 0.5),
          .by = c(site, region))

count(sample_num, region)

res_moran_full <- left_join(sample_num, res_moran_raw, by = c("site", "region"))  

res_moran_n <- reframe(res_moran_full,
                       k_num_sig = sum(p.value <= p_value_thread),
                       .by = c(site, region, sample_num, mean_vaf, sd_vaf, r50_vaf, idx)) %>% 
  pivot_wider(names_from = "idx",
              values_from = "k_num_sig",
              names_prefix = "sig_knum_") %>% 
  mutate(case = case_name,
         type = if_else(sig_knum_Moran  > k_nums_thread | sig_knum_Geary  > k_nums_thread ,
                        "Sig",
                        "notSig")) 



# pie plot ----
res <- select(res_moran_n, site, region, case, type) %>% 
  count(region, type) %>% 
  mutate(perc =  sprintf("%.2f%%\r\n(n=%d)", n/sum(n) * 100, n),
          sum_sample = sum(n),
          .by = c(region))

ggplot(res, aes(x =   'Content', y = n, fill = type)) +
  geom_bar(stat = 'identity', position = 'stack') +
  geom_text(aes(label = perc), size =8) +
  facet_wrap(~ region, scales = "free")+
  coord_polar(theta = 'y') +
  theme_minimal() +
  theme(strip.text = element_text(size = 26),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  scale_fill_manual(values = c(Sig = "darkred", notSig = "steelblue")) +
  ggtitle("Spatial Autocorrelation")


res <- bind_rows(res1, res2)



# sanki plot ----
library(ggalluvial)

data_plot <- dplyr::filter(res_moran_n, sample_num >= sample_num_thread) %>% 
  select(site, region, type) %>% 
  pivot_wider(names_from = "region", values_from = "type") %>% 
  mutate(type = case_when(is.na(metastasis) & !is.na(primary) ~ "primary_specific",
                          !is.na(metastasis) & is.na(primary) ~ "metastasis_specific",
                          !is.na(metastasis) & !is.na(primary) ~ "shared"),
         across(c(metastasis, primary), function(x){x[is.na(x)] = "notExist"; return(x) }) 
  ) %>% 
  arrange(type, metastasis, primary) %>% 
  mutate(Freq = 1)





is_alluvia_form(data_plot, weight ="Freq")


df_lodes <- to_lodes_form(data_plot,
                          key ="x", 
                          value = "stratum", 
                          id = "alluvium",
                          axes =2:3) %>% 
  mutate(stratum = factor(stratum, levels = c("Sig", "notSig", "notExist")),
         x = factor(x, levels = c("primary", "metastasis")))

is_lodes_form(df_lodes,key = "x",value = "stratum",id = "alluvium",weight ="Freq")



ggplot(df_lodes,aes(x = x, stratum =stratum, alluvium = alluvium,
                    fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0.32, knot.pos = 0.2) +
  geom_stratum(alpha = .9, color="grey20", width = 1/3) +
  geom_text(stat = "stratum", size =6,color="black") +
  scale_fill_manual(values = c(Sig = "darkred", 
                               notSig = "steelblue", 
                               notExist = "darkgreen")) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),
        axis.ticks =element_blank(),
        axis.text.y =element_blank(),
        axis.text.x = element_text(size = 18))+
  guides(fill = FALSE)

