setwd("~/project/crc/simulation/")
library(data.table)
library(parallel)
library(tidyverse)
source("~/software/script_tools/R_color.R")
source("~/software/script_tools/function_lib.R")


file_list <- list.files("res_simu/res_simulation_intro_tumor/", pattern = "Netural", full.names = T)
file_name <- file_list[10]

# get data ---------------------------------------------------------------------

data <- mclapply(file_list, mc.cores = 24, function(file_name){
  data <- fread(file_name)

  data <- lapply(1:nrow(data), function(i){
    temp <- data[i,]
    idx <- str_split(temp$cellmutation, ",", simplify = T) %>% as.vector()
    
    data.frame(x = temp$cellpos_x, y = temp$cellpos_y, idx)
  }) %>% 
    do.call(what = rbind)  
  
  data <- mutate(data, dataset = basename(file_name),
                 idx = if_else(idx %in% c("", NA), "0", idx),
                 migration_rate = str_extract(file_name, "(?<=_m-).*?(?=_)"),
                 mutation_rate = str_extract(file_name, "(?<=_u-).*(?=\\.txt)"))
  
  return(data)
}) %>% 
  do.call(what = rbind) %>% 
  mutate(across(contains("rate"), as.numeric))


head(data); dim(data)
unique(data$dataset) %>% length()

lapply(unique(data$mutation_rate), function(rate_m){
  temp <- filter(data, mutation_rate == rate_m)
  
  p1 <- ggplot(data = temp, aes(x = x, y = y, color = idx)) +
    geom_jitter(size = 0.25) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_wrap(~dataset, nrow = 4, ncol = 5) +
    coord_fixed()
  
  ggsave(plot = p1, 
         filename = paste0("res/mutation_",rate_m, "_mig.pdf"),
         width = 28, height = 25)
  })


# lineage by snvs --------------------------------------------------------------
data_plot <- filter(data,
                    migration_rate %in% c(0, 0.05, 0.15, 0.19),
                    mutation_rate %in% c(0.03, 0.06, 0.09)) %>% 
  arrange(idx) %>%
  group_by(migration_rate, mutation_rate, x, y) %>% 
  reframe(lineage = paste0(idx, collapse = "_"))    %>% 
  group_by(migration_rate, mutation_rate, lineage) %>%
  mutate(cell_num = length(x)) %>% 
  group_by(migration_rate, mutation_rate) %>% 
  arrange(desc(cell_num)) %>% 
  mutate(lineage_idx = factor(cell_num,
                              levels = unique(cell_num),
                              labels = paste0("r", 1:length(unique(cell_num))))) %>%
  arrange(migration_rate, mutation_rate, lineage_idx) %>% 
  group_by(migration_rate, mutation_rate) %>%
  mutate(cell_idx = 1:length(x))  

temp <- select(data_plot, migration_rate, mutation_rate, lineage_idx, cell_idx) %>% 
  group_by(migration_rate, mutation_rate,lineage_idx, ) %>% 
  mutate(cell_idx = min(cell_idx)) %>% 
  distinct %>% 
  group_by(migration_rate, mutation_rate) %>% 
  select(-cell_idx) %>% 
  nest(lineage_list = c(lineage_idx)  )


data_plot <- left_join(data_plot, temp, by = c("migration_rate", "mutation_rate")) %>% 
  filter(lineage_idx %in% unlist(lineage_list))


ggplot(data = data_plot, aes(x = x, y = y, color = lineage_idx)) +
  geom_jitter(size = 0.3) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(mutation_rate ~ migration_rate) +
  coord_fixed() +
  scale_color_manual(values = c(mycol22, mycol22_2))+
  theme(strip.text = element_text(size = 16),
        axis.text = element_text(size = 14))





