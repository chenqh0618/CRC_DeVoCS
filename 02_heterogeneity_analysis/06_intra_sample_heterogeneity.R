library(tidyverse)
library(parallel)
library(data.table)
library(agricolae)
library(ggpubr)
source("function_lib/3D_plot.R")
source("function_lib/calculate_moran_and_geary.R")

calculate_shannon_vaf <- function(vafs, method = "binning", breaks = seq(0, 1, 0.05)) {
  vafs <- vafs[!is.na(vafs)]  
  
  if (method == "unique") {
    counts <- table(vafs)
  } else if (method == "binning") {
    vaf_bins <- cut(vafs, breaks = breaks, include.lowest = TRUE, right = FALSE)
    counts <- table(vaf_bins)
  } 
  
  proportions <- counts / sum(counts)
  proportions <- proportions[proportions > 0] # AVOID log(0)
  
  if (length(proportions) <= 1) {
    return(0.0)
  }
  
  shannon_index <- -sum(proportions * log(proportions)) 
  
  return(shannon_index)
}


calculate_vaf_sd <- function(vafs) {
  vafs <- vafs[!is.na(vafs)] 
  if (length(vafs) < 2) {
    return(NA)
  }
  return(stats::sd(vafs))
}



# case1 ----
## calculate ITH
data_raw <- read_rds("/public/data/cqh_project/crc/amp_data/filtered_CASE1.rds")

res_vaf_ith <- select(data_raw, sample, starts_with("site")) %>% 
  pivot_longer(cols = starts_with("site"), names_to = "site", values_to = "vaf") %>% 
  reframe(vaf_sd = calculate_vaf_sd(vaf),
		  vaf_shannon = calculate_shannon_vaf(vaf, 
                                              method = "binning",
                                              breaks = seq(0, 1, 0.05)),
          .by = c(sample, region, tumor, x, y, z))  

ggplot(data = res_vaf_ith, aes(x = region, y = vaf_shannon, fill = tumor)) +
  geom_boxplot()  

## ANOVA
anova_result <- aov(vaf_shannon ~ region, data = res_vaf_ith)
lsd_result <- LSD.test(anova_result, "region", p.adj = "none",alpha = 0.01)
print(lsd_result)



## spatial autocorrelation analysis

data <- filter(data_raw,!region %in% c("LM1", "LM2", "LM3")) 

test_tab <- expand.grid(tumor = c("metastasis", "primary"),
                        k = 4:12,
                        stringsAsFactors = F)

res_moran_raw <- mclapply(1:nrow(test_tab), mc.cores = 36, function(i){
  temp <- dplyr::filter(data, 
                        site == test_tab$site[i],
                        tumor == test_tab$tumor[i])
  res <- NULL
  try({
    res <-   calculate_spatial_stats_3d(temp, 
                                        feature_col = "vaf_shannon", 
                                        k = test_tab$k[i]) %>% 
      mutate(k = test_tab$k[i],
             case = test_tab$case[i],
             tumor = test_tab$tumor[i],
             site = test_tab$site[i],
             type = case_when(p.value <= 0.01 & statistic > expected ~ "concentrated",
                              p.value <= 0.01 & statistic < expected ~ "spread",
                              .default = "nosig"))
  })
  return(res)
}) %>% 
  bind_rows()

res <- filter(res_moran_raw, idx == "Moran") %>% 
  arrange(type) %>% 
  count(type, tumor, case)



## 3D plot
plot_data <- list(occurrens = list(df = res_vaf_ith,
								  size_col = NULL,
								  color_col = "vaf_shannon",
								  point_default_color = NULL,
								  colorscale = "RdBu",
								  marker_opacity = 0.8,
								  marker_size_min = 5,
								  marker_size_max = 7,
								  show_colorbar = T,
								  show_legend_points = T,
								  point_trace_name = "occurrens"))  

fig <- plot_wireframe_sphere(sphere_list = list(LM1 = c(2.5, 2.5, 15, 2, sphere_diag),
												LM2 = c(7, 3, 18, 2, sphere_diag),
												LM3 = c(3.5, 7.5, 21.5, 3, sphere_diag),
												LM4 = c(12, 12, 23, 5, sphere_diag),
												primary = c(7.1, 6.8, 3, 5.7, sphere_diag)),
							line_color = "#BCBCBC",
							line_width = 1, 
							points_df = plot_data,
							title = site_name) %>% 
							  layout(scene = list(
								camera = list(
								  up = list(x = 0, y = 0, z = 1),       
								  center = list(x = 0, y = 0, z = 0),   
								  eye = list(x = -2, y = 2, z = 0) 
								)
							  ))




# case2 ----
## calculate ITH
data_raw <- read_rds("public/data/cqh_project/crc/amp_data/filtered_CASE2.rds")
data_vaf_ith <- select(data_raw, sample, starts_with("site")) %>% 
  pivot_longer(cols = starts_with("site"), names_to = "site", values_to = "vaf") %>% 
  reframe(vaf_sd = calculate_vaf_sd(vaf),
		  vaf_shannon = calculate_shannon_vaf(vaf, 
                                              method = "binning",
                                              breaks = seq(0, 1, 0.05)),
          .by = c(sample, region, tumor, x, y, z))  
ggplot(data = res_vaf_ith, aes(x = region, y = vaf_shannon, fill = tumor)) +
  geom_boxplot()  


## ANOVA
anova_result <- aov(vaf_shannon ~ region, data = res_vaf_ith)
lsd_result <- LSD.test(anova_result, "region", p.adj = "none",alpha = 0.01)
print(lsd_result)



## spatial autocorrelation analysis
test_tab <- expand.grid(tumor = c("metastasis", "primary"),
                        k = 4:12,
                        stringsAsFactors = F)

res_moran_raw <- mclapply(1:nrow(test_tab), mc.cores = 36, function(i){
  temp <- dplyr::filter(data_raw, 
                        site == test_tab$site[i],
                        tumor == test_tab$tumor[i])
  res <- NULL
  try({
    res <-   calculate_spatial_stats_3d(temp, 
                                        feature_col = "vaf_shannon", 
                                        k = test_tab$k[i]) %>% 
      mutate(k = test_tab$k[i],
             case = test_tab$case[i],
             tumor = test_tab$tumor[i],
             site = test_tab$site[i],
             type = case_when(p.value <= 0.01 & statistic > expected ~ "concentrated",
                              p.value <= 0.01 & statistic < expected ~ "spread",
                              .default = "nosig"))
  })
  return(res)
}) %>% 
  bind_rows()

res <- filter(res_moran_raw, idx == "Moran") %>% 
  arrange(type) %>% 
  count(type, tumor, case)

## 3D plot
plot_data <- list(occurrens = list(df = data_vaf_ith,
								  size_col = NULL,
								  color_col = "vaf_shannon",
								  point_default_color = NULL,
								  colorscale = "RdBu",
								  marker_opacity = 0.8,
								  marker_size_min = 5,
								  marker_size_max = 7,
								  show_colorbar = T,
								  show_legend_points = T,
								  point_trace_name = "occurrens"))

fig <- plot_wireframe_sphere(sphere_list = list(primary = c(3, 3, 3, 5, sphere_diag),
												metastasis = c(3, 3, 16, 4.7, sphere_diag)),
												
							line_color = "#BCBCBC",
							line_width = 0.65,
							points_df = plot_data,
							title = site_name) %>%
							  layout(scene = list(
								camera = list(
								  up = list(x = 0, y = 0, z = 1),        
								  center = list(x = 0, y = 0, z = 0),   
								  eye = list(x = -2, y = 2, z = 0)  
								)
							  )

