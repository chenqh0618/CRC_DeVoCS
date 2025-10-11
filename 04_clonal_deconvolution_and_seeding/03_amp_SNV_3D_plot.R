library(tidyverse)
library(data.table)
library(parallel)
library(vcfR)

setwd("~/project/crc/combined/")
source("~/software/script_tools/R_color.R")
vaf_thred = 0.05

# case1 ------------------
## read data ----
amp_path <- "res/res_amp_snv/case1_processed/s1_DepthAndNA_filtered_data.rds"
xyz_path <- "res/res_amp_snv/case1_processed/case1_combined.rds"



data_amp <- read_rds(amp_path) %>% 
  mutate(vaf =alt / total) %>% 
  mutate(site = paste0("chr", site)) %>% 
  filter(!sample %in% c("Normal", "normal1", "normal2", "X0", "X1")) %>% 
  filter(mean(total) > 1000, .by = c(site)) %>%
  arrange(region)
data_xyz <- read_rds(xyz_path) %>% 
  distinct(sample, region, tumor, x, y, z )  

data_amp <- left_join(data_amp, data_xyz) 



# site idx
region_tab <- dplyr::count(data_amp, region) %>% 
  mutate(idx = 2^(1:nrow(.)),
         n = n/length(unique(data_amp$site))) %>% 
  select(region, idx, region_sample_num = n)


data <- filter(data_amp, vaf > vaf_thred) %>%
  na.omit() %>% 
  dplyr::count(tumor, region, site) %>%
  rename( region_snv_sample_num = n) %>% 
  filter( region_snv_sample_num >= 5) %>%
  filter(length(unique(tumor))> 1, .by = c("site")) %>%
  left_join(region_tab, by = c("region")) 

site_tab <- reframe(data, idx = sum(idx), .by = site) %>% 
  arrange(idx)

  

## heatmap ----
vaf_thred = 0.05
data_plot <- filter(data_amp, 
                    # region %in% c("PI", "PO", "LM4"),
                    site %in% site_tab$site) %>% 
  select(site, sample, vaf) %>% 
  pivot_wider(names_from = sample, values_from = vaf) %>% 
  column_to_rownames("site")
data_plot <- (data_plot > vaf_thred) * 1
pheatmap::pheatmap(data_plot[total_snv_list$site,], 
                   show_rownames = T, show_colnames = T,
                   clustering_method = "ward.D2",
                   cluster_cols = F, cluster_rows = T,
                   treeheight_row = 0,
                   annotation_col = distinct(data_amp, sample, region) %>% 
                     column_to_rownames("sample"),
                   # annotation_row = select(total_snv_list, site, type) %>% 
                   #   column_to_rownames("site"),
                   annotation_colors = list(region = case1_color_list,
                                            type = c(type1 = "#295788",
                                                     type2 = "#0F4928",
                                                     type3 = "#A81F24")), 
                   # gaps_row = c(8, 28),
                   color = colorRampPalette(c( "white", "darkred"))(50))




site1 <- c("1_216595546")
site2 <- c("12_120575106","14_94927140", "19_48247924","7_1097597")
site3 <- c("10_91460743")
site4 <- c("X_139867384", "16_29755618")




## 3D point ----
source("script/function/plot_wireframe_sphere.R")

library(plotly)
site_name = "chr7_72171210" #case1,selected
site_name = "chr1_93619289" #case1, selected


sphere_diag <- 12
point_size = 10
vaf_thred <- 0.05

df <- list(bg = list(df = filter(data_amp, 
                                 site == site_name,
                                 vaf <= vaf_thred),
                     size_col = NULL,
                     color_col = NULL,
                     point_default_color = "grey",
                     colorscale = NULL,
                     marker_opacity = 0.6,
                     marker_size_min = point_size-1,
                     marker_size_max = point_size+1,
                     show_colorbar = T,
                     show_legend_points = T),
           oc1 = list(df = filter(data_amp,
                                 site == site_name,
                                  vaf > vaf_thred),
                     size_col = NULL,
                      color_col = NULL,
                      point_default_color = "darkred",
                      colorscale = NULL,
                      marker_opacity = 1,
                      marker_size_min = point_size-1,
                      marker_size_max = point_size+1,
                      show_colorbar = T,
                      show_legend_points = T) ) 
try({
  fig <- plot_wireframe_sphere(sphere_list = list(LM1 = c(2.5, 2.5, 15, 2, sphere_diag),
                                                  LM2 = c(7, 3, 18, 2, sphere_diag),
                                                  LM3 = c(3.5, 7.5, 21.5, 3, sphere_diag),
                                                  LM4 = c(12, 12, 23, 5, sphere_diag),
                                                  # PI = c(6, 6, 3, 5, sphere_diag),
                                                  # PO = c(8.5, 8, 3, 5, sphere_diag),
                                                  primary = c(7.1, 6.8, 3, 5.7, sphere_diag)
  ),
  line_color = "#BCBCBC",
  line_width = 1, 
  points_df = df,
  title = site_name) %>% 
    layout(scene = list(
      camera = list(
        up = list(x = 0, y = 0, z = 1),        
        center = list(x = 0, y = 0, z = 0),   
        eye = list(x = -2, y = 2, z = 0) 
      )
    ))
  # fig
  file_path <- "res_case1_snv1.png"
  plotly::save_image(p = fig, file = file_path, height = 1200, width = 900)
   })


# case2 ---------


## read data ----
amp_path <- "res/res_amp_snv/case2_processed/s1_DepthAndNA_filtered_data.rds"
xyz_path <- "res/res_amp_snv/case2_processed/case2_combined.rds"



data_amp <- read_rds(amp_path) %>% 
  mutate(vaf =alt / total) %>% 
  mutate(site = paste0("chr", site)) %>% 
  filter(!sample %in% c("Normal","normal", "normal1", "normal2", "X0", "X1")) %>% 
  filter(mean(total) > 1000, .by = c(site)) %>%
  arrange(region)
data_xyz <- read_rds(xyz_path) %>% 
  distinct(sample, region, tumor, x, y, z )  

data_amp <- left_join(data_amp, data_xyz)



# site idx
region_tab <- dplyr::count(data_amp, region) %>% 
  mutate(idx = 2^(1:nrow(.)),
         n = n/length(unique(data_amp$site))) %>% 
  select(region, idx, region_sample_num = n)


data <- filter(data_amp, vaf > vaf_thred) %>% 
  na.omit() %>%
  dplyr::count(tumor, region, site) %>%
  rename( region_snv_sample_num = n) %>% 
  filter( region_snv_sample_num >= 5) %>%
  filter(length(unique(tumor))> 1, .by = c("site")) %>%
  left_join(region_tab, by = c("region")) 

site_tab <- reframe(data, idx = sum(idx), .by = site) %>% 
  arrange(idx)




## heatmap ----
vaf_thred = 0.05
data_plot <- filter(data_amp, 
                    site %in% site_tab$site) %>% 
  select(site, sample, vaf) %>% 
  pivot_wider(names_from = sample, values_from = vaf) %>% 
  column_to_rownames("site")
data_plot <- (data_plot > vaf_thred) * 1
pheatmap::pheatmap(data_plot[site_tab$site,], 
                   show_rownames = T, show_colnames = T,
                   clustering_method = "ward.D",
                   cluster_cols = F, cluster_rows = T,
                   treeheight_row = 0,
                   annotation_col = distinct(data_amp, sample, region) %>% 
                     column_to_rownames("sample"),
                   # annotation_row = select(total_snv_list, site, type) %>% 
                   #   column_to_rownames("site"),
                   annotation_colors = list(region = case2_color_list,
                                            type = c(type1 = "#295788",
                                                     type2 = "#0F4928",
                                                     type3 = "#A81F24")), 
                   # gaps_row = c(8, 28),
                   color = colorRampPalette(c( "white", "darkred"))(50))





## 3D point ----
source("script/function/plot_wireframe_sphere.R")
library(plotly)

site_name = "chr8_53534127" #case2,selected
site_name <- "chr14_42076739" #case2, selected


sphere_diag <- 12
point_size = 12
vaf_thred <- 0.05

df <- list(bg = list(df = filter(data_amp, 
                                  site == site_name,
                                  vaf <= vaf_thred),
                      size_col = NULL,
                      color_col = NULL,
                      point_default_color = "grey",
                      colorscale = NULL,
                      marker_opacity = 0.6,
                      marker_size_min = point_size-1,
                      marker_size_max = point_size+1,
                      show_colorbar = T,
                      show_legend_points = T),
           oc1 = list(df = filter(data_amp, 
                                  site == site_name,
                                  vaf > vaf_thred),
                      size_col = NULL,
                      color_col = NULL,
                      point_default_color = "darkred",
                      colorscale = NULL,
                      marker_opacity = 1,
                      marker_size_min = point_size-1,
                      marker_size_max = point_size+1,
                      show_colorbar = T,
                      show_legend_points = T)
           ) 
try({
  fig <- plot_wireframe_sphere(sphere_list = list(primary = c(3, 3, 3, 5, sphere_diag),
                                                  metastasis = c(3, 3, 16, 4.7, sphere_diag)
  ),
  line_color = "#BCBCBC",
  line_width = 1, 
  points_df = df,
  title = site_name) %>% 
    layout(scene = list(
      camera = list(
        up = list(x = 0, y = 0, z = 1),       
        center = list(x = 0, y = 0, z = 0),   
        eye = list(x = -2, y = 2, z = 0) 
      )
    ))
   # fig
  file_path <- "res_case2_snv1.png"
  plotly::save_image(p = fig, file = file_path, height = 1200, width = 1000)
})



