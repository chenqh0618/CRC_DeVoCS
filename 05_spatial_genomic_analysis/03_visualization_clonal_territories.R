library(tidyverse)
library(data.table)
library(parallel)
library(vcfR)
library(reticulate)
library(plotly)

setwd("~/project/crc/combined/")
source("~/software/script_tools/R_color.R")
vaf_thred = 0.05


# case1 ----
## data input ----
amp_path <- "res/res_amp_snv/case1_processed/s1_DepthAndNA_filtered_data.rds"
xyz_path <- "res/res_amp_snv/case1_processed/case1_combined.rds"

data_amp <- read_rds(amp_path) %>% 
  mutate(vaf =alt / total) %>% 
  mutate(site = paste0("chr", site)) %>% 
  filter(!sample %in% c("Normal", "normal1", "normal2", "X0", "X1")) %>% 
  # filter(median(total) > 100, .by = c(site)) %>%
  arrange(region)
data_xyz <- read_rds(xyz_path) %>% 
  distinct(sample, region, tumor, x, y, z )  

data_amp <- left_join(data_amp, data_xyz) 



region_tab <- dplyr::count(data_amp, region) %>% 
  mutate(idx = 2^(1:nrow(.)),
         n = n/length(unique(data_amp$site))) %>% 
  select(region, idx, region_sample_num = n)


data <- filter(data_amp, vaf > vaf_thred) %>%
  na.omit() %>% 
  dplyr::count(tumor, region, site) %>%
  rename( region_snv_sample_num = n) %>% 
  filter( region_snv_sample_num >= 5) %>%
  left_join(region_tab, by = c("region")) 

site_tab <- reframe(data, idx = sum(idx), .by = site) %>% 
  arrange(idx) %>% 
  filter(idx >= 64)
table(site_tab$idx)



## heatmap ----
vaf_thred = 0.05
data_plot <- filter(data_amp, 
                    region %in% c("LM4"),
                    site %in% site_tab$site) %>% 
  select(site, sample, vaf) %>% 
  pivot_wider(names_from = sample, values_from = vaf) %>% 
  column_to_rownames("site")
data_plot <- (data_plot > vaf_thred) * 1


temp <- pheatmap::pheatmap(data_plot[site_tab$site,], 
                   show_rownames = T, show_colnames = T,
                   clustering_method = "ward.D",
                   cluster_cols = T, cluster_rows = T,
                   treeheight_row = 0, treeheight_col = 0,
                   # annotation_col = distinct(data_amp, sample, region) %>% 
                   #   column_to_rownames("sample"),
                   annotation_row = select(site_tab, site, idx) %>%
                     mutate(idx = as.character(idx)) %>% 
                     column_to_rownames("site"),
                   annotation_colors = list(region = case1_color_list,
                                            type = c(type1 = "#295788",
                                                     type2 = "#0F4928",
                                                     type3 = "#A81F24")), 
                   # gaps_row = c(8, 28),
                   color = colorRampPalette(c( "white", "darkred"))(50))


rownames(data_plot)

site_list <- temp$tree_row$labels[temp$tree_row$order] 


head(site_tab)
res_tab <- lapply(2^(1:6), function(i){
  contains_number(site_tab$idx, i)*1
}) %>% 
  do.call(what=cbind) %>% 
  `colnames<-`(region_tab$region) %>%
  `rownames<-`(site_tab$site)
res_tab <- res_tab[site_list,]

corrplot::corrplot(res_tab, is.corr = F, 
                   method = "circle", 
                   type = "full", 
                   order = "original")

## plot 3D ----

source("script/function/plot_wireframe_sphere.R")
res_dir <- "res/res_amp_vaf_distribution/"
library(plotly)

site_list <- c("chr19_49661058","chrX_139867384", "chr1_111738563", "chr1_216595546")
color_list <- c("steelblue", "darkred",  "#49A822", "#FFA500")
names(color_list) <- site_list
site_name <- "chr19_48247924"


sphere_diag <- 12
point_size = 12
vaf_thred <- 0.05
lapply(site_list, function(site_name){
  df <- list(bg = list(df = filter(data_amp,
                                   site == site_name,
                                   vaf <= vaf_thred,
                                   region %in% c("LM4")),
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
                                    vaf > vaf_thred,
                                    region %in% c("LM4")),
                        size_col = NULL,
                        # color_col = "vaf",
                        # colorscale = "Reds",
                        point_default_color = color_list[site_name],
                        marker_opacity = 1,
                        marker_size_min = point_size-1,
                        marker_size_max = point_size+1,
                        show_colorbar = T,
                        show_legend_points = T) ) 
  
  
  try({
    fig <- plot_wireframe_sphere(
      sphere_list = list(LM4 = c(12, 12, 23, 5, sphere_diag)),
      line_color = "#BCBCBC",
      line_width = 1,
      points_df = df,
      title = NULL
      ) %>%
      layout(
        scene =  list(
          camera = list(
            up = list(x = 0, y = 0, z = 1),
            # Z 轴朝上
            center = list(x = 0, y = 0, z = 0),
            # 相机看向原点
            eye = list(x = -2, y = 2, z = 0) # 另一个可能的45度视角（从第二象限看）
            # eye = list(x = 2.1, y = -2.1, z = 0) # 另一个可能的45度视角（从第四象限看）
          )
        ),
        margin = list(
          l = 1,
          r = 1,
          b = 1,
          t = 1
        ),
        # 调整边距
        legend = list(
          orientation = "h",
          # 水平图例
          xanchor = "center",
          x = 0.5,
          y = -0.1,
          # 图例位置在底部中间
          font = list(size = 10)
        )
        # 调整图例字体大小
      )
                                 
    # fig
    file_path <- paste0(res_dir, "/case1_single_site_", site_name, "_intratumor_spatial.jpg")
    plotly::save_image(p = fig, file = file_path, height = 800, width = 800)
  })
  
  
})


## plot 2D ----
data_bg <- filter(data_amp, site %in% site_list, region %in% c("LM4")) %>%
  filter(max(vaf, na.rm = T) <= vaf_thred, .by = sample) %>% 
  select(region, tumor, sample, x, y, z) %>% 
  distinct() %>%
  mutate(vaf = 0, site = "bg") %>% 
  na.omit() 
  
jetter_distance = 0
data_plot <- filter(data_amp, 
                     site %in% site_list,
                     region %in% c("LM4")) %>% 
  select(-c(temp, total, alt)) %>% 
  filter(vaf > vaf_thred) %>% 
  bind_rows(data_bg) %>% 
  na.omit() %>% 
  mutate(z = factor(z, 
                    levels= unique(.$z) %>% sort(decreasing = T),
                    labels = paste0("layer", sprintf("%02d", 1:length(unique(.$z))))),
         x = case_when(site == site_list[2] ~ x+jetter_distance,
                       site == site_list[3] ~ x-jetter_distance,
                       site == site_list[4] ~ x-jetter_distance,
                       .default = x),
         y = case_when(site == site_list[2] ~ y+jetter_distance,
                       site == site_list[3] ~ y+jetter_distance,
                       site == site_list[4] ~ y-jetter_distance,
                       .default = y))   %>% 
  filter(!z %in% c("layer01", "layer14"))

filter(data_plot, duplicated(sample))

site_list <- c("chr19_49661058","chrX_139867384", "chr1_111738563", "chr1_216595546")
color_list <- c("steelblue", "darkred",  "darkgreen", "#FFA500")
names(color_list) <- site_list
color_list <- c(color_list, bg = "grey")

plot_case1_2d <- ggplot(data = data_plot, 
       aes(x=x, y = y)) +
  geom_jitter(aes(color = site ), size = 4.8) +
  geom_circle(aes(x0 = 11.1, y0 = 12, r = 4.4),
              linetype =2)+
  facet_wrap(~z) +
  theme_bw() +
  scale_color_manual(values = color_list) +
  ylim(c(16.4, 7)) +
  xlim(c(16, 6.65)) +
  coord_fixed(ratio = 1) +
  theme(axis.text = element_text(size = 16),
        strip.text = element_text(size = 21),
        # legend.text = element_text(size = 16)
  ) 
plot_case1_2d





# case2 ----
## data input ----
amp_path <- "res/res_amp_snv/case2_processed/s1_DepthAndNA_filtered_data.rds"
xyz_path <- "res/res_amp_snv/case2_processed/case2_combined.rds"



data_amp <- read_rds(amp_path) %>% 
  mutate(vaf =alt / total) %>% 
  mutate(site = paste0("chr", site)) %>% 
  filter(!sample %in% c("Normal", "normal1", "normal2", "X0", "X1")) %>% 
  arrange(region)
data_xyz <- read_rds(xyz_path) %>% 
  distinct(sample, region, tumor, x, y, z )  

data_amp <- left_join(data_amp, data_xyz) 


region_tab <- dplyr::count(data_amp, region) %>% 
  mutate(idx = 2^(1:nrow(.)),
         n = n/length(unique(data_amp$site))) %>% 
  select(region, idx, region_sample_num = n)


data <- filter(data_amp, vaf > vaf_thred) %>%
  na.omit() %>% 
  dplyr::count(tumor, region, site) %>%
  rename( region_snv_sample_num = n) %>% 
  filter( region_snv_sample_num >= 3) %>%
  left_join(region_tab, by = c("region")) 

site_tab <- reframe(data, idx = sum(idx), .by = site) %>% 
  arrange(idx) %>% 
  filter(idx >= 32)
table(site_tab$idx)



## heatmap ----
vaf_thred = 0.05
data_plot <- filter(data_amp, 
                    region %in% c("LM"),
                    site %in% site_tab$site) %>% 
  select(site, sample, vaf) %>% 
  pivot_wider(names_from = sample, values_from = vaf) %>% 
  column_to_rownames("site")
data_plot <- (data_plot > vaf_thred) * 1
temp <- pheatmap::pheatmap(data_plot[site_tab$site,], 
                           show_rownames = T, show_colnames = T,
                           clustering_method = "ward.D2",
                           cluster_cols = T, cluster_rows = T,
                           treeheight_row = 0, treeheight_col = 0,
                           # annotation_col = distinct(data_amp, sample, region) %>% 
                           #   column_to_rownames("sample"),
                           annotation_row = select(site_tab, site, idx) %>%
                             mutate(idx = as.character(idx)) %>% 
                             column_to_rownames("site"),
                           annotation_colors = list(type = c(type1 = "#295788",
                                                             type2 = "#0F4928",
                                                             type3 = "#A81F24")), 
                           # gaps_row = c(8, 28),
                           color = colorRampPalette(c( "white", "darkred"))(50))


rownames(data_plot)

site_list <- temp$tree_row$labels[temp$tree_row$order] 

res_tab <- lapply(2^(1:5), function(i){
  contains_number(site_tab$idx, i)*1
}) %>% 
  do.call(what=cbind) %>% 
  `colnames<-`(region_tab$region) %>%
  `rownames<-`(site_tab$site)
res_tab <- res_tab[site_list,]

corrplot::corrplot(res_tab, is.corr = F, 
                   method = "circle", 
                   type = "full", 
                   order = "original")

## plot 3D ----
source("script/function/plot_wireframe_sphere.R")
res_dir <- "res/res_amp_vaf_distribution/"
library(plotly)
site_name = "chr11_56408815" #case1,selected 
site_name = "chr18_24039766" #case1, selected
site_name = "chr8_53534127" #case1, selected
site_name = "chr4_158094116" #case1, selected 

site_list <- c("chr11_56408815","chr18_24039766", "chr20_31669258", "chr4_158094116")
color_list <- c("steelblue", "darkred",  "#49A822", "#FFA500")
names(color_list) <- site_list
site_name <- "chr19_48247924"


sphere_diag <- 12
point_size = 10
vaf_thred <- 0.05
lapply(site_list, function(site_name){
  df <- list(bg = list(df = filter(data_amp,
                                   site == site_name,
                                   vaf <= vaf_thred,
                                   region %in% c("LM")),
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
                                    vaf > vaf_thred,
                                    region %in% c("LM")),
                        size_col = NULL,
                        # color_col = "vaf",
                        # colorscale = "Reds",
                        point_default_color = color_list[site_name],
                        marker_opacity = 1,
                        marker_size_min = point_size-1,
                        marker_size_max = point_size+1,
                        show_colorbar = T,
                        show_legend_points = T) ) 
  
  
  try({
    fig <- plot_wireframe_sphere(
      sphere_list = list(LM = c(3, 3, 16.5, 4.2, sphere_diag)),
      line_color = "#BCBCBC",
      line_width = 1,
      points_df = df,
      title = NULL
    ) %>%
      layout(
        scene =  list(
          camera = list(
            up = list(x = 0, y = 0, z = 1),
            # Z 轴朝上
            center = list(x = 0, y = 0, z = 0),
            # 相机看向原点
            eye = list(x = -2, y = 2, z = 0) # 另一个可能的45度视角（从第二象限看）
            # eye = list(x = 2.1, y = -2.1, z = 0) # 另一个可能的45度视角（从第四象限看）
          )
        ),
        margin = list(
          l = 1,
          r = 1,
          b = 1,
          t = 1
        ),
        # 调整边距
        legend = list(
          orientation = "h",
          # 水平图例
          xanchor = "center",
          x = 0.5,
          y = -0.1,
          # 图例位置在底部中间
          font = list(size = 10)
        )
        # 调整图例字体大小
      )
    
    # fig
    file_path <- paste0(res_dir, "/case2_single_site_", site_name, "_intratumor_spatial.jpg")
    plotly::save_image(p = fig, file = file_path, height = 800, width = 800)
  })
  
  
})

## plot 2D ----
library(PieGlyph)
data_bg <- filter(data_amp, site %in% site_list, region %in% c("LM")) %>%
  filter(max(vaf, na.rm = T) <= vaf_thred, .by = sample) %>% 
  select(region, tumor, sample, x, y, z) %>% 
  distinct() %>%
  mutate(vaf = 0, site = "bg") %>% 
  na.omit() 

data_plot <- filter(data_amp, 
                    site %in% site_list,
                    region %in% c("LM")) %>% 
  select(-c(temp, total, alt, ref)) %>% 
  filter(vaf > vaf_thred) %>% 
  bind_rows(data_bg) %>% 
  na.omit() %>% 
  mutate(z = factor(z, 
                    levels= unique(.$z) %>% sort(decreasing = T),
                    labels = paste0("layer", sprintf("%02d", 1:length(unique(.$z))))),
         vaf = 1)   %>% 
  filter(!z %in% c("layer01", "layer14", "layer15")) %>% 
  pivot_wider(names_from = site, values_from = vaf, values_fill = 0)

filter(data_plot, duplicated(sample))

site_list <- c("chr11_56408815","chr18_24039766", "chr20_31669258", "chr4_158094116")
color_list <- c("steelblue", "darkred",  "darkgreen", "#FFA500")
names(color_list) <- site_list
color_list <- c(color_list, bg = "grey")

plot_case2_2d <- ggplot(data = data_plot, 
                        aes(x=x, y = y)) +
  geom_pie_glyph(slices = site_list) +
  geom_circle(aes(x0 = 3, y0 = 3, r = 4.2),
              linetype =2)+
  facet_wrap(~z) +
  theme_bw() +
  scale_fill_manual(values = color_list) +
  ylim(c(7.2, -1.2)) +
  coord_fixed(ratio = 1) +
  theme(axis.text = element_text(size = 16),
        strip.text = element_text(size = 21),
  ) 
plot_case2_2d
