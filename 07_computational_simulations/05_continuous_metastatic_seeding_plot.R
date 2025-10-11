
pjdir <- "/public/data/cqh_project/crc/amp_case1/"
setwd(pjdir)
library(tidyverse)
library(data.table)
library(parallel)

source("~/software/script_tools/R_color.R")
source("~/software/script_tools/function_lib.R")

res_dir <- "res_plot//inter_tumor_migration_plot"
dir.create(res_dir, recursive = T)

# function the migration process -----------------------------------------------

migration <- function(id, mig_r, grow_r, t_ini, t_curt){
  t = t_curt - t_ini
  if(t <= 0) {return(NULL) }
  mother_id = id
  mother_mig_r = mig_r
  mother_grow_r = grow_r
  
  Nt = K / (1.0 + (K - 1.0) * exp(-mother_grow_r * t)) 
  Nt = min(1e9, Nt)
  
  Nm = rbinom(n = 1, size = as.integer(Nt^(2/3)), prob = mother_mig_r)
  # Nm = rbinom(1, as.integer(Nt), mother_mig_r*Nt^(-1/3))
  
  temp <- NULL
  if(Nm > 0){
    temp <- data.frame(id = paste0(mother_id, "__t", t_curt,"_mig", 1:Nm),
                       t_ini = t_curt,
                       grow_r = rexp(n = Nm, 
                                     rate = 1/mother_grow_r),
                       mig_r = mother_mig_r) 
  }
  return(temp)
}


#parameter setting -------------------------------------------------------------
grow_r = 5e-3
mig_r = 1e-6
rep = 1

keywd <- paste0( "grow_",grow_r, "_mig_", mig_r, "_rep" ,rep)



mig_r1 = mig_r2 = mig_r
K = 2e9
N_surgery = 1e9

Tmax <- -1.0 / grow_r * log((K - N_surgery) / (N_surgery * K - N_surgery))  
Tmax <- ceiling(Tmax) 

# simulation the migration process ---------------------------------------------
primary <- data.frame(id = "primary",
                      t_ini = 0,
                      grow_r = grow_r,
                      mig_r = mig_r1) %>% 
  mutate(wave = 0)

first_wave <- mclapply(1:Tmax, mc.cores = 12,
                       migration,
                       id = "primary", 
                       mig_r = mig_r1,
                       grow_r =  grow_r,
                       t_ini = 1 ) %>% 
  do.call(what = rbind) %>% 
  filter(grow_r < 0.05)  %>%
  mutate(wave = 1)


second_wave <- mclapply(1:Tmax, mc.cores = 12, function(tc){
  temp <- as.data.frame(first_wave) %>% 
    filter(tc - t_ini > 0)
  if(nrow(temp) == 0) return(NULL)
  
  lapply(1:nrow(temp), function(i){
    migration(id = temp$id[i],
              t_ini = temp$t_ini[i],
              mig_r = mig_r2,
              grow_r = temp$grow_r[i],
              t_curt = tc)
  }) %>% 
    do.call(what = rbind)
})  %>% 
  do.call(what = rbind) %>% 
  filter(grow_r < 0.05)  %>%
  mutate(wave = 2)


res <- rbind(primary, first_wave, second_wave)
write_rds(res, paste0(res_dir, "/simu_",keywd, ".rds"))


# plot the process -------------------------------------------------------------
temp <- res; dim(temp)
Tmax_plot <- 5000
data_plot <- mclapply(1:nrow(temp), mc.cores = 12, function(i){
  id = temp$id[i]
  t_ini = temp$t_ini[i]
  t_list = seq(t_ini, Tmax_plot, length = min((Tmax_plot-t_ini +1), 200))
  grow_r = temp$grow_r[i]
  wave = str_count(id, "mig")

  data.frame(id = id,
             wave = wave,
             t = t_list,
             Nt =  K / (1.0 + (K - 1.0) * exp(-grow_r * (t_list - t_ini))))
}) %>% do.call(what = rbind)



p1 <- ggplot(data = data_plot, aes(x = t, y = Nt, group = id)) +
  geom_line(data = filter(data_plot, wave == 0), color = "darkgreen", linewidth = 2) +
  # geom_line(data = filter(data_plot, wave == 2), color = "darkred", linewidth = 0.5) +
  geom_line(data = filter(data_plot, wave == 1), color = "pink", linewidth = 1) +
  geom_line(data = filter(data_plot, wave == 2), color = "darkred", linewidth = 0.5) +
  geom_hline(yintercept = 1e7, linetype = 2, linewidth = 2, color = "steelblue") +
  geom_hline(yintercept = 1e9, linetype = 2, linewidth = 2, color = "darkred") +
  scale_y_log10() +
  theme_classic() +
  theme(axis.text = element_text(size = 16))
p1

ggsave(plot = p1, filename = paste0(res_dir, "/plot_process_", keywd, ".pdf"),
       width = 16, height = 12)

