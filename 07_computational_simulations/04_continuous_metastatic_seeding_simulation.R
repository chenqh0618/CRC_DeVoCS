setwd("~/project/crc/simulation/")
.libPaths("/home/data/cqh/R/x86_64-pc-linux-gnu-library/4.4")
library(data.table)
library(tidyverse)
library(parallel)

res_dir <- "res_simu//res_differ_mig_grow_v1d"
dir.create(res_dir, recursive = T)

args = commandArgs(trailingOnly=TRUE)
# grow_r = 1e-3
# mig_r = 1e-6
# rep = 1
grow_r = args[1] %>% as.numeric() 
mig_r = args[2] %>% as.numeric()
rep = args[3]

keywd <- paste0( "grow_",grow_r, "_mig_", mig_r, "_rep" ,rep)




mig_r1 = mig_r2 = mig_r
K = 2e9
N_surgery = 1e9

Tmax <- -1.0 / grow_r * log((K - N_surgery) / (N_surgery * K - N_surgery))  
Tmax <- ceiling(Tmax) 


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


# first waved results ----------------------------------------------------------
primary <- data.frame(id = "primary",
                      t_ini = 0,
                      grow_r = grow_r,
                      mig_r = mig_r1,
                      wave = 0)



first_wave <- mclapply(1:Tmax, mc.cores = 8, function(tc){
  migration(id = "primary", 
            mig_r = mig_r1,
            grow_r =  grow_r,
            t_ini = 1, 
            t_curt = tc)
}) %>% 
  do.call(what = rbind) 

if (is.null(first_wave)){
  res_firstwave <- rbind(primary, first_wave) 
  
  write_rds(res_firstwave, file = paste0(res_dir,"/firstMig_",  keywd,  ".rds"))
  
  
  res_summary <- mutate(res_firstwave, 
                        Nsize = K / (1.0 + (K - 1.0) * exp(- grow_r * (Tmax-t_ini))),
                        label = if_else(Nsize >= 1e7, "oberseved", "unobserved")  ) %>% 
    group_by(wave, label) %>% 
    reframe(n = length(t_ini), 
            time = min(t_ini))
  write_rds(res_summary, file = paste0(res_dir, "/resSummary_", keywd, ".rds" ))
  
  stop("no migration happened!")
}


first_wave <- as.data.frame(first_wave) %>% 
  mutate(wave = 1) %>% 
  filter(grow_r < 0.05)   
res_firstwave <- rbind(primary, first_wave) 
write_rds(res_firstwave, file = paste0(res_dir,"/firstMig_",  keywd,  ".rds"))


res_summary <- mutate(res_firstwave, Nsize = K / (1.0 + (K - 1.0) * exp(- grow_r * (Tmax-t_ini))),
                      label = if_else(Nsize >= 1e7, "oberseved", "unobserved")  ) %>% 
  group_by(wave, label) %>% 
  reframe(n = length(t_ini), 
          time = min(t_ini))



# second wave results ----------------------------------------------------------
res_second_wave <- lapply(
  min(first_wave$t_ini+1):Tmax,
  # mc.cores = 8, 
  function(tc){
    seed_cell <- filter(first_wave, t_ini < tc)
    
    temp_res <- mclapply(1:nrow(seed_cell), mc.cores = 32, function(i){
      temp <- migration(id = seed_cell$id[i],
                        t_ini = seed_cell$t_ini[i],
                        mig_r = mig_r2,
                        grow_r = seed_cell$grow_r[i],
                        t_curt = tc)  
      
      if(is.null(temp))  return(NULL) 
      temp <- filter(temp, grow_r < 0.05)
      if(nrow(temp) == 0 )  return(NULL)
      
      temp_res_summary <- mutate(temp,
                                 Nsize = K / (1.0 + (K - 1.0) * exp(- grow_r * (Tmax-t_ini))),
                                 label = if_else(Nsize >= 1e7, "oberseved", "unobserved")  ) %>%
        group_by(label) %>% 
        reframe(n = length(t_ini), 
                time = min(t_ini))
    }) %>%
      do.call(what = rbind)  
    
    if(is.null(temp_res))  return(NULL) 
    if(nrow(temp_res) == 0 )  return(NULL)
    temp_res <- mutate(temp_res, wave = 2, time=tc)
    
    return(temp_res)
  }) %>% 
  do.call(what = rbind)

if(!is.null(res_second_wave)){
  res_second_wave <- group_by(res_second_wave, label) %>% 
    reframe(n = sum(n),
            time = min(time)) %>% 
    mutate(wave = 2)
  
  res_summary <- rbind(res_summary, res_second_wave) %>% 
    mutate(mig_r = mig_r,
           grow_r = grow_r,
           rep = rep)
  } 



write_rds(res_summary, 
          file = paste0(res_dir, "/resSummary_", keywd, ".rds" ))














