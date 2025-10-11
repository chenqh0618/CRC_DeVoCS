
na_filter <- function(data, type, n_thod = 0.8 ){
  temp <- as.data.frame(data) %>% 
    rownames_to_column("row") %>%
    pivot_longer(cols = -row, 
                 names_to = "col", 
                 values_to = "value")
  
  if(type == "row"){
    ncol = ncol(data)
    nna_idx <- filter(temp, !is.na(value)) %>% 
      group_by(row) %>% 
      count() %>% 
      mutate(n = n/ncol) %>% 
      filter(n > n_thod) %>% 
      pull(row)
    
    res <- filter(temp, row %in% nna_idx)
  }
  
  if(type == "col"){
    nrow = nrow(data)
    nna_idx <- filter(temp, !is.na(value)) %>% 
      group_by(col) %>% 
      count() %>% 
      mutate(n = n/nrow) %>% 
      filter(n > n_thod) %>% 
      pull(col)
    
    res <- filter(temp, 
                  col %in% nna_idx)
    
  }
  
  res <- pivot_wider(res, names_from = col, values_from = value) %>% 
    column_to_rownames("row")
  
  paste0(dim(res)[1], " rows and ", dim(res)[2], " cols is remain!") %>% 
    print()
  
  return(res)
}

