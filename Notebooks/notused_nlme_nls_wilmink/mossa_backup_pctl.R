#---------------------------------------------------
# Percentiles:


# df2_neg / out_neg2


out_neg2 %>% class()
neg2_parameters <- nlme_out_neg2 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))
neg2_parameters_pct <- neg2_parameters  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) 
  # %>% mutate(pctile = names(a)) %>% select(pctile, everything())

neg2_parameters <- nlme_out_neg2 %>% 
  summarise(across(everything(), mean))

neg2_parameters = list(mean = neg2_parameters, 
                       pctile = neg2_parameters_pct) %>% 
  enframe("metode", "value") %>% 
  unnest(value)

neg2_parameters_pct %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      aes(group = pctile) +
      aes(color = pctile) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()
