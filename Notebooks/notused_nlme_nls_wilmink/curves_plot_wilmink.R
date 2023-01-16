

# Maj Beldring Henningsen, majbh@sund.ku.dk

# curves created on herd level using mean, median or pecentiles

#-------------------------------------------------------
# Packages and settings:

library(tidyverse)
library(gridExtra)
library(data.table)
library(plotly)
library(GGally)
library(tidymodels)
library(nlstools) # for bootstrapping
library(nlme) # for nlslist
#library(ggExtra)
#library(ggalluvial)
Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

load("M:/PCR_data/wilmink_SCC_out.RData") # from NLSlist_running_crazy_outliers

#-------------------------------------------------------------------------------------------
# plot on herd level

# define wilmink eq 8:
logSCC_func <- function(DIM, a,b,k,d) 
  a + b * DIM + exp(-(exp(k)) * DIM)*d 

#------------------------------------------------------------------------------------------
# Percentile plot:

# Note, must remove BES_ID before possibe

### PARITY 2 POS and NEG
# POS
pos2_pct <- out_pos2 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

pos2_pct <- pos2_pct  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())

pos2_pct %>% 
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
      ggtitle("pctile of parametres - SCC Parity 2: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()

# NEG:
neg2_pct <- out_neg2 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

neg2_pct <- neg2_pct  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())

neg2_pct %>% 
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
      ggtitle("pctile of parametres - SCC Parity 2: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()




### PARITY 3: POS and NEG
# POS
pos3_pct <- out_pos3 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

pos3_pct <- pos3_pct  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())

pos3_pct %>% 
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
      ggtitle("pctile of parametres - SCC Parity 3: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()

# NEG:
neg3_pct <- out_neg3 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

neg3_pct <- neg3_pct  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())

neg3_pct %>% 
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
      ggtitle("pctile of parametres - SCC Parity 3: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()



### PARITY 4: POS and NEG
# POS
pos4_pct <- out_pos4 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

pos4_pct <- pos4_pct  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())

neg4_pct %>% 
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
      ggtitle("pctile of parametres - SCC Parity 4: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()

# NEG:
neg4_pct <- out_neg4 %>% 
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

neg4_pct <- neg4_pct  %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())

neg4_pct %>% 
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
      ggtitle("pctile of parametres - SCC Parity 4: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()



#-----------------------------------------------------------------------------------
# Plot with median:


#### PARITY 2:
# PARITY 2 - POS:
pos2_median <- out_pos2 %>% 
  summarise(across(everything(), median))

pos2_median %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 2: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()


# PARITY 2 - NEG:
neg2_median <- out_neg2 %>% 
  summarise(across(everything(), median))

neg2_median %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 2: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()




#### PARITY 3:
# PARITY 3 - POS:
pos3_median <- out_pos3 %>% 
  summarise(across(everything(), median))

pos3_median %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 3: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()


# PARITY 3 - NEG:
neg3_median <- out_neg3 %>% 
  summarise(across(everything(), median))

neg3_median %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 3: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()



#### PARITY 4:
# PARITY 4 - POS:
pos4_median <- out_pos4 %>% 
  summarise(across(everything(), median))

pos4_median %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 4: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()


# PARITY 4 - NEG:
neg4_median <- out_neg4 %>% 
  summarise(across(everything(), median))

neg4_median %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 4: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()







#-----------------------------------------------------------------------------------
# Plot with MEAN
# (doesn't make sense until outliers have been removed)


#### PARITY 2:
# PARITY 2 - POS:
pos2_mean <- out_pos2 %>% 
  summarise(across(everything(), mean))

pos2_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("mean of parametres - SCC Parity 2: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()


# PARITY 2 - NEG:
neg2_mean <- out_neg2 %>% 
  summarise(across(everything(), mean))

neg2_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("mean of parametres - SCC Parity 2: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()




#### PARITY 3:
# PARITY 3 - POS:
pos3_mean <- out_pos3 %>% 
  summarise(across(everything(), mean))

pos3_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("median of parametres - SCC Parity 3: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()


# PARITY 3 - NEG:
neg3_mean <- out_neg3 %>% 
  summarise(across(everything(), mean))

neg3_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("mean of parametres - SCC Parity 3: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()



#### PARITY 4:
# PARITY 4 - POS:
pos4_mean <- out_pos4 %>% 
  summarise(across(everything(), mean))

pos4_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("mean of parametres - SCC Parity 4: PCR POS") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()


# PARITY 4 - NEG:
neg4_mean <- out_neg4 %>% 
  summarise(across(everything(), mean))

neg4_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), logSCC_func)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("mean of parametres - SCC Parity 4: PCR NEG") +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()

