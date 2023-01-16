


# Packages and settings ---------------------------------------------


library(tidyverse)
library(ggpubr) # p values in boxplot
library(gridExtra) # gridarrange
ggplot2::theme_set(ggplot2::theme_bw())  # globally sets ggplot2 theme to theme_bw
library(GGally) # for ggpairs
library(nls.multstart)
library(nlme) # for nlslist

Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


# Loading data and preparing data -------------------------------------------



# load("M:/PCR_data/curve6_SCC.RData") # from PCR crypted container
load("R:/PCR_data/curve6_nlme_wilmink.RData")
# save.image("M:/PCR_data/curve6_nlme_wilmink.RData")


# create nls parameters:----------------------------------------------


f_wilmink <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}


# nls multistart (not grouped herd level):
# repeat for all 6 diff. datasets
nls_pos4 <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                         data = df4_pos,
                                         lower=c(a=0, b=0, k=-5, d=0),
                                         upper=c(a=9, b=1.5, k=0, d=5),
                                         start_lower = c(a=0, b=0, k=-5, d=0),
                                         start_upper = c(a=8, b=1, k=-0.01, d=4),
                                         iter = 500,
                                         supp_errors = "Y")


coef(nls_pos4)


# nlme wilmink  ---------------------------------------


# repeat for all 6 different data sets
nlme_neg2 <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                  data=df2_neg,
                  fixed=a+b+k+d~1,
                  random=a+b+k+d~1,
                  groups=~BES_ID,
                  start=c(a = 3.5, b = 0.0032, k = -2.3, d = 1.9),
                  na.action=na.exclude,
                  control = list(maxIter = 50, msMaxIter = 50))
#'
#'
#'
good_herds <- df2_neg %>% 
  count(BES_ID) %>% 
  slice_max(n, n = 10)
good_herds

good_herds <- df2_neg %>% 
  semi_join(good_herds, by = "BES_ID")

good_herds %>% 
  naniar::miss_var_summary()

good_herds %>% summarise(range(logSCC))

good_herds %>% 
  mutate(SCC = exp(logSCC)) %>% 
  mutate(loglogSCC = log(logSCC)) %>% 
  
  pivot_longer(c(SCC, logSCC, loglogSCC)) %>%  
  mutate(name = name %>% fct_inorder()) %>% 
  
  identity() %>% {
    ggplot(.) + 
      aes(value, group = interaction(BES_ID, name)) + 
      
      geom_freqpoly(aes(y=  after_stat(density))) +
      # stat_bin(aes(y = after_stat(density)), geom = "step") + 
      # stat_ecdf() + 
      facet_wrap(~name, scales = "free") + 
      NULL
  }
#'
#'
df_model %>% 
  rename_with(tolower) %>% 
  ungroup() %>% 
  glimpse() %>% 
  filter(parity == 2, res_major == 0) %>% 
  semi_join(
    count(., bes_id) %>% slice_max(n, n = 10),
    by = "bes_id"
  ) -> nice_cows
#'
#'
nice_cows %>%
  rename(logSCC = logscc) %>% 
  mutate(SCC = exp(logSCC)) %>% 
  mutate(loglogSCC = log(logSCC)) -> 
  nice_cows_scc

nice_cows_scc %>% 
  select(everything(), SCC, logSCC, loglogSCC) %>% 
  # nest(bes_id, dyr_id
  select(bes_id, dyr_id, parity, dim,  SCC, logSCC, loglogSCC) %>% 
  group_by(bes_id, dyr_id, parity) %>% 
  nest(data = c(dim, SCC, logSCC, loglogSCC)) %>% 
  ungroup() %>% 
  #' some cows are observed 43 times in one lactation...
  # mutate(length_n = data %>% map_dbl(nrow)) %>%
  # slice_min(length_n, n = 10)
  filter(data %>% map_dbl(nrow) %>% {. > 1}) %>% 
  mutate(shapiro.test = data %>% map(~ summarise(.x, 
                                                 across(c(SCC, logSCC, loglogSCC), possibly(. %>% shapiro.test %>% `[[`("p.value"), 
                                                                                            otherwise = NA))))) ->
  nice_cows_scc_normality

nice_cows_scc_normality %>% 
  unnest_wider(shapiro.test) ->
  nice_cows_scc_normality
View(nice_cows_scc_normality)

nice_cows_scc %>% 
  select(-data) %>%
  unnest_wider(shapiro.test) %>% 
  unnest(data) %>% 
  identity()

nice_cows_scc %>% 
  # sample_n(size = 10)
  pivot_longer(c(SCC, logSCC, loglogSCC)) %>%  
  mutate(name = name %>% fct_inorder()) %>% 
  
  identity() %>% {
    ggplot(.) + 
      aes(value, group = interaction(bes_id, dyr_id, name)) + 
      
      # geom_freqpoly(aes(y=  after_stat(density))) +
      # stat_bin(aes(y = after_stat(density)), geom = "step") + 
      stat_ecdf() +
      facet_wrap(~name, scales = "free") + 
      NULL
  }
