

# Maj Beldring Henningsen, majbh@sund.ku.dk

# wilmink for SCC curves
# nlme

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
library(lme4) # for nlmer
#library(ggExtra)
#library(ggalluvial)
Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

load("M:/PCR_data/curves_SCC.RData") # from PCR_data script


#-----------------------------------------------------
# create nls parameters:
library(nls.multstart)

f_wilmink <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}

nls_neg4 <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                            data = df4_neg,
                                            lower=c(a=0, b=0, k=-5, d=0),
                                            upper=c(a=9, b=1.5, k=0, d=5),
                                            start_lower = c(a=0, b=0, k=-5, d=0),
                                            start_upper = c(a=8, b=1, k=-0.01, d=4),
                                            iter = 500,
                                            supp_errors = "Y")

coef(nls_neg4)

nls_pos4 <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                         data = df4_pos,
                                         lower=c(a=0, b=0, k=-5, d=0),
                                         upper=c(a=9, b=1.5, k=0, d=5),
                                         start_lower = c(a=0, b=0, k=-5, d=0),
                                         start_upper = c(a=8, b=1, k=-0.01, d=4),
                                         iter = 500,
                                         supp_errors = "Y")

coef(nls_pos4)

#-----------------------------------------------------
# test nlme wilmink:
# from: https://quantdev.ssri.psu.edu/sites/qdev/files/GMChp10_Tutorial.html

# treating like latent parameters

nlme_neg4 <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                  data=df4_neg,
                  fixed=a+b+k+d~1,
                  random=a+b+k+d~1,
                  groups=~BES_ID,
                  start=c(a = 3.9, b = 0.003, k = -1.92, d = 2.5),
                  na.action=na.exclude,
                  control=(msMaxIter=400))

out_nlme_neg4 <- coef(nlme_neg4)

# running with many warnings! Running without DIM/150
test_nlme <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                          data=df4_neg,
                          fixed=a+b+k+d~1,
                          random=a+b+k+d~1,
                          groups=~BES_ID,
                          start=c(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
                          na.action=na.exclude)

summary(test_nlme)
out_test_nlme <- coef(test_nlme)






# pos testing
# nls out: 
# not running with: a = 4.2, b = 0.003, k = -1.8, d = 3.0
test_nlme_pos <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                  data=df4_pos,
                  fixed=a+b+k+d~1,
                  random=a+b+k+d~1,
                  groups=~BES_ID,
                  start=c(a = 4.2, b = 0.003, k = -1.8, d = 3.0),
                  na.action=na.exclude)

summary(test2_nlme)
out_nlme4pos <- coef(test2_nlme)

# compare outliers with nls:
library(rstatix)
out_neg4 %>%
  identify_outliers("d")
out_nlme4neg %>%
  identify_outliers("d")



# compare with nls output by percentage:

out_neg4 %>% 
  dplyr::select(-BES_ID) %>%
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))

out_nlme4neg %>% 
  #dplyr::select(-BES_ID) %>%
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))


#---------------------------------------------------------------
# plot curve mean and percentiles:

# define wilmink eq 8:
logSCC_func <- function(DIM, a,b,k,d) 
  a + b * DIM + exp(-(exp(k)) * DIM)*d 

# nlme mean:PARITY 4 - NEG:




## POS 4:

pos4_mean <- out_nlme4pos %>% 
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



#-----------------------------------------------------------------
# identify on curve:

# minimum:


