

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
library(rstatix) # outliers
#library(lme4) # for nlmer
#library(ggExtra)
#library(ggalluvial)
Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

load("M:/PCR_data/PCR_merge.RData") 
rm(df_pcr, df_curve); gc() 

# pre-preparing data;remove dates and DYR_ID:
df <- df_model %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, PCR_TEST, RES_MAJOR)

df <- df %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 4) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, RES_MAJOR)

df$BES_ID <- factor(df$BES_ID) 
df <- df %>% filter(BES_ID == 5051112)
df$BES_ID <- factor(df$BES_ID) # removed unused levels by resetting

df_pos <- df %>%
  filter(RES_MAJOR == 1) %>%
  #group_by(BES_ID, PARITY) %>%
  #mutate(count = n()) %>%
  #filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC)

df_neg <- df %>%
  filter(RES_MAJOR == 0) %>%
  #group_by(BES_ID, PARITY) %>%
  #mutate(DIM = DIM/150) %>% # making it faster to run model
  #mutate(count = n()) %>%
  #filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC)


# herd 5051112
dplyr::n_distinct(df_pos$DYR_ID) # 160
dplyr::n_distinct(df_neg$DYR_ID) # 257


#-------------------------------------------------------



#----------------------------------------------------------
hill_nls_outcomes_summary <- hill_nls_outcomes %>%
  group_by(term) %>% 
  summarise(mean = mean(estimate), 
            median = median(estimate),
            sd = sd(estimate)) %>% 
  ungroup()

knitr::kable(hill_nls_outcomes_summary, digits = 3)




#-----------------------------------------------------
# test nlme wilmink:
# from: https://quantdev.ssri.psu.edu/sites/qdev/files/GMChp10_Tutorial.html


upper_bound <- c(9,1.5,0,5)
lower_bound <- c(0,0,-5,0)


nlme_neg <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                 data=df_neg,
                 fixed=a+b+k+d~1,
                 random=a+b+k+d~1, # not interest in why DYR is equal - can leave out residuals
                 groups=~DYR_ID,
                 start=c(a = 4.02, b = 0.0035, k = -1.88, d = 2.9),
                 na.action=na.exclude,
                 control=(nlmeControl(opt = "nlminb", 
                                      optCtrl = list(method = "nlminb",
                                                     maxit=10000,
                                                     iter.max=10000,  
                                                     eval.max=10000,
                                                     lower = lower_bound, 
                                                     upper = upper_bound))))



nlme_pos <- nlme(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150)*d,
                          data=df_pos,
                          fixed=a+b+k+d~1,
                          random=a+b+k+d~1, # not interest in why DYR is equal - can leave out residuals
                          groups=~DYR_ID,
                          start=c(a = 4.02, b = 0.0035, k = -1.88, d = 2.9),
                          na.action=na.exclude,
                          control=(msMaxIter=200))

nlme_neg <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                 data=df_neg,
                 fixed=a+b+k+d~1,
                 #random=a+b+k+d~1,
                 groups=~DYR_ID,
                 start=c(a = 3.88, b = 0.0030, k = -1.80, d = 2.7),
                 na.action=na.exclude,
                 control=(msMaxIter=1000)
                 )

f_pos <- nlme(logSCC ~ SSasymp(age, Asym, R0, lrc),
            data = df_pos,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))

hist(df_pos$DIM)

summary(test_nlme)
out_nlme4pos <- coef(test2_nlme)


#-------------------------------------------------------------------
# output:

# outliers:
out_nlme4neg %>% identify_outliers("d")


# percentage:
out_nlme4neg %>% 
  #dplyr::select(-BES_ID) %>%
  as_tibble() %>% 
  # summarise(across(everything(), mean))
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))


#--------------------------------------------------------------
# nls test

# STARTLIST defined and equation 8
st <- list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6)
f_nls <- logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | BES_ID

# test with upper and lower boundries:
upper <- c(8,1,0,4)
lower <- c(0,0,-4,0)

control = nlmeControl(opt = "nlminb", upper = upper_bounds)

nls_1 <- nlsList(f_nls, df4_neg, start = sapply(st, mean),
                  control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                 printEval = FALSE, warnOnly = TRUE), upper = upper, lower = lower)

out <- coef(nls_1)
out <- cbind(BES_ID = rownames(out), out)



#---------------------------------------------------------------
# plot curve mean and percentiles:

# define wilmink eq 8:
logSCC_func <- function(DIM, a,b,k,d) 
  a + b * DIM + exp(-(exp(k)) * DIM)*d 

# extract mean:
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


