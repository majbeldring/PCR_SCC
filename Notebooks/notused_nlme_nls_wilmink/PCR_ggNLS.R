

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Wilmink curves with nls for: PCR TEST POS OR NEG

# Script includes:
## 1) preparedata: df_pos & df_neg for each parity -> 6 datasets in total to run
## 2) create 3 curves (1 with two curves for each parity) using boot (not herd level)
## 3) create curves on herd level
## 4) ggpairs


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

load("M:/PCR_data/PCR_merge2.RData") 
rm(df_pcr, df_curve); gc() 

# pre-preparing data;remove dates and DYR_ID:
df_model <- df_model %>% 
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)


# descriptives on over all data:
dplyr::n_distinct(df_model$BES_ID) # 3474 herds
dplyr::n_distinct(df_model$DYR_ID) # 997.784


summary(df_model)

#---------------------------------------------------------------------
# choose analysis focus: here PCR

# Choose: 1) Test type 2) breed, 3) parity, 4) herdtype:
df <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 4) %>%
  filter(HERDTYPE == 1) %>%
  filter(DIM < 306) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, 
                MILK, RES_MAJOR)

dplyr::n_distinct(df_pos$BES_ID) # 1141
dplyr::n_distinct(df$DYR_ID) # 3075

#----------------------------------------------------------------------
# specify data

### RES MAJOR POS : only treated with max 200 observations:
df_pos <- df %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
df_pos$BES_ID <- factor(df_pos$BES_ID) # keep only used levels by resetting the variable


### RES MAJOR NEG : only treated with max 200 observations:
df_neg <- df %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
df_pos$BES_ID <- factor(df_pos$BES_ID) # keep only used levels by resetting the variable




#------------------------------------------------------


# nls in ggplot - ONLY with test data when facet_wrap is applied
p_pos <- ggplot(df_pos, aes(x = DIM, y = logSCC)) + 
  ggtitle("nls, logSCC~DIM") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_pos

# nls in ggplot:
p_neg <- ggplot(df_neg, aes(x = DIM, y = logSCC)) + 
  ggtitle("parity 3+ neg, logSCC~DIM") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_neg







