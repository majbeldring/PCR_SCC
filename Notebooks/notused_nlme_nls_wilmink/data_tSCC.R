

# Maj Beldring Henningsen, majbh@sund.ku.dk

# preparing data for Wilmink SCC curves
# note: data are called the same as for the MILK; 
# so don't run SCC and milk simultaneous


#-------------------------------------------------------
# Packages and settings:

library(tidyverse)
library(gridExtra)
library(data.table)
#library(plotly)
#library(GGally)
#library(tidymodels)
#library(nlstools) # for bootstrapping
#library(nlme) # for nlslist
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
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  mutate(tSCC = SCC * MILK) %>%
  mutate(logtSCC = log(tSCC)) %>%
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, logtSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)


# descriptives on over all data:
dplyr::n_distinct(df_model$BES_ID) # 3474 herds
dplyr::n_distinct(df_model$DYR_ID) # 997.784


summary(df_model)

#---------------------------------------------------------------------
# Specify data:
# conventionel, holstein, DIM < 306, Must have PCR test

# divided into 3 parity groups

#-------------------------------------------------------------
### PARITY 2: 
df2 <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 2) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, logtSCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df2_pos <- df2 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logtSCC)
df2_pos$BES_ID <- factor(df2_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df2_neg <- df2 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logtSCC)
df2_neg$BES_ID <- factor(df2_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df2_pos$BES_ID) # 1141
summary(df2_pos)


#-------------------------------------------------------------
### PARITY 3: 
df3 <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 3) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, logtSCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df3_pos <- df3 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logtSCC)
df3_pos$BES_ID <- factor(df3_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df3_neg <- df3 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logtSCC)
df3_neg$BES_ID <- factor(df3_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df3_pos$BES_ID) # 416

#-------------------------------------------------------------
### PARITY 4: 
df4 <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 4) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, logtSCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df4_pos <- df4 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logtSCC)
df4_pos$BES_ID <- factor(df4_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df4_neg <- df4 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logtSCC)
df4_neg$BES_ID <- factor(df4_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df4_pos$BES_ID) # 303



# save.image("M:/PCR_data/curves_tSCC.RData") # not saved
