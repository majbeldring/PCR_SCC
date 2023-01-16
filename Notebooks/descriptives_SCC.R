

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Descriptives on data to paper

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
df_count <- df_model %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, SCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB)

# DIM: day 6 to day 305 (already cleaned in data_SCC)

# descriptives on over all data:
dplyr::n_distinct(df_count$BES_ID) # 3474 herds
dplyr::n_distinct(df_count$DYR_ID) # 997.406

summary(df_count)
# SCC: mean = 254, median = 85

# Holstein, conventionel and PCR tested in Parity 2,3 and 3+ (in merge2 parities have been selected):
df_count <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1)

dplyr::n_distinct(df_count$BES_ID) # 1368
dplyr::n_distinct(df_count$DYR_ID) # 133934

summary(df_count)
# DIM: mean= 152.7 days, median = 152 days
# SCC: mean = 285, median = 93


#-------------------------------------------------------------
### PARITY 2: 
df2 <- df_count %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 2) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, SCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df2_pos <- df2 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC, DYR_ID, SCC, MILK)
df2_pos$BES_ID <- factor(df2_pos$BES_ID) 

### RES MAJOR NEG
df2_neg <- df2 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC, DYR_ID, SCC, MILK)
df2_neg$BES_ID <- factor(df2_neg$BES_ID) 

### PARITY 3: 
df3 <- df_count %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 3) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, SCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df3_pos <- df3 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC, DYR_ID, SCC, MILK)
df3_pos$BES_ID <- factor(df3_pos$BES_ID) 

### RES MAJOR NEG
df3_neg <- df3 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC, DYR_ID, SCC, MILK)
df3_neg$BES_ID <- factor(df3_neg$BES_ID) 

### PARITY 4: 
df4 <- df_count %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 4) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, SCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df4_pos <- df4 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC, DYR_ID, SCC, MILK)
df4_pos$BES_ID <- factor(df4_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df4_neg <- df4 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC, DYR_ID, SCC, MILK)
df4_neg$BES_ID <- factor(df4_neg$BES_ID) # keep only used levels by resetting the variable


## COUNTING FOR:
## COUNTS >200, POS OR NEG, PARITY DIVIDED
# Herds:
dplyr::n_distinct(df2_pos$BES_ID) # 551
dplyr::n_distinct(df2_neg$BES_ID) # 575
dplyr::n_distinct(df3_pos$BES_ID) # 429
dplyr::n_distinct(df3_neg$BES_ID) # 396
dplyr::n_distinct(df4_pos$BES_ID) # 327
dplyr::n_distinct(df4_neg$BES_ID) # 283

# animals:
dplyr::n_distinct(df2_pos$DYR_ID) # 41664
dplyr::n_distinct(df2_neg$DYR_ID) # 44701
dplyr::n_distinct(df3_pos$DYR_ID) # 23051
dplyr::n_distinct(df3_neg$DYR_ID) # 40580
dplyr::n_distinct(df4_pos$DYR_ID) # 12549
dplyr::n_distinct(df4_neg$DYR_ID) # 10728




