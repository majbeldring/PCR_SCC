

# Maj Beldring Henningsen, majbh@sund.ku.dk

# preparing data for Wilmink SCC curves
# 6 pathogens

# Packages and settings:---------------------------------------


library(tidyverse)
library(gridExtra)
library(data.table)

Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


# Loading data and preparing data ---------------------------------------


load("M:/PCR_data/curve6_merge.RData") 
rm(df_pcr, df_curve); gc() 

# pre-preparing data;remove dates and DYR_ID:
df_model <- df_model %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  # SCC added for descriptive statistics (not included in saved data)
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, SCC, logSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)


# descriptives on over all data:
dplyr::n_distinct(df_model$BES_ID) # 3474 herds
dplyr::n_distinct(df_model$DYR_ID) # 997.406


#summary(df_model)


# Divide into Parity groups, conventionel herds, Holstein -----------------------

### PARITY 2: 
df2 <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 2) %>%
  filter(HERDTYPE == 1) %>%
  # added SCC for descpritives
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, SCC, logSCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df2_pos <- df2 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df2_pos$BES_ID <- factor(df2_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df2_neg <- df2 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df2_neg$BES_ID <- factor(df2_neg$BES_ID) # keep only used levels by resetting the variable



### PARITY 3: 
df3 <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 3) %>%
  filter(HERDTYPE == 1) %>%
  # added SCC for descpritives
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, SCC, logSCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df3_pos <- df3 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df3_pos$BES_ID <- factor(df3_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df3_neg <- df3 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df3_neg$BES_ID <- factor(df3_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df3_pos$BES_ID) # 451
dplyr::n_distinct(df3_neg$BES_ID) # 343



### PARITY 4: 
df4 <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 4) %>%
  filter(HERDTYPE == 1) %>%
  # added SCC for descriptive
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, SCC, logSCC, 
                MILK, RES_MAJOR)

### RES MAJOR POS ; min 200 obs
df4_pos <- df4 %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df4_pos$BES_ID <- factor(df4_pos$BES_ID) # keep only used levels by resetting the variable

### RES MAJOR NEG
df4_neg <- df4 %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df4_neg$BES_ID <- factor(df4_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df4_pos$BES_ID) # 329
dplyr::n_distinct(df4_neg$BES_ID) # 220



save.image("M:/PCR_data/curve6_SCC.RData")
