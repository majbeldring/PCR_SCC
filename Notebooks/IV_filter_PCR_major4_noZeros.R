#' 
#' 
#' Maj Beldring, majbh@sund.ku.dk
#' UCPH, November 2021
#' 
#' Paper 1, PCR, wilmink: 
#' Filter data for Wilmink modelling in regards to PCR result
#' Script #4 in PCR project
#'
#' redo to create df1_neg and df1_pos also
#'
# Packages and settings: ----------------------------------------


library(tidyverse)
library(gridExtra)
library(data.table)

Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


# Loading data and preparing data ---------------------------

load("K:/paperI/major4/III_merge_PCR_major4.RData") 
rm(df_lac, major, minor, pcr_full); gc() 


load("K:/paperI/major4/IV_filter_PCR_major4.RData") 

# df_model for paper I --------------------------------------

# filtering milk and SCC values
df1 <- df_all %>%
  filter(SCC > 0) %>%
  filter(SCC < 9999) %>%
  filter(MILK > 0) %>%
  filter(MILK < 100) 

dplyr::n_distinct(df_all$DYR_ID)   # 1551899
dplyr::n_distinct(df_all$BES_ID)   # 3843
dplyr::n_distinct(df1$DYR_ID)   # 1551847 , loosing unique 52 DYR_ID
dplyr::n_distinct(df1$BES_ID)   # 3843, loosing no herds


# lactation pahse 5-305
df2 <- df1 %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) 

dplyr::n_distinct(df_all$DYR_ID)   # 1551899
dplyr::n_distinct(df_all$BES_ID)   # 3843
dplyr::n_distinct(df2$DYR_ID)   # 1551847 , loosing unique 52 DYR_ID
dplyr::n_distinct(df2$BES_ID)   # 3843, loosing no herds


# must be PCR tested, holstein breed (1), Conventionel herdtype (1):
df3 <- df2 %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, PARITY, DIM, SCC, MILK, RES_MAJOR, RES_MINOR) # remove dates and others

dplyr::n_distinct(df3$DYR_ID)   # 214916 
dplyr::n_distinct(df3$BES_ID)   # 1486, loosing no herds


# create logSCC
df4 <- df3 %>%
  mutate(logSCC = log(SCC)) %>% 
  ungroup()


# removing SCC=1 (to avoid logSCC=0)
ggplot(df4, aes(x=SCC)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 colour="black", fill="white") +
  scale_x_continuous(trans="log10") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

# save data with SCC=1 (logSCC=0), before removing:
df_with_SCC1 <- df4

# # remove logSCC= 0. Do not impact model. Keep logSCC=0 (SCC=1)
# df5 <- df4 %>%
#   filter(logSCC > 0)

dplyr::n_distinct(df5$DYR_ID)   # 214916 , loosing none
dplyr::n_distinct(df5$BES_ID)   # 1486, loosing non!
# however loosing: df4_obs - df5_obs = 2753641 - 2752698 = 943 observations


# Parity grouping
# not needed. Do this directly when creting the 6 different dataframes
# df5 <- df4 %>%
#   filter(PARITY > 1) %>%
#   mutate(
#     PARITY = as.numeric(PARITY),
#     PARITY = replace(PARITY, PARITY > 3, 4),
#     PARITY = as.factor(PARITY))


# # create with herds with min 50 animals and 200 obs.
# This is instead done when creating the small data 
# df6 <- df5 %>%
#   group_by(BES_ID) %>%
#   filter(n() > 200)
# df6 <- df5 %>%
#   group_by(BES_ID, DYR_ID) %>%
#   filter(n() > 200)
# create coloumn counting herd occurences
# df <- df %>%
#   group_by(BES_ID, PARITY) %>%
#   mutate(count = n())


# factorize BES_ID. # don't do it with DYR_ID. Not needed for now. Data must be smaller for this
df6 <- df5 %>%
  mutate(BES_ID = factor(BES_ID))

df_model <- df6
summary(df_model)

rm(df1, df2, df3, df4, df5, df6); gc()


# Grouping data for modelling -------------------------------------------


# PARITY 2 - with SCC=1 (logSCC = 0) included ---------------------------

# creting a Parity 2 test data with SCC=1 (logSCC=0)
df_SCC_zero_P2_pos <- df_all %>%
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1) %>%
  filter(SCC > 0) %>%
  filter(SCC < 9999) %>%
  filter(MILK > 0) %>%
  filter(MILK < 100) %>%
  mutate(logSCC = log(SCC)) %>% 
  filter(PARITY == 2) %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>% # not sure if this should be added
  ungroup() %>%
  mutate(BES_ID = factor(BES_ID)) %>%
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)

df_SCC_zero_P2_neg <- df_all %>%
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1) %>%
  filter(SCC > 0) %>%
  filter(SCC < 9999) %>%
  filter(MILK > 0) %>%
  filter(MILK < 100) %>%
  mutate(logSCC = log(SCC)) %>% 
  filter(PARITY == 2) %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>% # not sure if this should be added
  ungroup() %>%
  mutate(BES_ID = factor(BES_ID)) %>%
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)



# PARITY 1 -------------------------------------------

# RES MAJOR POS ; min 200 obs
df1_pos <- df_model %>%
  ungroup() %>%
  filter(PARITY == 1) %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
df2_pos$BES_ID <- factor(df2_pos$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df1_pos$BES_ID) # 557
dplyr::n_distinct(df1_pos$DYR_ID) # 49107
summary(df1_pos)


# RES MAJOR POS ; min 200 obs
df1_neg <- df_model %>%
  ungroup() %>%
  filter(PARITY == 1) %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
df1_neg$BES_ID <- factor(df1_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df1_neg$BES_ID) # 643
dplyr::n_distinct(df1_neg$DYR_ID) # 77811
summary(df1_neg)




# all negative - only for Parity 2 --------------------------

# df2_ALL_NEG
df2_ALL_NEG <- df_model %>%
  ungroup() %>%
  filter(PARITY == 1) %>%
  filter(RES_MAJOR == 0) %>%
  filter(RES_MINOR == 0) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>% # not sure if this should be added
  ungroup() %>%
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
#dplyr::select(BES_ID, DYR_ID, DIM, logSCC)
df2_ALL_NEG$BES_ID <- factor(df2_ALL_NEG$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df2_ALL_NEG$BES_ID) # 41
dplyr::n_distinct(df2_ALL_NEG$DYR_ID) # 1731



# PARITY 2 ----------------------------------------------

# RES MAJOR POS ; min 200 obs
df2_pos <- df_model %>%
  ungroup() %>%
  filter(PARITY == 2) %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
df2_pos$BES_ID <- factor(df2_pos$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df2_pos$BES_ID) # 528
dplyr::n_distinct(df2_pos$DYR_ID) # 41136
summary(df2_pos)

# RES MAJOR NEG ; min 200 obs
df2_neg <- df_model %>%
  ungroup() %>%
  filter(PARITY == 2) %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
df2_neg$BES_ID <- factor(df2_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df2_neg$BES_ID) # 554
dplyr::n_distinct(df2_neg$DYR_ID) # 44210
summary(df2_neg)


# PARITY 3-------------------------------------------
# RES MAJOR POS ; min 200 obs
df3_pos <- df_model %>%
  ungroup() %>%
  filter(PARITY == 3) %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
df3_pos$BES_ID <- factor(df3_pos$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df3_pos$BES_ID) # 414
dplyr::n_distinct(df3_pos$DYR_ID) # 22722
summary(df3_pos)

# RES MAJOR NEG ; min 200 obs
df3_neg <- df_model %>%
  ungroup() %>%
  filter(PARITY == 3) %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC)
df3_neg$BES_ID <- factor(df3_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df3_neg$BES_ID) # 373
dplyr::n_distinct(df3_neg$DYR_ID) # 20063
summary(df3_neg)

# PARITY 4 ------------------------------------------
# RES MAJOR POS ; min 200 obs
df4_pos <- df_model %>%
  ungroup() %>%
  filter(PARITY > 3) %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC, PARITY)
df4_pos$BES_ID <- factor(df4_pos$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df4_pos$BES_ID) # 300
dplyr::n_distinct(df4_pos$DYR_ID) # 12027
summary(df4_pos)

# RES MAJOR NEG ; min 200 obs
df4_neg <- df_model %>%
  ungroup() %>%
  filter(PARITY > 3) %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  #dplyr::select(BES_ID, DIM, logSCC, PARITY)
  dplyr::select(BES_ID, DYR_ID, DIM, logSCC, SCC, PARITY)
df4_neg$BES_ID <- factor(df4_neg$BES_ID) # keep only used levels by resetting the variable

dplyr::n_distinct(df4_neg$BES_ID) # 265
dplyr::n_distinct(df4_neg$DYR_ID) # 10351
summary(df4_neg)


# Saving modelling data ----------------------------------------

save.image("K:/paperI/major4/IV_filter_PCR_major4.RData")
