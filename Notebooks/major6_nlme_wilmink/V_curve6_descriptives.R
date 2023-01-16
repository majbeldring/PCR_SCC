

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Descriptives for Wilmink SCC curves
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

# Data ----------------------------------------------

load("M:/PCR_data/curve6_data_SCC.RData")
# load("M:/PCR_data/curve6_nlme_wilmink.RData") # for results 




df2 %>%
  group_by(RES_MAJOR) %>%
  summarize(min = min(SCC),
            q1 = quantile(SCC, 0.25),
            median = median(SCC),
            mean = mean(SCC),
            q3 = quantile(SCC, 0.75),
            max = max(SCC))

pos2 <- df2 %>%
  filter(RES_MAJOR == 1)
summary(pos2)

library(Hmisc)
describe(pos2) 


# count for article -----------------------------------------

# we also need the animals, so wait untill animals is un-selected

# df2_pos and df2_neg:
summary(df2_pos)
dplyr::n_distinct(df2_pos$BES_ID) # 567
dplyr::n_distinct(df2_pos$DYR_ID) # 997.784

summary(df2_neg)
dplyr::n_distinct(df2_neg$BES_ID) # 522
dplyr::n_distinct(df2_neg$DYR_ID) # 997.784

# df3_pos and df3_neg:
summary(df3_pos)
dplyr::n_distinct(df3_pos$BES_ID) # 567
dplyr::n_distinct(df3_pos$DYR_ID) # 997.784

summary(df3_neg)
dplyr::n_distinct(df3_neg$BES_ID) # 522
dplyr::n_distinct(df3_neg$DYR_ID) # 997.784


# df4_pos and df4_neg:
summary(df4_pos)
dplyr::n_distinct(df4_pos$BES_ID) # 567
dplyr::n_distinct(df4_pos$DYR_ID) # 997.784

summary(df2_neg)
dplyr::n_distinct(df4_neg$BES_ID) # 522
dplyr::n_distinct(df4_neg$DYR_ID) # 997.784




