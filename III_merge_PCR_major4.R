#' 
#' 
#' Maj Beldring, majbh@sund.ku.dk
#' UCPH, November 2021
#' 
#' Paper 1, PCR, wilmink: 
#' Merging data
#' Script #3 in PCR project
#'
#'
# Packages and settings: ----------------------------------------

library(tidyverse)
library(lubridate) # for date wrangling
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit: size=56000


# Loading data: 
load("K:/paperI/major4/II_prepare_PCR_major4.RData")
rm(majordry, minor, pcr_full); gc()
rm(dryoff); gc()



# Merging step by step to production data ----------------------------------------

## df1: Production + herd type
df1 <- left_join(production, herd, by = "BES_ID", sort="TRUE",allow.cartesian=TRUE)
df1 <- df1 %>% 
  filter(DATO_FRA <= KONTROLDATO) %>% 
  filter(DATO_TIL >= KONTROLDATO) %>%
  dplyr::select(-DATO_FRA, -DATO_TIL) %>%
  mutate(HERDTYPE = case_when(HERD_TYPE == 'con' ~ 1, HERD_TYPE == 'eco' ~ 0)) %>%
  dplyr::select(-HERD_TYPE)

glimpse(df1) # HERDTYPE: 1=con, 0=eco

dplyr::n_distinct(production$DYR_ID)  # 2524081 unique DYR_ID production
dplyr::n_distinct(df1$DYR_ID)         # 2522761 unique DYR_ID production
dplyr::n_distinct(production$BES_ID)  # 3933
dplyr::n_distinct(df1$BES_ID)         # 3920

rm(herd, production); gc()



# BREED ------------------------------------------
## df2: + breed

df2 <- full_join(df1, breed, by = "DYR_ID", sort="TRUE",allow.cartesian=TRUE)
df2 <- df2 %>% 
  drop_na() %>%
  dplyr::select(-RACE) # 1=holstein, 2=jersey, 3=other dairy breeds

glimpse(df2) # BREED: 1=Holstein, 2=Jersey, 3= other
dplyr::n_distinct(df2$DYR_ID)   # 2520991
dplyr::n_distinct(df2$BES_ID)   # 3918

rm(df1, breed); gc()




# PARITY ------------------------------------------------
## df3 + df4: + calvings (create DIM + parity)

df3 <- full_join(df2, calvings, by = "DYR_ID", sort="TRUE",allow.cartesian=TRUE)

# df4: calving date will be start for each lactation
df4 <- df3 %>% 
  drop_na() %>%
  filter(CALVING_DATE < KONTROLDATO) %>%
  filter(KONTROLDATO - 350 < CALVING_DATE) %>% # setting max lact phase to 350 days
  arrange(DYR_ID, KONTROLDATO, desc(CALVING_DATE)) %>%
  distinct(DYR_ID, KONTROLDATO, .keep_all = TRUE)

glimpse(df4)
# data loss:
dplyr::n_distinct(df4$DYR_ID)   # 2379445
dplyr::n_distinct(df4$BES_ID)   # 3913

rm(df2, df3, calvings); gc()





# DIM ---------------------------------------------------

# add days in milk, 
df4$DIM <- as.Date(as.character(df4$KONTROLDATO), format="%Y-%m-%d")-
  as.Date(as.character(df4$CALVING_DATE), format="%Y-%m-%d")

df4 <- df4 %>% 
  mutate(DIM = as.numeric(DIM))

df5 <- df4 # since we skipped the dryoff step and want the same number of steps as in merge1


# temp saving
save.image("K:/paperI/major4/II_merge_TEMP1.RData") 




# df5 + df6: IMI post calving ---------------------------------------------

# create the IMI post coloumn
df5 <- df5 %>%
  group_by(DYR_ID, PARITY) %>%
  summarise(SCCpost = first(SCC)) %>%
  mutate(SCCpost = lead(SCCpost)) %>%
  left_join(df5, by = c('DYR_ID', 'PARITY')) %>%
  group_by(DYR_ID, PARITY) %>%
  mutate(SCCpost = replace(SCCpost, -1, NA)) %>%
  dplyr::select(-SCCpost, SCCpost)

# df6: now make a coloumn with all IMI or not
IMI <- df5 %>%
  dplyr::select(DYR_ID, PARITY, SCCpost) %>%
  drop_na()
df6 <- df5 %>%
  dplyr::select(-SCCpost)
df6 <- left_join(df6, IMI, by = c("DYR_ID", "PARITY"), sort="TRUE",allow.cartesian=TRUE)

# dropping all observations with no post SCC reported:
# NOTE: this will ONLY keep animals that continue next lactation phase!!!!
# ONLY apply if SCC post calving is needed...
df6 <- df6 %>%
  drop_na()


glimpse(df6) # 
dplyr::n_distinct(df5$DYR_ID)   # 2379445
dplyr::n_distinct(df6$DYR_ID)   # 1551899
dplyr::n_distinct(df5$BES_ID)   # 3913
dplyr::n_distinct(df6$BES_ID)   # 3843 herds;


# temp saving
save.image("K:/paperI/major4/II_merge_TEMP2.RData") 

rm(df4, df5, IMI); gc()
df_lac <- df6
# df6 /df_lac can be used for:
# SCC curves and milk curves without need for diagnostic test or treatment info
# additional info is the SCC status pre previous dryout


# prepare pcr and treatments before merging ------------------------------------
# Here other_treat should be modified to also include maybe..

# 1: create branch of master data frame with only calving date and last control date:
# The last control hereby functions as the an alternative dry off date
branch <- df6 %>% 
  arrange(DYR_ID, PARITY, desc(KONTROLDATO)) %>%
  distinct(DYR_ID, PARITY, .keep_all = TRUE)
branch <- branch %>%
  dplyr::select(DYR_ID, PARITY, KONTROLDATO, CALVING_DATE)


# 2: join each treatment data set with the branch
other_treat <- left_join(branch, other_treat, by = "DYR_ID", sort="TRUE",allow.cartesian=TRUE)
other_treat <- other_treat %>%
  group_by(DYR_ID, PARITY) %>%
  filter(OTHER_AB_DATE > CALVING_DATE |is.na(OTHER_AB_DATE)) %>%
  filter(OTHER_AB_DATE - 400 < CALVING_DATE |is.na(OTHER_AB_DATE)) %>%
  filter(OTHER_AB_DATE >= KONTROLDATO |is.na(OTHER_AB_DATE)) %>%
  filter(OTHER_AB_DATE - 50 < KONTROLDATO |is.na(OTHER_AB_DATE)) %>%
  dplyr::select(DYR_ID, PARITY, OTHER_AB, OTHER_AB_DATE)
other_treat <- other_treat %>%
  filter(OTHER_AB == "1") %>% # keep only AB treatments 
  drop_na()
other_treat <- other_treat %>% 
  arrange(DYR_ID, PARITY, desc(OTHER_AB_DATE)) %>%
  distinct(DYR_ID, PARITY, .keep_all = TRUE)



teat_treat <- left_join(branch, teat_treat, by = "DYR_ID", sort="TRUE",allow.cartesian=TRUE)
teat_treat <- teat_treat %>%
  group_by(DYR_ID, PARITY) %>%
  filter(TEAT_DATE >= KONTROLDATO |is.na(TEAT_DATE)) %>%
  filter(TEAT_DATE > CALVING_DATE |is.na(TEAT_DATE)) %>%
  filter(TEAT_DATE - 400 < CALVING_DATE |is.na(TEAT_DATE)) %>%
  filter(TEAT_DATE - 50 < KONTROLDATO |is.na(TEAT_DATE)) %>%
  dplyr::select(DYR_ID, PARITY, TEAT_TREAT, TEAT_DATE)
teat_treat <- teat_treat %>%
  drop_na()



dryoff_treat <- left_join(branch, dryoff_treat, by = "DYR_ID", sort="TRUE",allow.cartesian=TRUE)
dryoff_treat <- dryoff_treat %>%
  group_by(DYR_ID, PARITY) %>%
  filter(DRYTREAT_DATE >= KONTROLDATO |is.na(DRYTREAT_DATE)) %>%
  filter(DRYTREAT_DATE > CALVING_DATE |is.na(DRYTREAT_DATE)) %>%
  filter(DRYTREAT_DATE - 400 < CALVING_DATE |is.na(DRYTREAT_DATE))%>%
  filter(DRYTREAT_DATE - 50 < KONTROLDATO |is.na(DRYTREAT_DATE)) %>%
  dplyr::select(DYR_ID, PARITY, DRY_TREAT, DRYTREAT_DATE)
dryoff_treat <- dryoff_treat %>%
  drop_na()
 
 
# temp saving
save.image("K:/paperI/major4/II_merge_TEMP3.RData") 

  
pcr <- left_join(branch, pcr, by = "DYR_ID", sort="TRUE",allow.cartesian=TRUE)
pcr <- pcr %>%
  drop_na()
pcr <- pcr %>%
  filter(PCR_DATE > CALVING_DATE |is.na(PCR_DATE)) %>%
  filter(PCR_DATE - 400 < CALVING_DATE |is.na(PCR_DATE)) %>%
  filter(PCR_DATE - 50 < KONTROLDATO |is.na(PCR_DATE)) %>%
  filter(PCR_DATE + 36 > KONTROLDATO |is.na(PCR_DATE)) %>%
  dplyr::select(DYR_ID, PARITY, PCR_DATE, RES_MAJOR, RES_MINOR, PCR_TEST)



# 3: convert treats (and pcr tests to factors):
teat_treat <- teat_treat %>%
  mutate(TEAT_TREAT = factor(TEAT_TREAT))

dryoff_treat <- dryoff_treat %>%
  mutate(DRY_TREAT = factor(DRY_TREAT))

other_treat <- other_treat %>%
  mutate(OTHER_AB = factor(OTHER_AB))

pcr <- pcr %>%
  mutate(PCR_TEST = factor(PCR_TEST)) %>%
  mutate(RES_MINOR = factor(RES_MINOR)) %>%
  mutate(RES_MAJOR = factor(RES_MAJOR))


rm(branch); gc()




# Join treatments to master (df6) ---------------------------------------------

# Teat sealing : df7+8: + teat_treat
df7 <- left_join(df6, teat_treat, by = c("DYR_ID", "PARITY"), sort="TRUE",allow.cartesian=TRUE)
df8 <- df7 %>% 
  arrange(DYR_ID, KONTROLDATO, desc(TEAT_DATE)) %>%
  distinct(DYR_ID, KONTROLDATO, .keep_all = TRUE)

#convert all NA's in TEAT_TREAT to 0 (0= NOT TEat treated, or a treat treat not registered)
df8$TEAT_TREAT = factor(df8$TEAT_TREAT, levels=c(levels(df8$TEAT_TREAT), 0))
df8$TEAT_TREAT[is.na(df8$TEAT_TREAT)] = 0

glimpse(df8) # 
dplyr::n_distinct(df8$DYR_ID)   # 1551899 (# 1.551.121 in df6/df_curve)
dplyr::n_distinct(df8$BES_ID)   # 3843 (same as in df6)


# clean up:
rm(df6, df7, teat_treat); gc()




# Treatment at dry-off ----------------------------------------------------------

## df9+10: + dryoff_treat
df9 <- left_join(df8, dryoff_treat, by = c("DYR_ID", "PARITY"), sort="TRUE",allow.cartesian=TRUE)
# keep only treatments in dry-off period
df10 <- df9 %>% 
  arrange(DYR_ID, KONTROLDATO, desc(DRYTREAT_DATE)) %>%
  distinct(DYR_ID, KONTROLDATO, .keep_all = TRUE)

#convert all NA's in DRY_TEST to 0 (0= NOT PCR tested)
df10$DRY_TREAT = factor(df10$DRY_TREAT, levels=c(levels(df10$DRY_TREAT), 0))
df10$DRY_TREAT[is.na(df10$DRY_TREAT)] = 0
glimpse(df10) # 
dplyr::n_distinct(df10$DYR_ID)   # 1551899
dplyr::n_distinct(df10$BES_ID)   # 3843



# clean up:
rm(df8, df9, dryoff_treat); gc()

# temp saving
save.image("K:/paperI/major4/II_merge_TEMP4.RData") 




# Other AB treatments during dry-off -----------------------------------------------

## df11+12: + other_treat
df11 <- left_join(df10, other_treat, by = c("DYR_ID", "PARITY"), sort="TRUE",allow.cartesian=TRUE)

df12 <- df11 # so mathicng rest of script (had a distinct step here not needed)
 
#convert all NA's in OTHER_AB to 0 (0= NOT AB treated)
df12$OTHER_AB = factor(df12$OTHER_AB, levels=c(levels(df12$OTHER_AB), 0))
df12$OTHER_AB[is.na(df12$OTHER_AB)] = 0


dplyr::n_distinct(df12$DYR_ID)   # 1551899 (no data loss)
dplyr::n_distinct(df12$BES_ID)   # 3843


rm(df10, df11, other_treat); gc()


# PCR tests ---------------------------------------------------------

## df11: +pcr
df13 <- left_join(df12, pcr, by = c("DYR_ID", "PARITY"), sort="TRUE",allow.cartesian=TRUE)

df13 <- df13 %>% 
  arrange(DYR_ID, KONTROLDATO, desc(PCR_DATE)) %>%
  distinct(DYR_ID, KONTROLDATO, .keep_all = TRUE)


#convert all NA's in PCR_TEST to 0 (0= NOT AB treated)
df13$PCR_TEST = factor(df13$PCR_TEST, levels=c(levels(df13$PCR_TEST), 0))
df13$PCR_TEST[is.na(df13$PCR_TEST)] = 0

dplyr::n_distinct(df13$DYR_ID)   # 1551899
dplyr::n_distinct(df13$BES_ID)   # 3843

df_all <- df13
rm(pcr, df12, df13, dryoff); gc()




# saving data: ------------------------------------------------------

# df_lac: full lactation curve including MILK, SCC IMI post calving
# df_all: full lactation incl. MILK, SCC, treatments (teat, dry, other), PCR, IMI pre
# major, minor: Pathogen, date, DYR_ID, PCR Ct value
# pcr_full: major and mino combined without major/minor groups

# for modeliing in PCR projetc: use df_all

save.image("K:/paperI/major4/III_merge_PCR_major4.RData") 


