

# Maj Beldring Henningsen, majbh@sund.ku.dk

# woods for MILK curves

# nlslist for running model

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

# load("M:/PCR_data/curves_MILK.RData") # from data_MILK script - not saved, so must run this script

# note: SCC and milk data are named the same, so don't run woods and wilmink at the same time

#---------------------------------------------------------------------------------
# Woods NLS on herd level:


# woods eq 2 graessboell: 
MILK ~ a * (DIM ^ b) * exp(-(exp(k)) * DIM)

# where (from fathi 2008 - more can be found in article):
# a and b are parameters
## a: scaling factor to present yield in the beginnign of lactation
## b: inclining of slope
## c: declining of slope
# Y_0 = MILK_0 = initial milk yield (kg/day)  = 0
# DIM_m = t_m = time to peak (in days)        = b/c
# MILK_m = y_m = peak yield                   = a(a/b)^c) * exp^(-b)


# STARTLIST 
st <- list(a = 20, b = 0.2, k = -5.6)
f_nls <- MILK ~ a * (DIM ^ b) * exp(-(exp(k)) * DIM) | BES_ID


#---------------------------------------------------------------------------------
# Wilmink NLS on herd level with equation 8 graesboell:
# if parameters not retrieved from .../nls_output.RData


# nls model fit: Parity 2,3,4 and POS + NEg
fit_pos2 <- nlsList(f_nls, df2_pos, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE))
fit_neg2 <- nlsList(f_nls, df2_neg, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE))


fit_pos3 <- nlsList(f_nls, df3_pos, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE))
fit_neg3 <- nlsList(f_nls, df3_neg, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE))


fit_pos4 <- nlsList(f_nls, df4_pos, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE))
fit_neg4 <- nlsList(f_nls, df4_neg, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE))

#------------------------------------------------------------------------------------
# output coefficients / parameters:

# parameters:
out_neg2 <- coef(fit_neg2) %>%
  drop_na()
out_pos2 <- coef(fit_pos2) %>%
  drop_na()

out_neg3 <- coef(fit_neg3) %>%
  drop_na()
out_pos3 <- coef(fit_pos3) %>%
  drop_na()

out_neg4 <- coef(fit_neg4) %>%
  drop_na()
out_pos4 <- coef(fit_pos4) %>%
  drop_na()

# prepare output for plot
out_pos2 %>% class()
out_pos3 %>% class()
out_pos4 %>% class()
out_neg2 %>% class()
out_neg3 %>% class()
out_neg4 %>% class()


save.image("M:/PCR_data/woods_MILK_out.RData") 


#---------------------------------------------------------------------------------
# plotting intercorrelation between woods parameters:

out_pos2 %>% 
  ggpairs()
out_neg2 %>% 
  ggpairs()

out_pos3 %>% 
  ggpairs()
out_neg3 %>% 
  ggpairs()

out_pos4 %>% 
  ggpairs()
out_neg4 %>% 
  ggpairs()


