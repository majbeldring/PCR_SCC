

# Maj Beldring Henningsen, majbh@sund.ku.dk

# wilmink for SCC curves

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

# load("M:/PCR_data/curves_SCC.RData") # from data_tSCC script. Not saved, so run this script for data


#---------------------------------------------------------------------------------
# Wilmink NLS on herd level with equation 8 graesboell:


# STARTLIST defined and equation 8
st <- list(a = 7.9, b = 0.0018, k = -1.6, d = 2.0)
f_nls <- logtSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | BES_ID


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


save.image("M:/PCR_data/wilmink_tSCC_out.RData") 



#--------------------------------------------------------------
# Parameters plot
# for wilmink SCC: only histograms (ggpairs intercorrelation not useful here)




