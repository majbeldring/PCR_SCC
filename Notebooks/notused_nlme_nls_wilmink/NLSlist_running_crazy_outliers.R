

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

load("M:/PCR_data/curves_SCC.RData") # from PCR_data script


#---------------------------------------------------------------------------------
# Wilmink NLS on herd level with equation 8 graesboell:


# STARTLIST defined and equation 8
st <- list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6)
f_nls <- logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | BES_ID


#--------------------------------------------------------------------------------
# test with upper and lower boundries:
upper <- c(8,1,0,4)
lower <- c(0,0,-4,0)

f_neg4 <- nlsList(f_nls, df4_neg, start = sapply(st, mean),
                    control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                                   printEval = FALSE, warnOnly = TRUE), upper = upper, lower = lower)

out <- coef(f_neg4)
out <- cbind(BES_ID = rownames(out), out)

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
out_pos2 <- coef(fit_pos2)
out_neg2 <- coef(fit_neg2)
out_pos3 <- coef(fit_pos3)
out_neg3 <- coef(fit_neg3)
out_pos4 <- coef(fit_pos4) 
out_neg4 <- coef(fit_neg4)

out_pos2 <- cbind(BES_ID = rownames(out_pos2), out_pos2)
out_neg2 <- cbind(BES_ID = rownames(out_neg2), out_neg2)
out_pos3 <- cbind(BES_ID = rownames(out_pos3), out_pos3)
out_neg3 <- cbind(BES_ID = rownames(out_neg3), out_neg3)
out_pos4 <- cbind(BES_ID = rownames(out_pos4), out_pos4)
out_neg4 <- cbind(BES_ID = rownames(out_neg4), out_neg4)

# evaluating parameters:

summary(out_pos2) 
# mean vs median: 
# a+b: okay, 
# k: 10^4 diff, 
# d: 10^7 diff, 
# k+d: crazy low and high outliers
summary(out_neg2) 
# mean vs median: 
# a: *3, high outliers. b: ok, 
# k: 10^3 diff, crazy high outlier
# d: 10^5 diff: CRAZY high lower and higher outliers
summary(out_pos3) 
#mean vs median:
# a: diff ok, but high outlier (*30). b: ok
# k: diff 10^5 , crazy high and lower outliers. d: 10^8 diff, crazy high and low
summary(out_neg3) 
# mean vs median:
# a: diff ok, 10^2 outliers high and low. b: ok
# k: diff 10^4 , crazy high and lower outliers. d: 10^6 diff, crazy high and low
summary(out_pos4) 
# mean vs median:
# a: diff *2, igh outlier (10^3 ). b: *10 diff, low outlier
# k: diff 10^5 , crazy high and lower outliers. d: 10^7 diff, crazy high and low
summary(out_neg4) 
# mean vs median:
# a: diff ok, but low and high outlier (+/- *30). b: ok
# d: diff 10^5 , crazy high and lower outliers. d: 10^8 diff, crazy high and low

# overall: 
# a: not crazy, but some differences should be investigated
# b: overall okay
# k: crazy in all output
# d: crazy in all output


# prepare output for plot
out_pos2 %>% class()
out_pos3 %>% class()
out_pos4 %>% class()
out_neg2 %>% class()
out_neg3 %>% class()
out_neg4 %>% class()


save.image("M:/PCR_data/wilmink_SCC_out.RData") 



#--------------------------------------------------------------
# Parameters plot
# for wilmink SCC: only histograms (ggpairs intercorrelation not useful here)




