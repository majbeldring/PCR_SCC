
# Maj Beldring Henningsen, majbh@sund.ku.dk



# follow:
# fitting nlme with nls output:
# https://www.granvillematheson.com/post/nonlinear-modelling-using-nls-nlme-and-brms/
# 
# https://quantdev.ssri.psu.edu/sites/qdev/files/GMChp10_Tutorial.html

# nls grouped, with upper and lower (witput nlslist, just nls):
# https://stackoverflow.com/questions/27570758/r-how-to-use-bounds-and-the-port-algorithm-in-nlslist


#-------------------------------------------------------
# Packages
library(tidyverse)
library(nls.multstart)
library(nlme) # for nlslist
library(brms)
library(broom)
library(hrbrthemes)
library(broom)
library(viridis)
library(purrr)

# settings:
colourcodes <- c("#d4a665", "#d27fff", "#7fd9ff")
colourpal <- c(NLS="#d4a665", NLME="#d27fff", MCMC="#7fd9ff")
theme_set(hrbrthemes::theme_ipsum_rc())

# Sys.setlocale("LC_ALL","English") # date formatting
# memory.size()            # Checking your memory size
# memory.limit()           # Checking the set limit
# memory.limit(size=56000) # suggest for 64 bit
# options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

load("M:/PCR_data/PCR_merge.RData") 
rm(df_pcr, df_curve); gc() 

# one herd (holstein, conventional, Parity 2, herd no. 5051112)
df <- df_model %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1) %>%
  filter(PARITY == 2) %>%
  filter(BES_ID == 5051112) %>%
  #ungroup() %>%
  group_by(DYR_ID) %>%
  mutate(count = n()) %>%
  filter(count > 8) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC, RES_MAJOR)

df_pos <- df %>%
  filter(RES_MAJOR == 1) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC)

df_neg <- df %>%
  filter(RES_MAJOR == 0) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC)

df$DYR_ID <- factor(df$DYR_ID) # removed unused levels by resetting
df$RES_MAJOR <- factor(df$RES_MAJOR) 
df_pos$DYR_ID <- factor(df_pos$DYR_ID)
df_neg$DYR_ID <- factor(df_neg$DYR_ID)

#--------------------------------------------------------
# 8 animals POS and 8 animals NEG :

# POS
# 1008163608 , 1008634327 , 1010747294 , 	1011189781 , 
# 1011520009 , 1011720847 , 1010687774 , 	1008671172

# NEG
# 1008376733 , 1008422731 , 1008534128 , 1008688710 , 
# 1009412981 , 	1009738457 , 1009911167 , 	1010141028

#--------------------------------------------------------
# check data

ggplot(df_pos, aes(x=DIM, y=logSCC)) +
  geom_point(size=1) +
  ylim(c(0, 10))

ggplot(df_neg, aes(x=DIM, y=logSCC)) +
  geom_point(size=1) +
  ylim(c(0, 10))



#---------------------------------------------------------
# define functions and startlists:

# nls: startlist and function
st <- list(a = 3.86, b = 0.004, k = -2.2, d = 1.25)
f_nls <- logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | DYR_ID

# wilmink function for nlme
f_wilmink <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}


# upper and lower boundaries for startlist
upper_bound <- c(9,1.5,0,5)
lower_bound <- c(0,0,-5,0)


# define wilmink eq 8 for recreating curves:
logSCC_func <- function(DIM, a,b,k,d) 
  a + b * DIM + exp(-(exp(k)) * DIM)*d 




#-------------------------------------------------------
# step 1: running nls with multistart using:
# https://www.granvillematheson.com/post/nonlinear-modelling-using-nls-nlme-and-brms/


# nls multistart: THIS IS NOT gouped by DYR_ID, we just do nls to get start parametres for nlme
nls_pos <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                            data = df_neg,
                                            lower=c(a=0, b=0, k=-5, d=0),
                                            upper=c(a=9, b=1.5, k=0, d=5),
                                            start_lower = c(a=0, b=0, k=-5, d=0),
                                            start_upper = c(a=8, b=1, k=-0.01, d=4),
                                            iter = 500,
                                            supp_errors = "Y")

summary(nls_pos)

plot_nls <- function(nls_object, data) {
  predframe <- tibble(DIM=seq(from=min(data$DIM), to=max(data$DIM), 
                               length.out = 1024)) %>%
    mutate(logSCC_nls = predict(nls_object, newdata = list(DIM=.$DIM)))
  ggplot(data, aes(x=DIM, y=logSCC)) +
    geom_point(aes(colour = factor(RES_MAJOR)), size=1) +
    ggtitle("herd 5051112, nls fit, not grouped") +
    geom_line(data = predframe, aes(x=DIM, y=logSCC_nls))
}

plot_nls(wilmink_nls, df)




# step 1.a: nls test with minpack.lm - here start parameters must be defined more carefully:
wilmink_nls2 <- minpack.lm::nlsLM(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                   data = df_pos,
                                   lower=c(a=0, b=0, k=-5, d=0),
                                   upper=c(a=9, b=1.5, k=0, d=5),
                                   start=c(a=4, b=1, k=-1.8, d=2))




# coefficients from multistart and minpack:

coef(wilmink_nls) # multstart
# a            b            k            d 
# 3.663939435  0.004052322 -2.686848246  1.124153656 

coef(wilmink_nls2) # not multistart
# a            b            k            d 
# 3.664002453  0.004052034 -2.686198146  1.124593718 




# 2.   visualize the nls output - parameter distrubution from multistart fit

wilmink_nls_fit_func <- function(nls_df) {
  nls.multstart::nls_multstart(logSCC ~ f_nlme(DIM, a, b, k, d),
                               data = df_pos,
                               lower=c(a=0, b=0, k=-5, d=0),
                               upper=c(a=9, b=1.5, k=0, d=5),
                               start_lower = c(a=0, b=0, k=-5, d=0),
                               start_upper = c(a=8, b=1, k=-0.01, d=4),
                               iter = 500,
                               supp_errors = "Y")
}


df_nls <- df %>%
  nest(res_data = c(DYR_ID, DIM, logSCC)) %>% 
  mutate(wilmink_nls = map(res_data, ~wilmink_nls_fit_func(.x)))

# df2$wilmink_nls[[1]]
wilmink_nls_out <- df_nls %>% 
  mutate(outpars = map(wilmink_nls, ~broom::tidy(.x))) %>% 
  select(-wilmink_nls) %>% 
  unnest(cols="outpars")

ggplot(wilmink_nls_out, aes(x=estimate, colour=term, fill=term)) +
  geom_density(alpha=0.5, fill=colourcodes[1], colour=colourcodes[1]) +
  facet_wrap(~DYR_ID, scales = "free")




# . table with parameter summaries : used for nlme startlist

wilmink_nls_out_summary <- wilmink_nls_out%>%
  group_by(term) %>% 
  summarise(mean = mean(estimate), 
            median = median(estimate),
            sd = sd(estimate)) %>% 
  ungroup()

knitr::kable(wilmink_nls_out_summary, digits = 3)

# |term |   mean| median|    sd|
#   |:----|------:|------:|-----:|
#   |a    |  3.749|  3.749| 0.341|
#   |b    |  0.004|  0.004| 0.000|
#   |d    |  0.926|  0.926| 0.964|
#   |k    | -2.881| -2.881| 0.385|




#-----------------------------------------------------------
# 5. fit per RES_MAJOR:

# RES_MAJOR
df_nls2 <- df2 %>%
  nest(res_data = c(DYR_ID, DIM, logSCC)) %>% 
  mutate(wilmink_nls = map(res_data, ~wilmink_nls_fit_func(.x)))

wilmink_nls_plots2 <- df2 %>% 
  nest(res_data = c(DYR_ID, DIM, logSCC)) %>% 
  select(RES_MAJOR, res_data) %>%
  unnest(res_data)

wilmink_predtimes2 <- tidyr::crossing(RES_MAJOR=df2$RES_MAJOR, 
                                     DIM=seq(min(wilmink_nls_plots2$DIM),
                                             max(wilmink_nls_plots2$DIM),
                                             length.out=128))

wilmink_nlspreds2 <- wilmink_predtimes2 %>% 
  group_by(RES_MAJOR) %>% 
  nest(preds = DIM) %>% 
  left_join(select(df_nls2, RES_MAJOR, wilmink_nls)) %>% 
  mutate(preds = map2(preds, wilmink_nls, ~broom::augment(.y, newdata=.x))) %>% 
  select(-wilmink_nls) %>% 
  ungroup() %>% 
  unnest(cols=preds)


ggplot(wilmink_nls_plots2, aes(x=DIM, y=logSCC)) +
  geom_point(aes(colour = factor(RES_MAJOR)), size=1) +
  geom_line(data=wilmink_nlspreds2, aes(y=.fitted),  size=0.8) +
  facet_wrap(~RES_MAJOR, ncol=2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#-------------------------------------------------------------------------------------
# EXTRA withdraw min values and values at day 100 on the graphs










#----------------------------------------------------------------------------------------
# nlme: doesn't work

# Error in nlme.formula(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150) *  : 
#                         Singularity in backsolve at level 0, block 1
# In addition: Warning message:
# In nlme.formula(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150) *  :
#   Iteration 1, LME step: nlminb() did not converge (code = 1). Do increase 'msMaxIter'!


wilmink_nlme_fit <- nlme(logSCC ~ f_nlme(DIM, a, b, k, d), 
                         data = df2,
                         fixed= a + b + k + d ~ 1, 
                         random = a + b + k + d ~ 1, 
                         groups = ~ RES_MAJOR, 
                         start = wilmink_nls_out_summary$median,
                         verbose = F,
                         na.action=na.exclude,
                         control=(msMaxIter=10000))

wilmink_nlme_fit <- nlme(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150)*d, 
                         data = df2,
                         fixed= a + b + k + d ~ 1, 
                         random = a + b + k + d ~ 1, 
                         groups = ~ DYR_ID, 
                         start = wilmink_nls_out_summary$median,
                         na.action=na.exclude,
                         control=(msMaxIter=10000))


summary(wilmink_nlme_fit)









