
# Maj Beldring Henningsen, majbh@sund.ku.dk

# testing one herd
### test with nls, nlme and brms (bayesian)
### compare the three

# steps:
## 1: preparea and visualize data
## 2: define wilmink function
## 3: run nls with upper and lower (not grouped)
## 4: visualize nls output
## 5: nlme based on nls output
## 6: visualise nlme output; parameter distribution
## 7: create curves based on nlme mean
## 8: min and 100 days values on nlme based curve
## 9: nls grouped with upper and lower - not nlslist
## 10: .md on all


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

Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

load("M:/PCR_data/PCR_merge.RData") 
rm(df_pcr, df_curve); gc() 

# pre-preparing data;remove dates and DYR_ID:
df <- df_model %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1) %>%
  dplyr::select(BES_ID, DYR_ID, PARITY, DIM, logSCC, RES_MAJOR)

# choosing only one random herd: 5051112
df$BES_ID <- factor(df$BES_ID) 
df <- df %>% filter(BES_ID == 5051112)
df <- df %>% 
  ungroup()%>%
  dplyr::select(-BES_ID)
df$DYR_ID <- factor(df$DYR_ID) # removed unused levels by resetting
df$RES_MAJOR <- factor(df$RES_MAJOR) 

# dividing into POS and NEG (all parities; 2,3,4)
df_pos <- df %>%
  filter(RES_MAJOR == 1) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC, PARITY)

df_neg <- df %>%
  filter(RES_MAJOR == 0) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC, PARITY)


#------------------------------------------------------
# visualize data:

# herd 5051112
dplyr::n_distinct(df$DYR_ID) # 337
dplyr::n_distinct(df_pos$DYR_ID) # 160
dplyr::n_distinct(df_neg$DYR_ID) # 257


# checking logSCC distibution:
df %>%
  ggplot( aes(x=logSCC)) +
  geom_histogram( binwidth=0.2, fill="#69b3a2", color="#e9ecef", alpha=0.7) +
  ggtitle("herd 5051112: logSCC distibution parity 2,3,3+") +
  theme(plot.title = element_text(size=10))
# geombar(aes(fill = factor(PARITY)))

# fitting a single curve
ggplot(df, aes(x=DIM, y=logSCC)) +
  geom_point(aes(colour = factor(PARITY)), size=1) +
  ylim(c(0, 10))



#---------------------------------------------------------
# define functions and startlists:

# nls: startlist and function
st <- list(a = 3.86, b = 0.004, k = -2.2, d = 1.25)
f_nls <- logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | BES_ID

# wilmink function for nlme
f_nlme <- function(DIM, a,b,k,d){
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
wilmink_nls <- nls.multstart::nls_multstart(logSCC ~ f_nlme(DIM, a, b, k, d),
                                            data = df,
                                            lower=c(a=0, b=0, k=-5, d=0),
                                            upper=c(a=9, b=1.5, k=0, d=5),
                                            start_lower = c(a=0, b=0, k=-5, d=0),
                                            start_upper = c(a=8, b=1, k=-0.01, d=4),
                                            iter = 500,
                                            supp_errors = "Y")

summary(wilmink_nls)

plot_nls <- function(nls_object, data) {
  predframe <- tibble(DIM=seq(from=min(data$DIM), to=max(data$DIM), 
                               length.out = 1024)) %>%
    mutate(logSCC_nls = predict(nls_object, newdata = list(DIM=.$DIM)))
  ggplot(data, aes(x=DIM, y=logSCC)) +
    geom_point(aes(colour = factor(PARITY)), size=1) +
    ggtitle("herd 5051112: nls fit, not grouped") +
    geom_line(data = predframe, aes(x=DIM, y=logSCC_nls))
}

plot_nls(wilmink_nls, df)




# step 1.a: nls test with minpack.lm - here start parameters must be defined more carefully:
wilmink_nls2 <- minpack.lm::nlsLM(logSCC ~ f_nlme(DIM, a, b, k, d),
                                   data = df,
                                   lower=c(a=0, b=0, k=-5, d=0),
                                   upper=c(a=9, b=1.5, k=0, d=5),
                                   start=c(a=4, b=1, k=-1.8, d=2))




# coefficients from multistart and minpack:

coef(wilmink_nls) # multstart
# a            b            k            d 
# 3.863562327  0.003982492 -2.198826825  1.251174979 

coef(wilmink_nls2) # not multistart
# a            b            k            d 
# 3.863572936  0.003982442 -2.198609374  1.251426471 




# 2.   visualize the nls output - parameter distrubution from multistart fit

wilmink_nls_fit_func <- function(nls_df) {
  nls.multstart::nls_multstart(logSCC ~ f_nlme(DIM, a, b, k, d),
                               data = nls_df,
                               lower=c(a=0, b=0, k=-5, d=0),
                               upper=c(a=9, b=1.5, k=0, d=5),
                               start_lower = c(a=0, b=0, k=-5, d=0),
                               start_upper = c(a=8, b=1, k=-0.01, d=4),
                               iter = 500,
                               supp_errors = "Y")
}


df_nls <- df %>%
  nest(parity_data = c(DYR_ID, DIM, logSCC, RES_MAJOR)) %>% 
  mutate(wilmink_nls = map(parity_data, ~wilmink_nls_fit_func(.x)))

# df$wilmink_nls[[1]]
wilmink_nls_out <- df_nls %>% 
  mutate(outpars = map(wilmink_nls, ~broom::tidy(.x))) %>% 
  select(-wilmink_nls) %>% 
  unnest(cols="outpars")

ggplot(wilmink_nls_out, aes(x=estimate, colour=term, fill=term)) +
  geom_density(alpha=0.5, fill=colourcodes[1], colour=colourcodes[1]) +
  facet_wrap(~term, scales = "free")




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
#   |a    |  3.896|  3.919| 0.222|
#   |b    |  0.004|  0.004| 0.000|
#   |d    |  2.628|  1.759| 2.079|
#   |k    | -1.910| -2.070| 0.867|




#-----------------------------------------------------------
# 5. fit per PARITY and RES_MAJOR:


# PARITY
wilmink_nls_plots <- df %>% 
  nest(parity_data = c(DYR_ID, DIM, logSCC, RES_MAJOR)) %>% 
  select(PARITY, parity_data) %>%
  unnest(parity_data)

wilmink_predtimes <- tidyr::crossing(PARITY=df$PARITY, 
                                  DIM=seq(min(wilmink_nls_plots$DIM),
                                           max(wilmink_nls_plots$DIM),
                                           length.out=128))

wilmink_nlspreds <- wilmink_predtimes %>% 
  group_by(PARITY) %>% 
  nest(preds = DIM) %>% 
  left_join(select(df_nls, PARITY, wilmink_nls)) %>% 
  mutate(preds = map2(preds, wilmink_nls, ~broom::augment(.y, newdata=.x))) %>% 
  select(-wilmink_nls) %>% 
  ungroup() %>% 
  unnest(cols=preds)


ggplot(wilmink_nls_plots, aes(x=DIM, y=logSCC)) +
  geom_point(aes(colour = factor(PARITY)), size=1) +
  geom_line(data=wilmink_nlspreds, aes(y=.fitted),  size=0.8) +
  facet_wrap(~PARITY, ncol=3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# RES_MAJOR
df_nls2 <- df %>%
  nest(res_data = c(DYR_ID, DIM, logSCC, PARITY)) %>% 
  mutate(wilmink_nls = map(res_data, ~wilmink_nls_fit_func(.x)))

wilmink_nls_plots2 <- df %>% 
  nest(res_data = c(DYR_ID, DIM, logSCC, PARITY)) %>% 
  select(RES_MAJOR, res_data) %>%
  unnest(res_data)

wilmink_predtimes2 <- tidyr::crossing(RES_MAJOR=df$RES_MAJOR, 
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










#-----------------------------------------------------
# fitting to nlme with mean values from nls output:

# testing
df_model <- df %>%
  select(PARITY, logSCC, DIM) %>%
  ungroup() 

df_model2 <- df %>%
  select(DYR_ID, logSCC, DIM) %>%
  ungroup()


wilmink_nlme_fit <- nlme(logSCC ~ f_nlme(DIM, a, b, k, d), 
                      data = df,
                      fixed= a + b + k + d ~ 1, 
                      random = a + b + k + d ~ 1, 
                      groups = ~ RES_MAJOR, 
                      start = wilmink_nls_out_summary$mean,
                      verbose = F,
                      na.action=na.exclude,
                      control=(msMaxIter=10000))

df_matt <- df %>%
  filter(PARITY %in% c("2","3")) %>%
  mutate(LactIndex = str_c(DYR_ID, "_", PARITY)) %>%
  group_by(LactIndex) %>%
  mutate(Nobs = n(), FirstObs = min(DIM)) %>%
  ungroup()

df_matt %>% count(LactIndex) %>% count(n)

model_data <- df_matt %>% filter(Nobs %in% 8:9, FirstObs < 20)
model_data <-model_data %>% filter(LactIndex %in% unique(model_data$LactIndex)[1:10])
ggplot(model_data, aes(x=DIM, group=LactIndex, col=LactIndex, y=logSCC)) +
  geom_line() +
  theme(legend.pos="none")

wilmink_nlme_fit <- nlme(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150)*d, 
                         data = model_data,
                         fixed= a + b + k + d ~ 1, 
                         #random = a + b + k + d ~ 1, 
                         random = a + b ~ 1, 
                         groups = ~ LactIndex, 
                         start = wilmink_nls_out_summary$median,
                         na.action=na.exclude,
                         control=(msMaxIter=400))

summary(wilmink_nlme_fit)

# following error can be ocervome by dividing DIM with 150
# error: model may be overfitted
# Error in chol.default((value + t(value))/2) : 
#   the leading minor of order 3 is not positive definite


# error when it's start to run:
# Error in MEestimate(nlmeSt, grpShrunk) : 
#   Singularity in backsolve at level 0, block 1
# In addition: Warning messages:
#   1: In nlme.formula(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150) *  :
#      Iteration 1, LME step: nlminb() did not converge (code = 1). Do increase 'msMaxIter'!
#   2: In nlme.formula(logSCC ~ a + b * DIM/150 + exp(-(exp(k)) * DIM/150) *  :
#      Iteration 2, LME step: nlminb() did not converge (code = 1). Do increase 'msMaxIter'!

summary(wilmink_nlme_fit)






