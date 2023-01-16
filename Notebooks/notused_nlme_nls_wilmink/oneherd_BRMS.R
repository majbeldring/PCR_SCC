
# Maj Beldring Henningsen, majbh@sund.ku.dk

# testing one herd
## brms (bayesian)
## compare with nls and nlme


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
# MCMC bayesian brms

wilmink_priors <- c(
  set_prior("normal(0, 1)", nlpar = "a", lb=0, ub=7),
  set_prior("normal(0, 1)", nlpar = "b", lb=0, ub=1),
  set_prior("normal(0, 1)", nlpar = "k", lb=-5, ub=0), 
  set_prior("normal(0, 1)", nlpar = "d", lb=0, ub=7),
  set_prior("normal(0, 1)", class="sigma"))

wilmink_bay_function <- bf(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                             # Nonlinear variables
                             a + b + k + d ~ 1,
                             # Nonlinear fit
                             nl = TRUE)

wilmink_bay_fit <- brm(wilmink_bay_function,
                       family=gaussian(), 
                       data = df,
                       prior = wilmink_priors )


