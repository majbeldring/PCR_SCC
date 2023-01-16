
# Maj Beldring Henningsen, majbh@sund.ku.dk
# 2021-06-23

# Nonlinear modelling:
# Running nlm using the output from nls

# This example is for the wilmink function used for the somatic cell count.


# Packages and settings -------------------------------------------------------


library(tidyverse)
library(ggpubr) # p values in boxplot
library(gridExtra) # gridarrange
ggplot2::theme_set(ggplot2::theme_bw()) 
library(GGally) # for ggpairs
library(nls.multstart)
library(nlme) 

Sys.setlocale("LC_ALL","English") # date formatting
options(stringsAsFactors = FALSE) # prevent factorizing caracters


# Loading data and preparing data -------------------------------------------------

# load("M:/PCR_data/curves_SCC.RData") 

# DATA with three coloumns
## DIM: time variable, continuous
## logSCC: Somatic cell count logtransformed, continuous
## BES_ID: herd ID, factor


# Wilmink function -------------------------------------------------------------------


f_wilmink <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}


# nls multistart (population level) ---------------------------------------------------


# nls.multstart: Robust Non-Linear Regression using AIC Score
nls_pos2 <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                         # population level
                                         data = df2_pos,
                                         # lower & upper boundaries for parameters estimates:
                                         lower=c(a=0, b=0, k=-5, d=0),
                                         upper=c(a=9, b=1.5, k=0, d=5),
                                         # lower & upper boundaries for start parameters:
                                         start_lower = c(a=0, b=0, k=-5, d=0),
                                         start_upper = c(a=8, b=1, k=-0.01, d=4),
                                         iter = 500,
                                         supp_errors = "Y")

nlme_start_pos2 <- coef(nls_pos2) %>% 
  as_tibble() 


# output stats
# coef(nlme_neg2) %>% 
#   as_tibble() %>% 
#   summarise(across(everything(), mean))



# nlme wilmink on herd level ---------------------------------------

nlme_pos2 <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                  data=df2_pos,
                  fixed=a+b+k+d~1,
                  random=a+b+k+d~1,
                  groups=~BES_ID,
                  #start=c(a = 3.5, b = 0.0032, k = -2.3, d = 1.9),
                  start = nlme_start_pos2$value,
                  # start = nlme_start_pos2$median, 
                  # start = nlme_start_pos2$mean,
                  na.action=na.exclude,
                  control = list(maxIter = 1200, msMaxIter = 1200))

out_nlme_pos2 <- coef(nlme_pos2)

# output stats
coef(nlme_pos2) %>% 
  as_tibble() %>% 
  summarise(across(everything(), mean))
  #summarise(across(everything(), median))
  #summarise(across(everything(), 
  #                 list(
  #                   "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))



# output distribtuions, parameter correlations ------------------------------


# Intercorrelation between parameters
gg_pos2 <- out_nlme_pos2 %>% 
  ggpairs()


# visualize nlme parameters boxplot:
p_a <- ggplot(out_nlme_pos2, aes(x=PARITY, y=a, fill=PCR)) + 
  geom_boxplot() +
  # include t.test (p values)
  stat_compare_means(method = "t.test") +
  labs(title="parameter a, PCR POS vs NEG",x="Parity", y = "a values") 



# Histogram
ha_pos2 <- ggplot(out_nlme_pos2, aes(x=a))+
#ha_neg2 <- ggplot(out_nlme_neg2, aes(x=a)) +
  geom_histogram(color="darkred", fill="coral") + # pos colous
  #geom_histogram(color="darkblue", fill="lightblue") + # neg colous
  labs(x="a: Parity 2 POS") +
  #labs(x="a: Parity 2 NEG") +
  xlim(2.0, 6.5)



# percentile curves based on wilmink ---------------------------------------


pos2_pct <- out_nlme_pos2 %>% 
  as_tibble() %>% 
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1))))) %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())


curve_pct_pos2 <- pos2_pct %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      aes(group = pctile) +
      aes(color = pctile) +
      geom_line() +
      
      #labs(caption = "PCR POS, parity 2") +
      ggtitle("SCC curve, pctile of parametres Wilmink: PCR POS parity 2") +
      #ylim(3.0, 11.0) +
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()




# plot mean curves based on wilmink ----------------------------------------


pos2_mean <- out_nlme_pos2 %>% 
  summarise(across(everything(), mean))

curve_mean_pos2 <- pos2_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      #geom_line(size = rel(1.5), colour = "#E69F00") + # orange for POS
      geom_line(size = rel(1.5), colour = "#56B4E9") + # blue for NEG
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("SCC curve NEG parity 2: Mean of Wilmink parametres") +
      #ylim(3.5, 10.5) +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()
  

# retrieve values from fitted curve --------------------------------------------


est_pos2_mean <- 
  pos2_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% 
  identity() %>% {
    bind_rows(
      slice_min(., logSCC, n = 1),
      filter(., DIM == 100),
      filter(., DIM == 150)
    )
  }




# THE END-------------------------------------------------------




