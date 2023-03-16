
# PCR project
# Correlation between parameters


# settings and data ------------------------------
library("tidyverse")
library(GGally)

(load("K:/paperI/major4/Ouput_nlme.RData"))



# ggpairs ----------------------------------

# NEG parity 2
nlme_out_neg2 %>% 
  ggpairs()

