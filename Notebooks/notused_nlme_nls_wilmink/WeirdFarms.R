

# Maj Beldring Henningsen, majbh@sund.ku.dk

# WeirdFarms:

# from list with outliers:

# create:
## 2 list for each parameter - one with upper outliers and one with lower
## take lower 10 percent
## 8 lists in total

## 1. create boxplot for every 10 day;
## 2. ecdf curve


#-------------------------------------------------------
# Packages and settings:

library(tidyverse)
library(gridExtra)
library(data.table)
library(plotly)
library(GGally)
#library(tidymodels)
#library(nlstools) # for bootstrapping
#library(nlme) # for nlslist
#library(ggExtra)
#library(ggalluvial)
Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

# load: 
load("M:/PCR_data/wilmink_SCC_out.RData") # loading parameter output from PCR_nlslist


#---------------------------------------------------------------------------------

# if needed:
out_4pos$BES_ID <- rownames(out_4pos)

#------------------------------------------------------------------------------------
# output

lines(predict(fit_pos4)~df_pos$DIM,col="purple")

# plot : NOT WORKING for any due to inverse zero error
plot_nls(fit_pos2) # Not working
tidy(fit_pos2) # not working
tidy(fit_pos4, conf.int=TRUE) # not working



# identify outliers in parameter output:
boxplot.stats(out_4pos$a)$out
boxplot.stats(out_4pos$b)$out
boxplot.stats(out_4pos$k)$out
boxplot.stats(out_4pos$d)$out

out_a1 <- boxplot.stats(out_4pos$a)$out
out_a <- which(out_4pos$a %in% c(out_a1)); out_a
out_b1 <- boxplot.stats(out_4pos$b)$out
out_b <- which(out_4pos$b %in% c(out_b1)); out_b
out_k1 <- boxplot.stats(out_4pos$k)$out
out_k <- which(out_4pos$k %in% c(out_k1)); out_k
out_d1 <- boxplot.stats(out_4pos$d)$out
out_d <- which(out_4pos$d %in% c(out_d1)); out_d



# manuel filter for outliers
filter_4pos <- out_4pos %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)



#-------------------------------------------------------
# plot the crazy parameters:



df_51 <- df4_pos %>%
  filter(BES_ID == "5181412")
p_51 <- ggplot(df_51, aes(x = DIM, y = logSCC)) + 
  ggtitle("5181412 pos Parity >3") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_point() +
  geom_smooth() 
p_51


# farm: 4942112 (243 obs for df_pos)
# pos
df_4942 <- df4_pos %>%
  filter(BES_ID == "4942112")
p_49 <- ggplot(df_4942, aes(x = DIM, y = logSCC)) + 
  ggtitle("4942112 pos Parity >3") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_49
# neg:
df2_4942 <- df4_neg %>%
  filter(BES_ID == "4942112")
p2_49 <- ggplot(df2_4942, aes(x = DIM, y = logSCC)) + 
  ggtitle("4942112 neg Parity >3") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p2_49



# farm: 4942112 (243 obs for df_pos) : THIS FARM HAS NO NEG tests!
# pos
df_322 <- df4_pos %>%
  filter(BES_ID == "3225812")
p_322 <- ggplot(df_322, aes(x = DIM, y = logSCC)) + 
  ggtitle("3225812 pos - NO NEG tests") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_322
# neg:
df2_322 <- df4_neg %>%
  filter(BES_ID == "3225812")
p2_322 <- ggplot(df2_322, aes(x = DIM, y = logSCC)) + 
  ggtitle("3225812 neg") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p2_322




# farm: 4942112 (243 obs for df_pos) : THIS FARM HAS NO NEG tests!
# pos
df_35 <- df4_pos %>%
  filter(BES_ID == "3503112")
p_35 <- ggplot(df_35, aes(x = DIM, y = logSCC)) + 
  ggtitle("3503112 pos") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_35
# neg:
df2_35 <- df4_neg %>%
  filter(BES_ID == "3503112")
p2_35 <- ggplot(df2_35, aes(x = DIM, y = logSCC)) + 
  ggtitle("3503112 neg") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p2_35




# plotting all the weird k parameter farms:
# pos
df_weird <- df4_pos %>%
  filter(BES_ID == "4942112" | BES_ID == "3225812" | BES_ID == "3503112"  | BES_ID == "3888412" |
           BES_ID == "2278512" | BES_ID == "2473112" | BES_ID ==  "4114512" | BES_ID == "4525012" | 
           BES_ID == "4196512"  | BES_ID == "5119312" | BES_ID == "4202812	" | BES_ID == "3369712" |
           BES_ID == "4965912	" | BES_ID == "3507612" | BES_ID == "3804112" | BES_ID == "4949012")
p_weird <- ggplot(df_weird, aes(x = DIM, y = logSCC)) + 
  ggtitle("weird pos") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_weird
# neg:
df2_weird <- df4_neg %>%
  filter(BES_ID == "4942112" | BES_ID == "3225812" | BES_ID == "3503112"  | BES_ID == "3888412" |
           BES_ID == "2278512" | BES_ID == "2473112" | BES_ID ==  "4114512" | BES_ID == "4525012" | 
           BES_ID == "4196512"  | BES_ID == "5119312" | BES_ID == "4202812	" | BES_ID == "3369712" |
           BES_ID == "4965912	" | BES_ID == "3507612" | BES_ID == "3804112" | BES_ID == "4949012")
p2_weird <- ggplot(df2_weird, aes(x = DIM, y = logSCC)) + 
  ggtitle("weird neg") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p2_weird


# testing outlier detection:
out_k <- which(out_4pos$k %in% c(out_k1))
out_k
