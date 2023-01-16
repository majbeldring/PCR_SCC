

# Maj Beldring Henningsen, majbh@sund.ku.dk

# outliers removal and identification for creating weird farms curves

### removal:
# check guide in tools -> outliers
# remove lower and higher qauntile and test again with 


### Outlier identification
# prep for visualizing weird farms

# from list with outliers:

# create:
## 2 list for each parameter - one with upper outliers and one with lower
## take lower 10 percent
## 8 lists in total

## 1. create boxplot for every 10 day;
## 2. ecdf curve


# GOAL:
## Create a weird farm curve for each data set for parameter k (and a and d?) = 12 curves

# Method:
## grepl; to retrive data for only these farms in each data: df_pos2 etc...
## boxplot per 10 DIM for each 

#-------------------------------------------------------
# Packages and settings:

library(tidyverse)

Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

# load: 
load("M:/PCR_data/wilmink_SCC_out.RData") # loading parameter output from PCR_nlslist
rm(fit_neg2, fit_neg3, fit_neg4, fit_pos2, fit_pos3, fit_pos4); gc()

# create DIM interval
df2_pos <- df2_pos %>%
  mutate(DIM_int = cut(DIM, 
                 #Add 1 to the maximum value in dim to make sure it is included in the categorization.
                 breaks = seq(min(DIM),max(DIM)+1,10),
                 #Set this to T to include the lowest value
                 include.lowest = T,
                 #To set labels as a sequence of integers
                 labels = F)) %>%
  mutate(DIM_int = factor(DIM_int))

df2_neg <- df2_neg %>%
  mutate(DIM_int = cut(DIM, 
                       #Add 1 to the maximum value in dim to make sure it is included in the categorization.
                       breaks = seq(min(DIM),max(DIM)+1,10),
                       #Set this to T to include the lowest value
                       include.lowest = T,
                       #To set labels as a sequence of integers
                       labels = F)) %>%
  mutate(DIM_int = factor(DIM_int))

df3_pos <- df3_pos %>%
  mutate(DIM_int = cut(DIM, 
                       breaks = seq(min(DIM),max(DIM)+1,10),
                       include.lowest = T,
                       labels = F)) %>%
  mutate(DIM_int = factor(DIM_int))

df3_neg <- df3_neg %>%
  mutate(DIM_int = cut(DIM, 
                       breaks = seq(min(DIM),max(DIM)+1,10),
                       include.lowest = T,
                       labels = F)) %>%
  mutate(DIM_int = factor(DIM_int))

df4_pos <- df4_pos %>%
  mutate(DIM_int = cut(DIM, 
                       breaks = seq(min(DIM),max(DIM)+1,10),
                       include.lowest = T,
                       labels = F)) %>%
  mutate(DIM_int = factor(DIM_int))

df4_neg <- df4_neg %>%
  mutate(DIM_int = cut(DIM, 
                       breaks = seq(min(DIM),max(DIM)+1,10),
                       include.lowest = T,
                       labels = F)) %>%
  mutate(DIM_int = factor(DIM_int))



#---------------------------------------------------------------------------------
# 1 percent extreme MAX outliers:

## Parity 2, POS, 1 pct max values
maxa_pos2 <- out_pos2 %>%  slice_max(out_pos2$a, prop = 0.01) # a 
maxk_pos2 <- out_pos2 %>%  slice_max(out_pos2$k, prop = 0.01) # k 
maxd_pos2 <- out_pos2 %>%  slice_max(out_pos2$d, prop = 0.01) # d

## Parity 2, NEG, 1 pct max values
maxa_neg2 <- out_neg2 %>% slice_max(out_neg2$a, prop = 0.01) # a
maxk_neg2 <- out_neg2 %>% slice_max(out_neg2$k, prop = 0.01) # k
maxd_neg2 <- out_neg2 %>% slice_max(out_neg2$d, prop = 0.01) # d

## Parity 3, POS, 1 pct max values
maxa_pos3 <- out_pos3 %>% slice_max(out_pos3$a, prop = 0.01) # a
maxk_pos3 <- out_pos3 %>% slice_max(out_pos3$k, prop = 0.01) # k
maxd_pos3 <- out_pos3 %>% slice_max(out_pos3$d, prop = 0.01) # d

## Parity 3, NEG, 1 pct max values
maxa_neg3 <- out_neg3 %>% slice_max(out_neg3$a, prop = 0.01) # a
maxk_neg3 <- out_neg3 %>% slice_max(out_neg3$k, prop = 0.01) # k
maxd_neg3 <- out_neg3 %>% slice_max(out_neg3$d, prop = 0.01) # d

## Parity 4, POS, 1 pct max values
maxa_pos4 <- out_pos4 %>% slice_max(out_pos4$a, prop = 0.01) # a
maxk_pos4 <- out_pos4 %>% slice_max(out_pos4$k, prop = 0.01) # k
maxd_pos4 <- out_pos4 %>% slice_max(out_pos4$d, prop = 0.01) # d

## Parity 4, NEG, 1 pct max values
maxa_neg4 <- out_neg4 %>% slice_max(out_neg4$a, prop = 0.01) # a
maxk_neg4 <- out_neg4 %>% slice_max(out_neg4$k, prop = 0.01) # k
maxd_neg4 <- out_neg4 %>% slice_max(out_neg4$d, prop = 0.01) # d


#-----------------------------------------------------------------------------------
# create dataframes only with outsiders:

maxa_pos2[,"BES_ID"] # "4528812" "1783912" "5861412" "3892412" "3308412"
maxa_pos2 <- dplyr::filter(df2_pos, grepl('4528812|1783912|5861412|3892412|3308412', BES_ID))
maxk_pos2[,"BES_ID"] # "4202812" "4980012" "5860612" "4984212" "2460812"
maxk_pos2 <- dplyr::filter(df2_pos, grepl('4202812|4980012|5860612|4984212|2460812', BES_ID))
maxd_pos2[,"BES_ID"] # "2263312" "4438812" "3077012" "6542912" "4202812"
maxd_pos2 <- dplyr::filter(df2_pos, grepl('2263312|4438812|3077012|6542912|4202812', BES_ID))

maxa_neg2[,"BES_ID"] # "3732812" "4628312" "3230912" "5089412" "4951412"
maxa_neg2 <- dplyr::filter(df2_neg, grepl('3732812|4628312|3230912|5089412|4951412', BES_ID))
maxk_neg2[,"BES_ID"] # "4457912" "4263112" "1931612" "3732812" "3483912"
maxk_neg2 <- dplyr::filter(df2_neg, grepl('4457912|4263112|1931612|3732812|3483912', BES_ID))
maxd_neg2[,"BES_ID"] # "4263112" "4457912" "4739812" "2111012" "1931612"
maxd_neg2 <- dplyr::filter(df2_neg, grepl('4263112|4457912|4739812|2111012|19316122', BES_ID))

maxa_pos3[,"BES_ID"] # "4579412" "5222812" "4209712" "3707012"
maxa_pos3 <- dplyr::filter(df3_pos, grepl('4579412|5222812|4209712|3707012', BES_ID))
maxk_pos3[,"BES_ID"] # "5965412" "2075212" "4158812" "4821712"
maxk_pos3 <- dplyr::filter(df3_pos, grepl('5965412|2075212|4158812|4821712', BES_ID))
maxd_pos3[,"BES_ID"] # "4386812" "5965412" "5731812" "5907712"
maxd_pos3 <- dplyr::filter(df3_pos, grepl('4386812|5965412|5731812|5907712', BES_ID))

maxa_neg3[,"BES_ID"] # "4735512" "1721212" "5660112"
maxa_neg3 <- dplyr::filter(df3_neg, grepl('4735512|1721212|5660112', BES_ID))
maxk_neg3[,"BES_ID"] # "2462112" "6496312" "2279712"
maxk_neg3 <- dplyr::filter(df3_neg, grepl('2462112|6496312|2279712', BES_ID))
maxd_neg3[,"BES_ID"] # "5860612" "3388512" "2279712"
maxd_neg3 <- dplyr::filter(df3_neg, grepl('5860612|3388512|2279712', BES_ID))

maxa_pos4[,"BES_ID"] # "5076512" "3750012" "5162212"
maxa_pos4 <- dplyr::filter(df4_pos, grepl('5076512|3750012|5162212', BES_ID))
maxk_pos4[,"BES_ID"] # "5181412" "5127112" "4158812"
maxk_pos4 <- dplyr::filter(df4_pos, grepl('5181412|5127112|4158812', BES_ID))
maxd_pos4[,"BES_ID"] # "5760812" "5181412" "5251012"
maxd_pos4 <- dplyr::filter(df4_pos, grepl('5760812|5181412|5251012', BES_ID))

maxa_neg4[,"BES_ID"] # "3315412" "3304412"
maxa_neg4 <- dplyr::filter(df4_neg, grepl('3315412|3304412', BES_ID))
maxk_neg4[,"BES_ID"] # "3267112" "4278212"
maxk_neg4 <- dplyr::filter(df4_neg, grepl('3267112|4278212', BES_ID))
maxd_neg4[,"BES_ID"] # "2587212" "4650012"
maxd_neg4 <- dplyr::filter(df4_neg, grepl('2587212|4650012', BES_ID))



#------------------------------------------------------------------------------------
# plot the new max dataframes as boxplot

# df2 pos:
# reference: all df
ggplot(df2_pos, aes(x=DIM_int, y=logSCC)) +
  geom_boxplot() + 
  labs(title="parity 2, POS, reference", x="10 days intervals, starting DIM=6", y = "logSCC")
# a, k, d
ggplot(maxa_pos2, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 2, POS, max parameter a",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxk_pos2, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 2, POS, max parameter k",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxd_pos2, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 2, POS, max parameter d",x="10 days intervals, starting DIM=6", y = "logSCC")


# df2 neg:
# reference: all df
ggplot(df2_neg, aes(x=DIM_int, y=logSCC)) +
  geom_boxplot() + 
  labs(title="parity 2, NEG, reference", x="10 days intervals, starting DIM=6", y = "logSCC")
# a, k, d
ggplot(maxa_neg2, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 2, NEG, max parameter a",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxk_neg2, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 2, NEG, max parameter k",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxd_neg2, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 2, NEG, max parameter d",x="10 days intervals, starting DIM=6", y = "logSCC")



# df3 pos:
# reference: all df
ggplot(df3_pos, aes(x=DIM_int, y=logSCC)) +
  geom_boxplot() + 
  labs(title="parity 3, POS, reference", x="10 days intervals, starting DIM=6", y = "logSCC")
# a, k, d
ggplot(maxa_pos3, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 3, POS, max parameter a",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxk_pos3, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 3, POS, max parameter k",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxd_pos3, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 3, POS, max parameter d",x="10 days intervals, starting DIM=6", y = "logSCC")


# df3 neg:
# reference: all df
ggplot(df3_neg, aes(x=DIM_int, y=logSCC)) +
  geom_boxplot() + 
  labs(title="parity 3, NEG, reference", x="10 days intervals, starting DIM=6", y = "logSCC")
# a, k, d
ggplot(maxa_neg3, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 3, NEG, max parameter a",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxk_neg3, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 3, NEG, max parameter k",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxd_neg3, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 3, NEG, max parameter d",x="10 days intervals, starting DIM=6", y = "logSCC")



# df4 pos:
# reference: all df
ggplot(df4_pos, aes(x=DIM_int, y=logSCC)) +
  geom_boxplot() + 
  labs(title="parity 4, POS, reference", x="10 days intervals, starting DIM=6", y = "logSCC")
# a, k, d
ggplot(maxa_pos4, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, POS, max parameter a",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxk_pos4, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, POS, max parameter k",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxd_pos4, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, POS, max parameter d",x="10 days intervals, starting DIM=6", y = "logSCC")


# df4 neg:
# reference: all df
ggplot(df4_neg, aes(x=DIM_int, y=logSCC)) +
  geom_boxplot() + 
  labs(title="parity 4, NEG, reference", x="10 days intervals, starting DIM=6", y = "logSCC")
# a, k, d
ggplot(maxa_neg4, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, NEG, max parameter a",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxk_neg4, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, NEG, max parameter k",x="10 days intervals, starting DIM=6", y = "logSCC")
ggplot(maxd_neg4, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, NEG, max parameter d",x="10 days intervals, starting DIM=6", y = "logSCC")





# some specific herds:
no_neg_pcr <- dplyr::filter(df4_pos, grepl('3503112', BES_ID))

ggplot(no_neg_pcr, aes(x=DIM_int, y=logSCC, fill=DIM_int)) +
  geom_boxplot() + 
  theme(legend.position="none") + 
  labs(title="parity 4, POS, no NEG PCR",x="10 days intervals, starting DIM=6", y = "logSCC")




## ecdf plot:
# crazy k and d herd 5181412 from df4_pos:
crazy_5181412 <- dplyr::filter(df4_pos, grepl('5181412', BES_ID))

ggplot(crazy_5181412, aes(logSCC)) + 
  stat_ecdf(geom = "step") + 
  xlim(0.0, 9.0) +
  theme(legend.position="none") + 
  theme_bw() + 
  labs(title="parity 4 pos, crazy k and d herd 5181412")

ggplot(df4, aes(logSCC)) + 
  stat_ecdf(geom = "step") + 
  xlim(0.0, 9.0) +
  theme(legend.position="none") + 
  theme_bw() + 
  labs(title="parity 4 - all")







#####################################################################################
# Other approaches


#------------------------------------------------------------------------------------
# output

#lines(predict(fit_pos4)~df_pos$DIM,col="purple")

# test with pos2:
# a tibble with mean, median, sd:
out_pos2 %>%
  as_tibble() %>% 
  tidy(out_pos2)
dplyr::n_distinct(df2_pos$BES_ID)


#---------------------------------------------------------------------------
# identify outliers with rstatix:
library(rstatix)


extreme_pos2 <- out_pos2 %>%
  identify_outliers("a")

extreme_pos2k <- out_pos2 %>%
  identify_outliers("b")

out_pos2 %>%
  identify_outliers("b")

# replace row number with the BES_ID
# create a list with the TRUE to extreme outlier
# save as data: a data frame for each parameter in each output = 24 data frames
# seperate into lower or higher outlier -> 48 data frames
# create curves for all these...

ggplot(out_pos3) +
  aes(x = "", y = b) +
  geom_boxplot(fill = "#0c4c8a") +
  theme_minimal()

out2_pos_extreme <- is_extreme(out_pos2$k)

out_pos2 %>%
  as_tibble() %>% 
  is_extreme(out_pos2$k)


#------------------------------------------------------------------------
boxplot.stats(out_pos2$a)$out

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


#-----------------------------------------------------------------
# manuel filter for removing outliers

filter_4pos <- out_4pos %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)



#-----------------------------------------------------------------
# testing outlier detection:

out_k <- which(out_4pos$k %in% c(out_k1))
out_k



#-------------------------------------------------------
# individual herd plots of the crazy parameters:


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



