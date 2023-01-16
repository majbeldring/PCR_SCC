

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Wilmink curves with nls for: PCR TEST POS OR NEG

# Script includes:
## 1) preparedata: df_pos & df_neg for each parity -> 6 datasets in total to run
## 2) create 3 curves (1 with two curves for each parity) using boot (not herd level)
## 3) create curves on herd level
## 4) ggpairs


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

load("M:/PCR_data/PCR_merge2.RData") 
rm(df_pcr, df_curve); gc() 

# pre-preparing data;remove dates and DYR_ID:
df_model <- df_model %>% 
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)


# descriptives on over all data:
dplyr::n_distinct(df_model$BES_ID) # 3474 herds
dplyr::n_distinct(df_model$DYR_ID) # 997.784


summary(df_model)

#---------------------------------------------------------------------
# choose analysis focus:

# Choose: 1) Test type 2) breed, 3) parity, 4) herdtype:
df <- df_model %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 4) %>%
  filter(HERDTYPE == 1) %>%
  filter(DIM < 306) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, 
                MILK, RES_MAJOR)

dplyr::n_distinct(df_pos$BES_ID) # 1141
dplyr::n_distinct(df$DYR_ID) # 3075

#----------------------------------------------------------------------
# specify data and bootstrapping
## NOTE: MUST RUN FOR EACH PARITY. 
## REMEMBER RUN THE CURVE LINE PLOT (line 194) AFTER EACH PRITY GROUP HAS BEEN FITTED

### RES MAJOR POS
# only treated with max 200 observations:
df_pos <- df %>%
  filter(RES_MAJOR == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
df_pos$BES_ID <- factor(df_pos$BES_ID) # keep only used levels by resetting the variable


# Boot fitting treated:
nls_pos <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_pos, 
             start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
             control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                            printEval = FALSE, warnOnly = TRUE), 
             algorithm = "port", trace = FALSE)

boo_pos <- nlsBoot(nls_pos, niter = 400)

# boxplot of parameters
p_boo_pos <- plot(boo_pos, type = "boxplot", ask = FALSE)

# Matrix with the bootstrapped parameter estimates
Theta_pos <- boo_pos$coefboot

# Model: logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])

# X axis (DIM alredy set to 0-305 days)
x_eval <- seq(min(df_pos$DIM), max(df_pos$DIM), length.out = 100)

# Matrix with the predictions
Pred_pos <- apply(Theta_pos, 1, function(theta) fun(x_eval, theta))

# Packing estimates pre plotting
Estims_pos <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_pos, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

# ggplot with wilmink curve created with output parameters
p_pos <- ggplot(data = Estims_pos, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  #geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 6.0) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink for: PCR POS, parity 2, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_pos




#--------------------------------------------------------------------------
### RES MAJOR NEG

# only treated with max 200 observations:
df_neg <- df %>%
  filter(RES_MAJOR == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
df_neg$BES_ID <- factor(df_neg$BES_ID) # keep only used levels by resetting the variable


# Boot fitting treated:
nls_neg <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_neg, 
               start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
               control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                              printEval = FALSE, warnOnly = TRUE), 
               algorithm = "port", trace = FALSE)

boo_neg <- nlsBoot(nls_neg, niter = 400)

# boxplot of parameters
plot(boo_neg, type = "boxplot", ask = FALSE)

# Matrix with the bootstrapped parameter estimates
Theta_neg <- boo_neg$coefboot

# Model: logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])

# X axis (DIM alredy set to 0-305 days)
x_eval <- seq(min(df_neg$DIM), max(df_neg$DIM), length.out = 100)

# prediction mtrix
Pred_neg <- apply(Theta_neg, 1, function(theta) fun(x_eval, theta))

# Packing estimates pre plotting
Estims_neg <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_neg, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

# ggplot with wilmink curve created with output parameters
p_neg <- ggplot(data = Estims_neg, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  #geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(1.5), colour = "darkred") + 
  ylim(3.5, 6.0) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink for: PCR NEG, parity 2, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_neg


#------------------------------------------------------
# combine curves, histograms etc:

# curves: create for each parity
grid.arrange(p_pos, p_neg, ncol=2, nrow=1)

# histograms of parameters: create for each parity


# Curves
p_4 <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 5.75) +
  theme_bw() + labs(title = "Parity 4: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_4

p_4_2 <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 6.0) +
  theme_bw() + labs(title = "Parity 2: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_4_2

p_4_3 <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 6.5) +
  theme_bw() + labs(title = "Parity 2: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_4_3

p_4_4 <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 7.0) +
  theme_bw() + labs(title = "Parity 2: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_4_4



grid.arrange(p_2_2, p_3_2, p_4_2, ncol=3, nrow=1)



#---------------------------------------------------------------------------------
# Wilmink NLS on herd level with equation 8 græsbøll:

# STARTLIST defined and equation 8
st <- list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6)
f_nls <- logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | BES_ID


# conventionel, holstein, parity 4, Major POS
fit_pos4 <- nlsList(f_nls, df_pos, start = sapply(st, mean),
               control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                              printEval = FALSE, warnOnly = TRUE))

out_4pos <- coef(fit_pos4) %>%
  drop_na()

plot_nls(fit_pos4)

tidy(fit_pos4)
# ... and confidence intervals
tidy(fit_pos4, conf.int=TRUE)

# identify outliers in parameter output:
boxplot.stats(out_4pos$a)$out
boxplot.stats(out_4pos$b)$out
boxplot.stats(out_4pos$k)$out
boxplot.stats(out_4pos$d)$out

out_a1 <- boxplot.stats(out_4pos$a)$out
out_a <- which(out_4pos$a %in% c(out_a1))
out_a

out_b1 <- boxplot.stats(out_4pos$b)$out
out_b <- which(out_4pos$b %in% c(out_b1))
out_b

out_k1 <- boxplot.stats(out_4pos$k)$out
out_k <- which(out_4pos$k %in% c(out_k1))
out_k

out_d1 <- boxplot.stats(out_4pos$d)$out
out_d <- which(out_4pos$d %in% c(out_d1))
out_d

out_4pos[out_k, ]


# manuel filter for outliers
filter_4pos <- out_4pos %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)


curve_pos4 <- out_4pos %>%
  summarise_at(c("a", "b", "k", "d"), mean, na.rm = TRUE)

plot(df_pos$DIM,df_pos$logSCC)
lines(df_pos$DIM,predict(fit_pos4))

#-------------------------------------------------------------
# plotting curves

fit_pos4 <- lm(logSCC ~ 4.627 + 0.00366 * DIM + exp(-(exp(-2.00955)) * DIM)*2.3088, data = df_pos)
plot(test_curve)

tidymodel_data <- df_pos %>%
  ungroup() %>% #was rowwise_df, so needed to ungroup(), https://stackoverflow.com/a/33292159/4927395
  select(c("BES_ID", "DIM", "logSCC")) %>%  #select only essential columns. Data to nest within CowID, because multiple observations for each CowID at each DIM
  nest(yields = c(logSCC, DIM)) %>%
  mutate(model = map(
    #         yields, ~ lm(MilkYieldPerDay ~ DIM, 
    yields, ~ nls(formula = logSCC ~ 4.627 + 0.00366 * DIM + exp(-(exp(-2.00955)) * DIM)*2.3088, data = .x)))



p <- ggplot(data = df_pos, 
            aes(x = DIM, y = logSCC)) + 
  geom_smooth()+
  ylab("logSCC")+
  xlab("DIM")+
  theme_bw()+  
  theme(legend.position = "none")
p


ggplot(df_pos, aes(DIM, logSCC)) +
  geom_smooth(method = lm, 
              formula = logSCC ~ 4.627 + 0.00366 * DIM + exp(-(exp(-2.00955)) * DIM)*2.3088)

#---------------------------------------------------------------------------------
# histograms based on wilmink nls on herd level

# Filter (based on output): Data prep so this can be avoided
out_m2 <- out_m2 %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)

m2_a <- qplot(out_m2$a,
              geom="histogram",
              binwidth = 0.1,  
              #main = "a: con, parity 2, holstein", 
              xlab = "Parity 2, a values",  
              fill=I("white"), col=I("green4"))

#----------------------------------------------------------------------------------
# ggpairs

out_4pos %>% 
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5) %>%
  ggpairs()



#----------------------------------------------------------------------
# constructing curves on a herd level







library(nlshelper)
plot_nls(out_4pos, df_pos)



# nls in ggplot - ONLY with test data when facet_wrap is applied
p_nls <- ggplot(df_pos, aes(x = DIM, y = logSCC)) + 
  ggtitle("nls, logSCC~DIM") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_nls

# nls in ggplot:
p_neg <- ggplot(df_neg, aes(x = DIM, y = logSCC)) + 
  ggtitle("parity 3+ neg, logSCC~DIM") +
  stat_smooth(method = 'nls', formula = y ~ a + b * x + exp(-(exp(k)) * x)*d, se = FALSE, 
              method.args = list(start=c(a = 4, b = 0.003, k = -2, d = 2.6), 
                                 control=nls.control(maxiter=200)), colour = "black")+
  geom_smooth() 
p_neg





#------------------------

