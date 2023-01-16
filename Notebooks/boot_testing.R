

# Maj Beldring Henningsen, majbh@sund.ku.dk

# NLS Curves using Bootstrapping


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
library(nlshelper) # for tidy(fit)
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
df <- df_model %>% 
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)

# eco, Holstein, Parity 2:
df <- df %>% 
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 2) %>%
  filter(HERDTYPE == 0) %>%
  filter(DIM < 306) %>%
  dplyr::select(BES_ID, DYR_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, 
                MILK, RES_MAJOR)

dplyr::n_distinct(df$BES_ID) # 160
dplyr::n_distinct(df$DYR_ID) # 7327

#----------------------------------------------------------------------
# RES major POS vs NEg

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
plot(boo_pos, type = "boxplot", ask = FALSE)

# Matrix with the bootstrapped parameter estimates
Theta_pos <- boo_pos$coefboot

# Model: logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])

# Points where to evaluate the model
x_eval <- seq(min(df_pos$DIM), max(df_pos$DIM), length.out = 100)

# Matrix with the predictions
Pred_pos <- apply(Theta_pos, 1, function(theta) fun(x_eval, theta))

# Pack the estimates for plotting
Estims_pos <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_pos, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

p_pos <- ggplot(data = Estims_pos, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(0.5), colour = "darkgreen") + 
  ylim(3.5, 5.3) +
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

# Defining x axis
x_eval <- seq(min(df_neg$DIM), max(df_neg$DIM), length.out = 100)

# prediction mtrix
Pred_neg <- apply(Theta_neg, 1, function(theta) fun(x_eval, theta))

# Pack the estimates for plotting
Estims_neg <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_neg, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

p_neg <- ggplot(data = Estims_neg, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  #stat_summary(geom="ribbon", fun.data=mean_cl_normal, width=0.1, conf.int=0.95, fill="lightblue") +
  #geom_ribbon(aes(ymin = ci_lower_est, ymax = ci_upper_est),    # shadowing cnf intervals 
  #fill = "grey") + 
  geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(0.5), colour = "darkred") + 
  ylim(3.5, 5.3) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink for: PCR NEG, parity 2, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_neg

grid.arrange(p_pos, p_neg, ncol=2, nrow=1)



