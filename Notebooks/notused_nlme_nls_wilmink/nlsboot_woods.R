

# Maj Beldring Henningsen, majbh@sund.ku.dk

# nlsboot function for recreating curves on national level


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
# Loading data 

load("M:/PCR_data/woods_MILK_out.RData") # from PCR_data script

#--------------------------------------------------------------------------
# nlsBoot

# remember to save plot before re-running for another Parity

## POS
nls_pos <- nls(formula = MILK ~ a * (DIM ^ b) * exp(-(exp(k)) * DIM), data = df4_pos, 
             start = list(a = 20, b = 0.2, k = -5.6),
             control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                            printEval = FALSE, warnOnly = TRUE), 
             algorithm = "port", trace = FALSE)

boo_pos <- nlsBoot(nls_pos, niter = 400)
p_boo_pos <- plot(boo_pos, type = "boxplot", ask = FALSE) # boxplot of parameters
Theta_pos <- boo_pos$coefboot # Matrix with the bootstrapped parameter estimates

# Model: MILK ~  a * (DIM ^ b) * exp(-(exp(k)) * DIM)
fun <- function(DIM, theta) 
  theta["a"] * (DIM ^ theta["b"]) * exp(-(exp(theta["k"])) * DIM)

x_eval <- seq(min(df4_pos$DIM), max(df4_pos$DIM), length.out = 100) # X axis (DIM alredy set to 0-305 days)
Pred_pos <- apply(Theta_pos, 1, function(theta) fun(x_eval, theta)) # Matrix with the predictions

# Packing estimates pre plotting
Estims_pos <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_pos, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)


## NEG
nls_neg <- nls(formula = MILK ~ a * (DIM ^ b) * exp(-(exp(k)) * DIM), data = df4_neg, 
               start = list(a = 20, b = 0.2, k = -5.6),
               control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                              printEval = FALSE, warnOnly = TRUE), 
               algorithm = "port", trace = FALSE)

boo_neg <- nlsBoot(nls_neg, niter = 400)
p_boo_neg <- plot(boo_neg, type = "boxplot", ask = FALSE) # boxplot of parameters
Theta_neg <- boo_neg$coefboot # Matrix with the bootstrapped parameter estimates

# Model: MILK ~  a * (DIM ^ b) * exp(-(exp(k)) * DIM)
fun <- function(DIM, theta) 
  theta["a"] * (DIM ^ theta["b"]) * exp(-(exp(theta["k"])) * DIM)

x_eval <- seq(min(df4_neg$DIM), max(df4_neg$DIM), length.out = 100) # X axis (DIM alredy set to 0-305 days)
Pred_neg <- apply(Theta_neg, 1, function(theta) fun(x_eval, theta)) # Matrix with the predictions

# Packing estimates pre plotting
Estims_neg <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_neg, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)


#------------------------------------------------------
# curves - national level, created with bootstrapping parameters


p2_milk <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(0, 75) +
  theme_bw() + labs(title = "Boot nls, Parity 2: Red = NEG, Green = POS", x = "DIM", y = "MILK")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)

p3_milk <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(0, 75) +
  theme_bw() + labs(title = "Boot nls, Parity 2: Red = NEG, Green = POS", x = "DIM", y = "MILK")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)

p4_milk <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(0, 75) +
  theme_bw() + labs(title = "Boot nls, Parity 2: Red = NEG, Green = POS", x = "DIM", y = "MILK")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)




grid.arrange(p2_milk, p3_milk, p4_milk, ncol=3, nrow=1)



