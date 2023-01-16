

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
# Loading preppet data

load("M:/PCR_data/woods_SCC_out.RData") # from PCR_data script


#-------------------------------------------------------------------------
# Define wilmink function (national level)

# Model: logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) 
  theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])


#--------------------------------------------------------------------------
# nlsBoot


# Boot fitting treated:
nls_pos <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_pos, 
             start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
             control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                            printEval = FALSE, warnOnly = TRUE), 
             algorithm = "port", trace = FALSE)

boo_pos <- nlsBoot(nls_pos, niter = 400) # runs 10-30 min.

p_boo_pos <- plot(boo_pos, type = "boxplot", ask = FALSE) # boxplot of parameters
Theta_pos <- boo_pos$coefboot # Matrix with the bootstrapped parameter estimates

x_eval <- seq(min(df_pos$DIM), max(df_pos$DIM), length.out = 100) # X axis (DIM alredy set to 0-305 days)
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



# NEG
nls_neg <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_neg, 
               start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
               control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                              printEval = FALSE, warnOnly = TRUE), 
               algorithm = "port", trace = FALSE)

boo_neg <- nlsBoot(nls_neg, niter = 400)

plot(boo_neg, type = "boxplot", ask = FALSE) # boxplot of parameters
Theta_neg <- boo_neg$coefboot # Matrix with the bootstrapped parameter estimates
x_eval <- seq(min(df_neg$DIM), max(df_neg$DIM), length.out = 100) # X axis (DIM alredy set to 0-305 days)
Pred_neg <- apply(Theta_neg, 1, function(theta) fun(x_eval, theta)) # prediction matrix

# Packing estimates pre plotting
Estims_neg <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_neg, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)


#-------------------------------------------------------------------------
# plot single parity (national level)


# ggplot with wilmink curve created with output parameters
p_pos <- ggplot(data = Estims_pos, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  #geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 6.0) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink for: PCR POS, parity 2, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_pos


# ggplot with wilmink curve created with output parameters
p_neg <- ggplot(data = Estims_neg, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  #geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(1.5), colour = "darkred") + 
  ylim(3.5, 6.0) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink for: PCR NEG, parity 2, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_neg

# curves: create for each parity
grid.arrange(p_pos, p_neg, ncol=2, top="Parity 2 - SCC")


#------------------------------------------------------
# Plot all parities


# Curves
p2_SCC <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 5.75) +
  theme_bw() + labs(title = "Parity 2: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)

p3_SCC <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 6.0) +
  theme_bw() + labs(title = "Parity 3: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)

p4_SCC <- ggplot() + 
  geom_line(data = Estims_neg, aes(x = x, y = median_est), size = rel(1.5), colour = "darkred") + 
  geom_line(data = Estims_pos, aes(x = x, y = median_est), size = rel(1.5), colour = "darkgreen") + 
  ylim(3.5, 6.5) +
  theme_bw() + labs(title = "Parity 2: Red = NEG, Green = POS", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)




grid.arrange(p2_SCC, p3_SCC, p4_SCC, ncol=3, nrow=1)



