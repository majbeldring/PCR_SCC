

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Model assessment


# Packages and settings ----------------------------------------

library(tidyverse)
library(ggpubr) # p values in boxplot
library(gridExtra) # gridarrange
ggplot2::theme_set(ggplot2::theme_bw())  # globally sets ggplot2 theme to theme_bw


Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


# Loading data and preparing data -------------------------------------

#load("K:/paperI/major4/VI_parameters_PCR_major4.RData")
load("K:/paperI/major4/VII_deliverance.RData") # only model output and grouped dataset. Notice: 4 -> >3 not changed in figures

# check output generated here:
# load("K:/paperI/major4/backup_old_pre9999/nlme_wilmink4.RData") 




# diagnostic plot ----------------------------------------------------

# check: https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/model-checking.html

# QQ

# Q-Q plots take your sample data, sort it in ascending order, 
# and then plot them versus quantiles calculated from a theoretical distribution. 
# The number of quantiles is selected to match the size of your sample data
# Use the qq quantile plot...


# NEGATIVE: #F8766D
# POSITIVE #00BFC4
qq_pos1 <- ggqqplot(residuals(nlme_pos1), shape=1, 
              xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", title = "Parity 1, positive", 
              font.x = c(14), font.y = c(14), 
              color = "#00BFC4")
qq_pos1$layers[[2]]$aes_params$colour <- "black" 
#qq_neg3

qq_all <- ggarrange(qq_neg1, qq_pos1, qq_neg2, qq_pos2, qq_neg3, qq_pos3, qq_neg4, qq_pos4, 
                    ncol=2, nrow=4, 
                    common.legend = TRUE, legend="right")
qq_all

ggsave("C:/Users/zjt234/PhD/PaperI_PCR/002_pvm_submission_revised/qq_all.tiff", width = 40, height = 40, units = "cm", dpi=300)
#ggsave("C:/Users/zjt234/PhD/PaperI_PCR_Wilmink/final_figures/qq_all.pdf")
#ggsave("C:/Users/zjt234/PhD/PaperI_PCR_Wilmink/final_figures/qq_neg2.tiff", width = 40, height = 40, units = "cm", dpi=300)
#ggsave("C:/Users/zjt234/PhD/PaperI_PCR_Wilmink/final_figures/qq_neg2.pdf")

# basic R (ylab is cut off...)
# qqnorm(resid(nlme_neg2), pch=1, cex.lab = 1.5, main=NULL) # quantiles
# qqline(resid(nlme_neg2), col = "red", lwd = 2)
# 
# tiff(file="C:/Users/zjt234/PhD/PaperI_PCR_Wilmink/final_figures/qq_neg2.tiff",
#      width=40, height=40, units="in", res=100)
# qqnorm(resid(nlme_neg2), pch=1, cex.lab = 2, main=NULL) # quantiles
# qqline(resid(nlme_neg2), col = "red", lwd = 2)
# dev.off()





# RESIDUALS Figure 5 paper (sample is used) -----------------------------------


# plot of fitted vs residuals: Remake as a ggplot or with neater 
plot(fitted(nlme_neg2), resid(nlme_neg2), pch=20); abline(h=0, col="red", lwd = 2) # Works fine.

# test with ggplot... Only a part of the data
library(broom.mixed)
broom::augment(nlme_pos4, data = df4_pos) -> 
  aug_nlme_pos4


# gg residuals plot with a larger sample of data:
res_pos3 <- ggplot(aug_nlme_pos3 %>% sample_n(25e3), aes(.fitted, .resid)) + 
  ggpubr::theme_classic2() + 
  #geom_point(size=1) +
  # geom_point(alpha = 1/10) +
  geom_point(shape=1, size=2, fill='white') +
  geom_hline(yintercept = 0, color = 'red', size=1) +
  #geom_smooth(se=FALSE, color = "red") +
  labs(x = "Fitted", y = "Residuals", title = "Parity 3, positive") +
  ylim(-6, 6) +
  xlim(2.5, 7) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
        #,legend.position = "top"
  ) +
  NULL

resid_sample <- ggarrange(res_neg1, res_pos1, res_neg2, res_pos2, res_neg3, res_pos3, res_neg4, res_pos4,  
                       ncol=2, nrow=4, 
                       common.legend = TRUE, legend="right")

resid_sample
ggsave("C:/Users/zjt234/PhD/PaperI_PCR/002_pvm_submission_revised/Residuals.tiff", width = 40, height = 40, units = "cm", dpi=300)



# Residuals plot, all data.
resid6 <- ggplot(aug_nlme_pos4, aes(.fitted, .resid)) + 
  ggpubr::theme_classic2() + 
  #geom_point(size=1) +
  # geom_point(alpha = 1/10) +
  geom_point(shape=1, size=1, fill='white') +
  geom_hline(yintercept = 0, color = "red", size=1) +
  # geom_smooth(se=FALSE, color = "red") +
  labs(x = "Fitted", y = "Residuals", title = "Parity 4, positive") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)
        #,legend.position = "top"
  ) +
  NULL

resid_all <- ggarrange(resid1, resid4, resid2, resid5, resid3, resid6, 
                    ncol=2, nrow=3, 
                    common.legend = TRUE, legend="right")


ggsave("C:/Users/zjt234/PhD/PaperI_PCR_Wilmink/final_figures/resid_all.tiff", width = 40, height = 40, units = "cm", dpi=300)
















# NOT used diagnostic - or NOT working.......................

library(car)
qqPlot(neg2_residuals)
qqplot(neg2_residuals)



qqnorm(nlme_neg2, datax=T)
qqline(x, datax=T, col="blue", lwd=2)



library(ggfortify) # fortify is not working with nlme. So can't make it into a dataframe..
ggplot(nlme_neg2, aes(x = .fitted, y = .resid)) + geom_point() # invalid, not a dataframe
# basic summary
summary(nlme_neg2)
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(nlme_neg2)  # Plot the model information

# ggplot, neg2. Not plotting the quantiles....
library(data.table)
neg2_residuals <- data.table(neg2_residuals = residuals(nlme_neg2))
neg2_residuals[, logSCC := neg2_residuals/sd(neg2_residuals)]
ggqqplot(data = neg2_residuals,
         x = "logSCC", color = "black", 
         #ggtheme = theme_bw()
)

shapiro.test(resid(nlme_neg2)[1:1000]) # not useful for big data

# QQ plot:
qqnorm(nlme_neg2, ~ ranef(., standard = TRUE)) # shows each parameter. Don't use
qqnorm(nlme_neg2, ~ resid(., type = "p")) 

plot(nlme_neg1) # not useful. Data is 

# test if relevant to compare with nls:
nls_neg1 <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                         data = df1_neg,
                                         lower=c(a=0, b=0, k=-5, d=0),
                                         upper=c(a=9, b=1.5, k=0, d=5),
                                         start_lower = c(a=0, b=0, k=-5, d=0),
                                         start_upper = c(a=8, b=1, k=-0.01, d=4),
                                         iter = 500,
                                         supp_errors = "Y")

nls_out_neg1 <- coef(nls_neg1) %>% 
  as_tibble() 

qqnorm(nls_neg1, ~ resid(., type = "p")) 



# CI ---------------------------------------------------------------

library(tidyverse)
plyr::ddply(all_out_nlme, c("PARITY", "PCR"), summarise,
               N    = length(d),
               mean = mean(d),
               median = median(d),
               sd   = sd(d),
               se   = sd / sqrt(N),
               lower = mean - qt(1 - (0.05 / 2), N - 1) * se,
               upper = mean + qt(1 - (0.05 / 2), N - 1) * se)



# other methods: 

# library(Rmisc) # list as upper, mean, lower
# CI(nlme_out_neg2$a, ci=0.95) # 3.667126 3.636400 3.605674

# #manual calculation of 95% CI:
# sample.mean <- mean(nlme_out_neg1$a)
# sample.n <- length(nlme_out_neg1$a)
# sample.sd <- sd(nlme_out_neg1$a)
# sample.se <- sample.sd/sqrt(sample.n)

# degrees.freedom = sample.n - 1
# t.score = qt(p=0.05/2, df=degrees.freedom,lower.tail=F)

# margin.error <- t.score * sample.se
# lower.bound <- sample.mean - margin.error
# upper.bound <- sample.mean + margin.error

# library(MASS)
# confint(nlme_neg1, level = 0.95) # needs MASS to be installed



# Outliers parameters -------------------------------------------

box_param <- ggplot(nlme_out_pos4, aes(class, hwy)) + 
  geom_boxplot(outlier.colour = "red", outlier.shape = 1)




# SAVE -------------------------------------------------

save.image("K:/paperI/major4/resubmit.RData")
