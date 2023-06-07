

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Model assessment


# PCR project
# Residual plots

# Figure 5


# settings and data ------------------------------
library("tidyverse")
library("nlme")
library(gridExtra) # gridarrange
library(broom.mixed) # augment
library(ggpubr) # gridarrange 
ggplot2::theme_set(ggplot2::theme_bw())  # globally sets ggplot2 theme to theme_bw

(load("K:/paperI/major4/Ouput_nlme.RData"))




# RESIDUALS Figure 5 paper (sample is used) -----------------------------------

cbbPalette <- c("#56B4E9", "#E69F00", "#009E73","#CC79A7", "#F0E442", "#0072B2", "#D55E00")


# test with ggplot... Only a part of the data
broom::augment(nlme_pos1, data = df1_pos) -> 
  aug_nlme_pos1

# gg residuals plot with a larger sample of data:
ggplot(aug_nlme_pos1 %>% sample_n(25e3), aes(.fitted, .resid)) + 
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

ggplot(aug_nlme_pos4, aes(.fitted, .resid)) + 
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



# NEW RESIDUAL DENSITY PLOT ----------------------------


tibble(Fitted = fitted(nlme_pos1), Residual = resid(nlme_pos1)) |>
  mutate(Fitted = round(Fitted, 1), Residual = round(Residual, 1)) |>
  count(Fitted, Residual) |>
  ggplot(aes(x=Fitted, y=Residual, z=n)) +
  geom_contour()

res1neg <- tibble(Fitted = fitted(nlme_neg1), Residual = resid(nlme_neg1)) |>
  ggplot(aes(x=Fitted, y=Residual)) +
  geom_density2d(color = "#56B4E9") +
  stat_smooth(se=FALSE, color = "#56B4E9")  + 
  labs(x = "Fitted", y = "Residuals", title = "Parity 1, negative") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))

resid_all <- ggarrange(res1neg, res1pos, res2neg, res2pos, res3neg, res3pos, res4neg, res4pos,  
                          ncol=2, nrow=4, 
                          common.legend = TRUE, legend="right")

resid_all
ggsave("C:/Users/zjt234/PhD/PaperI_PCR/003_Animals_Submission/new_resid.tiff", width = 40, height = 20, units = "cm", dpi=300)



# Boxplot. add real data
data_with_residual <- df1_neg |> mutate(Residual = resid(nlme_neg1))
ggplot(data_with_residual, aes(x=DIM, y=Residual, group=DIM)) +
  geom_boxplot()


