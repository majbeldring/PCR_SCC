

# Matt Denwood / Maj Beldring

## Objective:  determine se/sp of SCC compared to PCR based on model output


# Data and packages ------------------------------------------

library("tidyverse")
library("nlme")
library("pbapply")

(load("K:/paperI/major4/Ouput_nlme.RData"))


# 1. Extract prediction and residual for each model ---------------
expand_grid(DIM = 1:305, PCR = c("negative","positive"), Parity=1:4) |>
  group_by(PCR, Parity) |>
  group_split() |>
  pblapply(function(x){
    mod <- get(str_c("nlme_", str_sub(x$PCR[1], 1, 3), x$Parity[1]))
    x |>
      mutate(Prediction = predict(mod, newdata=x, level=0)) |>
      mutate(ResidualSD = summary(mod)$sigma, Nobs = length(predict(mod)))
  }) |>
  bind_rows() ->
  predictions

# 2. Calculate NPV and PPV based on normal distribution and different SCC thresholds:
predictions |>
  expand_grid(SCC_Threshold = c(150,200,400)) |>
  mutate(Parameter = case_when(
    PCR=="negative" ~ "NPV",
    PCR=="positive" ~ "PPV"
  )) |>
  mutate(Value = case_when(
    Parameter=="NPV" ~ pnorm(log(SCC_Threshold), Prediction, ResidualSD),
    Parameter=="PPV" ~ pnorm(log(SCC_Threshold), Prediction, ResidualSD, lower.tail=FALSE)
  )) |>
  identity() ->
  pred_vals

# 3. Convert to Se and Sp:
pred_vals |>
  select(-PCR, -Prediction, -ResidualSD) |>
  pivot_wider(names_from="Parameter", values_from=c("Value","Nobs")) |>
  rename(NPV=Value_NPV, PPV=Value_PPV, Nneg=Nobs_NPV, Npos=Nobs_PPV) |>
  mutate(ObsPrev = Npos/(Nneg+Npos)) |>
  mutate(Se = (PPV*ObsPrev) / ((PPV*ObsPrev) + (1-NPV)*(1-ObsPrev))) |>
  mutate(Sp = (NPV*(1-ObsPrev)) / ((1-PPV)*ObsPrev + NPV*(1-ObsPrev))) |>
  pivot_longer(cols=c(NPV,PPV,Se,Sp), names_to="Parameter", values_to="Value") ->
  estimates

# 4. plot
ggplot(estimates |> filter(Parameter %in% c("Se","Sp")), aes(x=DIM, y=Value*100, col=factor(SCC_Threshold))) +
  geom_line(size = rel(0.5)) +
  facet_grid(str_c("Parity: ", Parity) ~ Parameter) +
  theme_light() +
  theme(legend.pos="bottom") +
  labs(y=NULL, color="SCC Threshold:  ") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

ggsave("C:/Users/zjt234/PhD/PaperI_PCR/002_pvm_submission_revised/Se_Sp.tiff", width = 40, height = 40, units = "cm", dpi=300)


### Notes

# We can go from PPV, NPV, obs prev to Se, Sp using the same formulae as to go from Se, Sp, true prev to PPV, NPV:
TP <- runif(1)
FP <- runif(1)
TN <- runif(1)
FN <- runif(1)

PPV = TP/(TP+FP)
NPV = TN/(FN+TN)
Se = TP/(TP+FN)
Sp = TN/(TN+FP)

ObsPrev = (TP+FP) / (TP+FP+TN+FN)
TruePrev = (TP+FN) / (TP+FP+TN+FN)

PPV; (Se*TruePrev) / ((Se*TruePrev) + (1-Sp)*(1-TruePrev))
NPV; (Sp*(1-TruePrev)) / ((1-Se)*TruePrev + Sp*(1-TruePrev))

Se; (PPV*ObsPrev) / ((PPV*ObsPrev) + (1-NPV)*(1-ObsPrev))
Sp; (NPV*(1-ObsPrev)) / ((1-PPV)*ObsPrev + NPV*(1-ObsPrev))

