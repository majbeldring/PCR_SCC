## Updated 2023-03-13 to correct a mistake in the previous code/graph (ignoring herd)

library("tidyverse")
library("nlme")
library("pbapply")

(load("K:/paperI/major4/Ouput_nlme.RData"))

## Merging conditional modes:
expand_grid(PCR = c("negative","positive"), Parity=1:4) |>
  group_by(PCR, Parity) |>
  group_split() |>
  pblapply(function(x){
    m <- get(str_c("nlme_", str_sub(x$PCR[1], 1, 3), x$Parity[1]))
    coef(m) |>
      as.data.frame() |>
      rownames_to_column("BES_ID") |>
      bind_cols(x)
  }) |>
  bind_rows() |>
  pivot_longer(cols=a:d, names_to="Parameter", values_to="Estimate") ->
estimates_by_herd

# Now you can look for changes within herd eg:
estimates_by_herd |>
  pivot_wider(names_from="PCR", values_from=Estimate) |>
  mutate(Ratio = positive / negative) |>
  ggplot(aes(x=Ratio, col=factor(Parity))) +
  stat_ecdf() +
  facet_wrap(~Parameter, scales="free_x") +
  xlab("Within-herd ratio of PCR positive to negative") +
  geom_vline(xintercept=1) +
  theme_light()


## Objective:  determine se/sp of SCC compared to PCR based on model output

## 0. Find herds with both negative and positive PCR tests per parity
expand_grid(PCR = c("negative","positive"), Parity=1:4) |>
  group_by(Parity) |>
  group_split() |>
  pblapply(function(x){
    m1 <- get(str_c("nlme_", str_sub(x$PCR[1], 1, 3), x$Parity[1]))
    m2 <- get(str_c("nlme_", str_sub(x$PCR[2], 1, 3), x$Parity[2]))
    inner_join(
      tibble(BES_ID = dimnames(ranef(m1))[[1]]),
      tibble(BES_ID = dimnames(ranef(m2))[[1]]),
      by="BES_ID"
    ) |>
      mutate(Parity = x$Parity[1])
  }) |>
  bind_rows() ->
  herds

# 1. Extract prediction and residual for each model and herd
expand_grid(DIM = 1:400, PCR = c("negative","positive"), Parity=1:4) |>
  full_join(herds, by="Parity", multiple="all") |>
  group_by(PCR, Parity) |>
  group_split() |>
  pblapply(function(x){
    mod <- get(str_c("nlme_", str_sub(x$PCR[1], 1, 3), x$Parity[1]))
    x |>
      mutate(Prediction = predict(mod, newdata=x, level=1)) |>
      mutate(ResidualSD = summary(mod)$sigma) |>
      ## NOTE: I don't have the full data here, so am using number of SCC obs
      ## by herd as a proxy for number of animals per herd
      inner_join(
        mod$groups |> count(BES_ID, name="Nobs"),
        by="BES_ID"
      )
  }) |>
  bind_rows() ->
  predictions
stopifnot(all(!is.na(predictions)))

# 2. Calculate NPV and PPV based on normal distribution and different SCC thresholds:
predictions |>
  expand_grid(SCC_Threshold = c(100,150,200,400)) |>
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
  all_estimates

# Summarise herds:
all_estimates |>
  mutate(SCC_Threshold = factor(SCC_Threshold)) |>
  group_by(Parameter, Parity, DIM, SCC_Threshold) |>
  summarise(Mean = mean(Value), LCI = quantile(Value, 0.025), UCI = quantile(Value, 0.975), .groups="drop") ->
  estimates

# 4. plot of Se/Sp (probability of PCR conditional on SCC)
ggplot(estimates |> filter(Parameter %in% c("Se","Sp")), aes(x=DIM, y=Mean*100, ymin=LCI*100, ymax=UCI*100, col=SCC_Threshold, fill=SCC_Threshold)) +
  # CI are optional:
  #geom_ribbon(alpha=0.25, size=0) +
  geom_line() +
  facet_grid(str_c("Parity: ", Parity) ~ Parameter) +
  theme_light() +
  theme(legend.pos="bottom") +
  labs(y=NULL, color="SCC Threshold:  ", fill="SCC Threshold:  ")


# 5. plot of PPV/NPV (probability of SCC conditional on PCR)
ggplot(estimates |> filter(Parameter %in% c("PPV","NPV")), aes(x=DIM, y=Mean*100, ymin=LCI*100, ymax=UCI*100, col=SCC_Threshold, fill=SCC_Threshold)) +
  # CI are optional:
  geom_ribbon(alpha=0.25, size=0) +
  geom_line() +
  facet_grid(str_c("Parity: ", Parity) ~ Parameter) +
  theme_light() +
  theme(legend.pos="bottom") +
  labs(y=NULL, color="SCC Threshold:  ", fill="SCC Threshold:  ")
  #scale_color_brewer()


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
