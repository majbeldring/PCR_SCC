
# NLS, NLME, BRMS


#### Libraries and settings ####

library(tidyverse)
library(nls.multstart)
library(nlme)
library(brms)
library(hrbrthemes)
library(broom)
library(viridis)
colourcodes <- c("#d4a665", "#d27fff", "#7fd9ff")
colourpal <- c(NLS="#d4a665", NLME="#d27fff", MCMC="#7fd9ff")

theme_set(hrbrthemes::theme_ipsum_rc())

#### Loading data and preparing data ####

load("M:/PCR_data/PCR_merge.RData") 
rm(df_pcr, df_curve); gc() 

# one herd (holstein, conventional, Parity 2, herd no. 5051112)
df <- df_model %>% 
  filter(DIM < 306) %>%
  filter(DIM > 5) %>%
  filter(PCR_TEST == 1) %>%
  filter(BREED == 1) %>%
  filter(HERDTYPE == 1) %>%
  filter(PARITY == 2) %>%
  filter(BES_ID == 5051112) %>%
  #ungroup() %>%
  group_by(DYR_ID) %>%
  mutate(count = n()) %>%
  filter(count > 8) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC, RES_MAJOR)

df_pos <- df %>%
  filter(RES_MAJOR == 1) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC)

df_neg <- df %>%
  filter(RES_MAJOR == 0) %>%
  ungroup() %>%
  dplyr::select(DYR_ID, DIM, logSCC)

df$DYR_ID <- factor(df$DYR_ID) # removed unused levels by resetting
df$RES_MAJOR <- factor(df$RES_MAJOR) 
df_pos$DYR_ID <- factor(df_pos$DYR_ID)
df_neg$DYR_ID <- factor(df_neg$DYR_ID)


#### NLS ####

# nonlinear least squares:
# NLS makes use of gradient descent to find the set of parameters, 
# which maximise the likelihood of the data under them.
# Issues with local minimal - multistart avoid this 
# multistart uses Levenberg-Marquardt algorithm

# visualize data:

ggplot(df, aes(x=DIM, y=logSCC)) +
  geom_point(size=1) +
  ylim(c(1, 9))

# fitting

wilminkfunc <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}

wilmink_nls_fit <- nls.multstart::nls_multstart(logSCC ~ wilminkfunc(DIM, a, b, k, d),
                                             data = df,
                                             lower=c(a=0, b=0, k=-5, d=0),
                                             upper=c(a=9, b=1.5, k=0, d=5),
                                             start_lower = c(a=0, b=0, k=-5, d=0),
                                             start_upper = c(a=8, b=1, k=-0.01, d=4),
                                             iter = 500,
                                             supp_errors = "Y")

summary(wilmink_nls_fit)

plot_nls <- function(nls_object, data) {
  predframe <- tibble(DIM=seq(from=min(data$DIM), to=max(data$DIM), 
                               length.out = 1024)) %>%
    mutate(SCC = predict(nls_object, newdata = list(DIM=.$DIM)))
  ggplot(data, aes(x=DIM, y=logSCC)) +
    geom_point(size=1) +
    geom_line(data = predframe, aes(x=DIM, y=SCC))
}

plot_nls(wilmink_nls_fit, df)




# fitting with minpck (parametres most be better for this):

wilmink_nls_fit2 <- minpack.lm::nlsLM(logSCC ~ wilminkfunc(DIM, a, b, k, d),
                                   data = df,
                                   lower=c(a=0, b=0, k=-5, d=0),
                                   upper=c(a=9, b=1.5, k=0, d=5),
                                   start=c(a=4, b=0.004, k=-2.3, d=1.8))


# comparing multstart and minpack:
coef(wilmink_nls_fit) # multstart
coef(wilmink_nls_fit2) # not multstart 


#### NLS: Fitting the curves ####

wilmink_nls_fit_func <- function(pf_df) {
  nls.multstart::nls_multstart(logSCC ~ wilminkfunc(DIM, a, b, k, d),
                               data = pf_df,
                               lower=c(a=0, b=0, k=-5, d=0),
                               upper=c(a=9, b=1.5, k=0, d=5),
                               start_lower = c(a=0, b=0, k=-5, d=0),
                               start_upper = c(a=8, b=1, k=-0.01, d=4),
                               iter=500,
                               supp_errors = "Y")
} 

res_data <- df %>% 
  nest(df_res = c(DYR_ID, DIM, logSCC)) %>% 
  mutate(wilmink_nls_fit = map(df_res ~wilmink_nls_fit_func(.x)))


plot_nls(dfa$wilmink_nls_fit[[3]], dfa$pf[[3]])
plot_nls(dfa$wilmink_nls_fit[[8]], dfa$pf[[8]])
plot_nls( dfa$wilmink_nls_fit[[12]], dfa$pf[[12]])


# checking parameter distribution


wilmink_nls_outcomes <- dfa %>% 
  mutate(outpars = map(wilmink_nls_fit, ~broom::tidy(.x))) %>% 
  select(-pf, -wilmink_nls_fit) %>% 
  unnest(cols="outpars")

ggplot(wilmink_nls_outcomes, aes(x=estimate, colour=term, fill=term)) +
  geom_density(alpha=0.5, fill=colourcodes[1], colour=colourcodes[1]) +
  facet_wrap(~term, scales = "free")

wilmink_nls_outcomes_summary <- wilmink_nls_outcomes %>%
  group_by(term) %>% 
  summarise(mean = mean(estimate), 
            median = median(estimate),
            sd = sd(estimate)) %>% 
  ungroup()

knitr::kable(wilmink_nls_outcomes_summary, digits = 3)

#### NLS: checing the fits ####


wilmink_nls_plots <- dfa %>% 
  select(PET, pf) %>%
  unnest(pf)

wilmink_predDIMs <- tidyr::crossing(PET=dfa$PET, 
                                  DIM=seq(min(wilmink_nls_plots$DIM),
                                           max(wilmink_nls_plots$DIM),
                                           length.out=128))

wilmink_nlspreds <- wilmink_predDIMs %>% 
  group_by(PET) %>% 
  nest(preds = DIM) %>% 
  left_join(select(dfa, PET, wilmink_nls_fit)) %>% 
  mutate(preds = map2(preds, wilmink_nls_fit, ~broom::augment(.y, newdata=.x))) %>% 
  select(-wilmink_nls_fit) %>% 
  ungroup() %>% 
  unnest(cols=preds)


ggplot(wilmink_nls_plots, aes(x=DIM, y=logSCC)) +
  geom_point() +
  geom_line(data=wilmink_nlspreds, aes(y=.fitted), colour=colourcodes[1], size=0.7) +
  facet_wrap(~PET, ncol=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### NLME : Fitting using frequentist multilevel modelling ####

# nlme:
# multilevel modelling approach, 
# still using maximum likelihood


# prepare data in long unnested form:
pf_modeldata <- dfa %>% 
  select(PET:pf) %>% 
  unnest(cols="pf")

head(pf_modeldata, n = 12)


# nlme fit:
wilmink_nlme_fit <- nlme(logSCC ~ wilminkfunc(DIM, a, b, k, d), 
                      data = pf_modeldata,
                      fixed=a + b + c ~ 1, 
                      random = a + b + c ~ 1, 
                      groups = ~ PET, 
                      start = wilmink_nls_outcomes_summary$mean,
                      verbose = F)

summary(wilmink_nlme_fit)

# parameters:

nlme_coef = as_tibble(coef(wilmink_nlme_fit), rownames = 'PET')
nlme_coef

#### nlme - checking the fit ####

wilmink_nlmepreds <- wilmink_predDIMs %>% 
  mutate(.fitted=predict(wilmink_nlme_fit, newdata=wilmink_predDIMs))

ggplot(pf_modeldata, aes(x=DIM, y=logSCC)) +
  geom_point() +
  geom_line(data=wilmink_nlmepreds, aes(y=.fitted), colour=colourcodes[2], size=0.7) +
  facet_wrap(~PET, ncol=4)


# compare to nls output:

nlme_coef_tidy <- nlme_coef %>% 
  gather(Parameter, Estimate, -PET) %>% 
  mutate(Model = "NLME")

nls_coef_tidy <- wilmink_nls_outcomes %>% 
  select(PET, Parameter=term, Estimate=estimate) %>% 
  mutate(Model = "NLS")

nls_nlme_comparison <- full_join(nls_coef_tidy, nlme_coef_tidy)

ggplot(nls_nlme_comparison, aes(x=Estimate, colour=Model, fill=Model)) +
  geom_density(alpha=0.3) +
  scale_colour_manual(values=colourpal) +
  scale_fill_manual(values=colourpal) +
  facet_wrap(~Parameter, scales="free")



#### brms - bayesian multilevel ####

# brms:
# using Markov Chain Monte Carlo (MCMC) instead of maximum likelihood


# brms - modelling a single curve:
wilminkprior <- c(
  set_prior("normal(0.2, 0.1)", nlpar = "a", lb=0, ub=1),
  set_prior("normal(2, 1)", nlpar = "b", lb=1),
  set_prior("normal(7, 3)", nlpar = "c", lb=0), 
  set_prior("normal(0.05, 0.2)", class="sigma"))

wilmink_bayes_fit_formula <- bf(logSCC ~ 1 - ( ( (1-a) * DIM^b) / 
                                                      ( 10^c + (DIM)^b ) ),
                             # Nonlinear variables
                             a + b + c ~ 1,
                             # Nonlinear fit
                             nl = TRUE)

wilmink_bayes_fit <- brm(
  wilmink_bayes_fit_formula,
  family=gaussian(), 
  data = df,
  prior = wilminkprior )



#### brms - checking the fit ####

summary(wilmink_bayes_fit)
plot(wilmink_bayes_fit)
pairs(wilmink_bayes_fit)

# degree of correlation:

predDIMs <- unique(wilmink_predDIMs$DIM)

wilmink_bayes_fitted <- fitted(wilmink_bayes_fit, 
                            newdata=list(DIM = predDIMs)) %>% 
  as_tibble()

wilmink_bayes_pred <- predict(wilmink_bayes_fit,
                           newdata=list(DIM = predDIMs)) %>%  
  as_tibble()


wilmink_bayes_ribbons <- tibble(
  DIM = predDIMs,
  logSCC=wilmink_bayes_fitted$Estimate,
  Estimate = wilmink_bayes_fitted$Estimate,
  pred_lower = wilmink_bayes_pred$Q2.5,
  pred_upper = wilmink_bayes_pred$Q97.5,
  fitted_lower = wilmink_bayes_fitted$Q2.5,
  fitted_upper = wilmink_bayes_fitted$Q97.5)

ggplot(df, aes(x=DIM, y=logSCC)) +
  geom_point(size=3) +
  geom_ribbon(data=wilmink_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.2, fill=colourcodes[3]) +
  geom_ribbon(data=wilmink_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.5, fill=colourcodes[3]) +
  geom_line(data=wilmink_bayes_ribbons, aes(y=Estimate), colour=colourcodes[3], 
            size=1)

     
summary(wilmink_nls_fit)
summary(wilmink_bayes_fit)


# showing wilmink_stan code:

wilminkstan <- "
  real wilmink_stan(real DIM, real a, real b, real c) {
  
    real pred;
    
    pred = 1 - ( ( (1-a) * DIM^b) / ( 10^c + (DIM)^b ) );
    
    return(pred);
  }
"

get_prior(bf(logSCC ~ wilmink_stan(DIM, a, b, k, d),
             # Nonlinear variables
             a + b + c ~ 1 + (1|k|PET),
             # Nonlinear fit
             nl = TRUE), data=pf_modeldata)


wilminkprior_multilevel <- c(
  set_prior("normal(0.2, 0.1)", nlpar = "a", lb=0, ub=1),
  set_prior("normal(2, 1)", nlpar = "b", lb=1),
  set_prior("normal(7, 3)", nlpar = "c", lb=0), 
  set_prior("normal(0.05, 0.02)", class="sigma"),
  set_prior("normal(0.03, 0.02)", class="sd", nlpar="a"),
  set_prior("normal(0.3, 0.1)", class="sd", nlpar="b"),
  set_prior("normal(0.7, 0.2)", class="sd", nlpar="c"),
  set_prior("lkj(2)", class = "cor"))

wilmink_multilevelbayes_formula <- bf(logSCC ~ wilmink_stan(DIM, a, b, k, d),
                                   # Nonlinear variables
                                   a + b + c ~ 1 + (1|k|PET),
                                   # Nonlinear fit
                                   nl = TRUE)


# check STAN code before fitting (showing what we DON'T have to write when using STAN)

make_stancode(wilmink_multilevelbayes_formula,
              family=gaussian(), 
              data = pf_modeldata,
              prior = wilminkprior_multilevel,
              stanvars = stanvar(scode = wilminkstan, 
                                 block="functions"))

# fit STAN

wilmink_multilevelbayes_fit <- brm(
  wilmink_multilevelbayes_formula,
  family=gaussian(), 
  data = pf_modeldata,
  prior = wilminkprior_multilevel,
  stanvars = stanvar(scode = wilminkstan, 
                     block="functions"), 
  cores = 4)


# check the fit:

summary(wilmink_multilevelbayes_fit)
plot(wilmink_multilevelbayes_fit, ask = FALSE)


# individual fits:

expose_functions(wilmink_multilevelbayes_fit, vectorize = TRUE)

wilmink_mlbayes_pred <- predict(wilmink_multilevelbayes_fit, 
                             newdata=wilmink_predDIMs) %>% 
  as_tibble()
wilmink_mlbayes_fitted <- fitted(wilmink_multilevelbayes_fit, 
                              newdata=wilmink_predDIMs) %>% 
  as_tibble()

wilmink_mlbayes_ribbons <- tibble(
  PET = wilmink_predDIMs$PET,
  DIM = wilmink_predDIMs$DIM,
  logSCC = wilmink_mlbayes_fitted$Estimate,
  pred_lower = wilmink_mlbayes_pred$Q2.5,
  pred_upper = wilmink_mlbayes_pred$Q97.5,
  fitted_lower = wilmink_mlbayes_fitted$Q2.5,
  fitted_upper = wilmink_mlbayes_fitted$Q97.5)

ggplot(pf_modeldata, aes(x=DIM, y=logSCC)) +
  geom_point() +
  geom_line(data=wilmink_mlbayes_ribbons, aes(y=logSCC), colour=colourcodes[3], 
            size=1)  +
  geom_ribbon(data=wilmink_mlbayes_ribbons, alpha=0.2, aes(ymin=pred_lower, 
                                                        ymax=pred_upper), 
              fill=colourcodes[3]) +
  geom_ribbon(data=wilmink_mlbayes_ribbons, alpha=0.5, aes(ymin=fitted_lower,
                                                        ymax=fitted_upper), 
              fill=colourcodes[3]) +
  facet_wrap(~PET, ncol=4)




#### AVERAGE curves ####
probnames <- c(20, 50, 80, 95)
probs <- c(40, 60, 25, 75, 10, 90, 2.5, 97.5)/100

probtitles <- probs[order(probs)]*100
probtitles <- paste("Q", probtitles, sep="")


wilmink_mlbayes_avgpred <- predict(wilmink_multilevelbayes_fit, 
                                newdata=list(DIM=wilmink_predDIMs$DIM),
                                re_formula = NA,
                                probs=probs) %>% 
  as_tibble() %>% 
  mutate(Curve = "Prediction Intervals",
         Effects = "Fixed")

wilmink_mlbayes_avgfitted <- fitted(wilmink_multilevelbayes_fit, 
                                 newdata=list(DIM=wilmink_predDIMs$DIM),
                                 re_formula = NA,
                                 probs=probs) %>% 
  as_tibble() %>% 
  mutate(Curve = "Credible Intervals",
         Effects = "Fixed")

wilmink_mlbayes_avgfitted_ns <- fitted(wilmink_multilevelbayes_fit, 
                                    newdata=list(DIM=wilmink_predDIMs$DIM, 
                                                 PET=rep("new", nrow(wilmink_predDIMs))),
                                    probs=probs, allow_new_levels=TRUE) %>% 
  as_tibble() %>% 
  mutate(Curve = "Credible Intervals",
         Effects = "Fixed + Random")

wilmink_mlbayes_avgpred_ns <- predict(wilmink_multilevelbayes_fit, 
                                   newdata=list(DIM=wilmink_predDIMs$DIM, 
                                                PET=rep("new", nrow(wilmink_predDIMs))),
                                   probs=probs, allow_new_levels=TRUE) %>% 
  as_tibble() %>% 
  mutate(Curve = "Prediction Intervals",
         Effects = "Fixed + Random")

wilmink_mlbayes_aribbons <- bind_rows(wilmink_mlbayes_avgpred, 
                                   wilmink_mlbayes_avgfitted,
                                   wilmink_mlbayes_avgpred_ns,
                                   wilmink_mlbayes_avgfitted_ns) %>% 
  mutate(DIM = rep(wilmink_predDIMs$DIM, 4))


avg_pal <- viridis::plasma(n=4)
names(avg_pal) <- paste(probnames, "%", sep="")

ggplot(wilmink_mlbayes_aribbons, aes(x=DIM, y=Estimate)) +
  geom_ribbon(aes_string(ymin=probtitles[1], ymax=probtitles[2]), 
              fill=avg_pal[1], alpha=0.6) +
  geom_ribbon(aes_string(ymin=probtitles[7], ymax=probtitles[8]), 
              fill=avg_pal[1], alpha=0.6) +
  geom_ribbon(aes_string(ymin=probtitles[2], ymax=probtitles[3]), 
              fill=avg_pal[2], alpha=0.6) +
  geom_ribbon(aes_string(ymin=probtitles[6], ymax=probtitles[7]), 
              fill=avg_pal[2], alpha=0.6) +
  geom_ribbon(aes_string(ymin=probtitles[3], ymax=probtitles[4]), 
              fill=avg_pal[3], alpha=0.6) +
  geom_ribbon(aes_string(ymin=probtitles[5], ymax=probtitles[6]), 
              fill=avg_pal[3], alpha=0.6) +
  geom_ribbon(aes_string(ymin=probtitles[4], ymax=probtitles[5]), 
              fill=avg_pal[4], alpha=0.6) +
  facet_wrap(Effects~Curve, ncol=2)



#### COMPARING ALL the fits ####

wilmink_mlbayes_arrays <- coef(wilmink_multilevelbayes_fit)

wilmink_mlbayes_outcomes <- rbind(a=wilmink_mlbayes_arrays$PET[, , 1],
                               b=wilmink_mlbayes_arrays$PET[, , 2],
                               c=wilmink_mlbayes_arrays$PET[, , 3]) %>% 
  as_tibble(rownames='PET') %>% 
  mutate(Model = 'MCMC',
         Parameter = rep(c('a', 'b', 'c'), each=nrow(dfa))) %>% 
  select(PET, Parameter, Estimate, Model)

model_outcome_comparison <- bind_rows(nls_nlme_comparison, wilmink_mlbayes_outcomes)

ggplot(model_outcome_comparison, aes(x=Estimate, colour=Model, fill=Model)) +
  geom_density(alpha=0.3) +
  scale_fill_manual(values=colourpal) +
  scale_colour_manual(values=colourpal) +
  facet_wrap(~Parameter, scales="free")


preddata <- tibble(
  PET = wilmink_nlspreds$PET,
  DIM = wilmink_nlspreds$DIM,
  NLS = wilmink_nlspreds$.fitted,
  NLME = wilmink_nlmepreds$.fitted,
  MCMC = wilmink_mlbayes_fitted$Estimate
) %>% 
  gather(Model, logSCC, -PET, -DIM) %>% 
  mutate(Model = fct_inorder(Model)) %>% 
  arrange(PET, DIM, Model)


ggplot(pf_modeldata, aes(x=DIM, y=logSCC)) +
  geom_point() +
  geom_line(data=preddata, aes(y=logSCC, colour=Model), 
            size=0.7) +
  facet_wrap(~PET, ncol=2) +
  scale_colour_manual(values=colourpal)

