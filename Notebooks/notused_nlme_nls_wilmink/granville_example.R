
# NLS, NLME, BRMS

# amazing script copied from:
# https://www.granvillematheson.com/post/nonlinear-modelling-using-nls-nlme-and-brms/





#### Libraries and settings ####

#remotes::install_github("mathesong/kinfitr")
library(tidyverse)
library(kinfitr)
library(nls.multstart)
library(nlme)
library(brms)
library(hrbrthemes)
library(broom)
library(viridis)
colourcodes <- c("#d4a665", "#d27fff", "#7fd9ff")
colourpal <- c(NLS="#d4a665", NLME="#d27fff", MCMC="#7fd9ff")

theme_set(hrbrthemes::theme_ipsum_rc())

#### Data from kinfitr package ####

data(pbr28)
pbr28$jsondata[[1]]$Metabolite

xtract_metabolite <- function(jsondata) {
  suppressMessages(
    parentFraction <- as_tibble(jsondata$Metabolite$Data$Values,
                                .name_repair = "unique")
  )
  names(parentFraction) = jsondata$Metabolite$Data$Labels
  parentFraction$Time = with(parentFraction,
                             sampleStartTime + 0.5*sampleDuration)
  parentFraction <- select(parentFraction, Time, parentFraction)
  return(parentFraction)
}

pfdata <- pbr28 %>% 
  mutate(pf = map(jsondata, extract_metabolite)) %>% 
  select(PET, Subjname, PETNo, Genotype, pf)


#### NLS ####

# nonlinear least squares:
# NLS makes use of gradient descent to find the set of parameters, 
# which maximise the likelihood of the data under them.
# Issues with local minimal - multistart avoid this 
# multistart uses Levenberg–Marquardt algorithm

# visualize data:

pfdat <- pfdata$pf[[1]]

ggplot(pfdat, aes(x=Time, y=parentFraction)) +
  geom_point(size=3) +
  ylim(c(0, 1))

# fitting

hillfunc <- function(Time, a, b, c) {
  1 - ( ( (1-a) * Time^b) / ( c + (Time)^b ) )
} 

hill_nls_fit <- nls.multstart::nls_multstart(parentFraction ~ hillfunc(Time, a, b, c),
                                             data = pfdat,
                                             lower=c(a=0, b=1, c=0),
                                             upper=c(a=1, b=Inf, c=1e12),
                                             start_lower = c(a=0, b=1, c=0),
                                             start_upper = c(a=0.5, b=500, c=1e12),
                                             iter = 500,
                                             supp_errors = "Y")

summary(hill_nls_fit)

plot_nls <- function(nls_object, data) {
  predframe <- tibble(Time=seq(from=min(data$Time), to=max(data$Time), 
                               length.out = 1024)) %>%
    mutate(ParentFrac = predict(nls_object, newdata = list(Time=.$Time)))
  ggplot(data, aes(x=Time, y=parentFraction)) +
    geom_point(size=3) +
    geom_line(data = predframe, aes(x=Time, y=ParentFrac))
}

plot_nls(hill_nls_fit, pfdat)


# re-fitting

hillfunc <- function(Time, a, b, c) {
  1 - ( ( (1-a) * Time^b) / ( 10^c + (Time)^b ) )
} 

hill_nls_fit <- nls.multstart::nls_multstart(parentFraction ~ hillfunc(Time, a, b, c),
                                             data = pfdat,
                                             lower=c(a=0, b=1, c=0),
                                             upper=c(a=1, b=Inf, c=12),
                                             start_lower = c(a=0, b=1, c=0),
                                             start_upper = c(a=0.5, b=500, c=12),
                                             iter = 500,
                                             supp_errors = "Y")

summary(hill_nls_fit)

plot_nls(hill_nls_fit, pfdat)


# fitting with minpck (parametres most be better for this):

hill_nls_fit2 <- minpack.lm::nlsLM(parentFraction ~ hillfunc(Time, a, b, c),
                                   data = pfdat,
                                   lower=c(a=0, b=1, c=0),
                                   upper=c(a=1, b=Inf, c=12),
                                   start=c(a=0.5, b=5, c=5))


# comparing multstart and minpack:
coef(hill_nls_fit) # multstart
coef(hill_nls_fit2) # not multstart 


#### NLS: Fitting the curves ####

hill_nls_fit_func <- function(pf_df) {
  nls.multstart::nls_multstart(parentFraction ~ hillfunc(Time, a, b, c),
                               data = pf_df,
                               lower=c(a=0, b=1, c=0),
                               upper=c(a=1, b=Inf, c=12),
                               start_lower = c(a=0, b=1, c=0),
                               start_upper = c(a=0.5, b=500, c=12),
                               iter=500,
                               supp_errors = "Y")
} 

pfdata <- pfdata %>% 
  mutate(hill_nls_fit = map(pf, ~hill_nls_fit_func(.x)))


plot_nls( pfdata$hill_nls_fit[[3]], pfdata$pf[[3]])
plot_nls( pfdata$hill_nls_fit[[8]], pfdata$pf[[8]])
plot_nls( pfdata$hill_nls_fit[[12]], pfdata$pf[[12]])


# checking parameter distribution


hill_nls_outcomes <- pfdata %>% 
  mutate(outpars = map(hill_nls_fit, ~broom::tidy(.x))) %>% 
  select(-pf, -hill_nls_fit) %>% 
  unnest(cols="outpars")

ggplot(hill_nls_outcomes, aes(x=estimate, colour=term, fill=term)) +
  geom_density(alpha=0.5, fill=colourcodes[1], colour=colourcodes[1]) +
  facet_wrap(~term, scales = "free")

hill_nls_outcomes_summary <- hill_nls_outcomes %>%
  group_by(term) %>% 
  summarise(mean = mean(estimate), 
            median = median(estimate),
            sd = sd(estimate)) %>% 
  ungroup()

knitr::kable(hill_nls_outcomes_summary, digits = 3)

#### NLS: checing the fits ####


hill_nls_plots <- pfdata %>% 
  select(PET, pf) %>%
  unnest(pf)

hill_predtimes <- tidyr::crossing(PET=pfdata$PET, 
                                  Time=seq(min(hill_nls_plots$Time),
                                           max(hill_nls_plots$Time),
                                           length.out=128))

hill_nlspreds <- hill_predtimes %>% 
  group_by(PET) %>% 
  nest(preds = Time) %>% 
  left_join(select(pfdata, PET, hill_nls_fit)) %>% 
  mutate(preds = map2(preds, hill_nls_fit, ~broom::augment(.y, newdata=.x))) %>% 
  select(-hill_nls_fit) %>% 
  ungroup() %>% 
  unnest(cols=preds)


ggplot(hill_nls_plots, aes(x=Time, y=parentFraction)) +
  geom_point() +
  geom_line(data=hill_nlspreds, aes(y=.fitted), colour=colourcodes[1], size=0.7) +
  facet_wrap(~PET, ncol=4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### NLME : Fitting using frequentist multilevel modelling ####

# nlme:
# multilevel modelling approach, 
# still using maximum likelihood


# prepare data in long unnested form:
pf_modeldata <- pfdata %>% 
  select(PET:pf) %>% 
  unnest(cols="pf")

head(pf_modeldata, n = 12)


# nlme fit:
hill_nlme_fit <- nlme(parentFraction ~ hillfunc(Time, a, b, c), 
                      data = pf_modeldata,
                      fixed=a + b + c ~ 1, 
                      random = a + b + c ~ 1, 
                      groups = ~ PET, 
                      start = hill_nls_outcomes_summary$mean,
                      verbose = F)

summary(hill_nlme_fit)

# parameters:

nlme_coef = as_tibble(coef(hill_nlme_fit), rownames = 'PET')
nlme_coef

#### nlme - checking the fit ####

hill_nlmepreds <- hill_predtimes %>% 
  mutate(.fitted=predict(hill_nlme_fit, newdata=hill_predtimes))

ggplot(pf_modeldata, aes(x=Time, y=parentFraction)) +
  geom_point() +
  geom_line(data=hill_nlmepreds, aes(y=.fitted), colour=colourcodes[2], size=0.7) +
  facet_wrap(~PET, ncol=4)


# compare to nls output:

nlme_coef_tidy <- nlme_coef %>% 
  gather(Parameter, Estimate, -PET) %>% 
  mutate(Model = "NLME")

nls_coef_tidy <- hill_nls_outcomes %>% 
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
hillprior <- c(
  set_prior("normal(0.2, 0.1)", nlpar = "a", lb=0, ub=1),
  set_prior("normal(2, 1)", nlpar = "b", lb=1),
  set_prior("normal(7, 3)", nlpar = "c", lb=0), 
  set_prior("normal(0.05, 0.2)", class="sigma"))

hill_bayes_fit_formula <- bf(parentFraction ~ 1 - ( ( (1-a) * Time^b) / 
                                                      ( 10^c + (Time)^b ) ),
                             # Nonlinear variables
                             a + b + c ~ 1,
                             # Nonlinear fit
                             nl = TRUE)

hill_bayes_fit <- brm(
  hill_bayes_fit_formula,
  family=gaussian(), 
  data = pfdat,
  prior = hillprior )



#### brms - checking the fit ####

summary(hill_bayes_fit)
plot(hill_bayes_fit)
pairs(hill_bayes_fit)

# degree of correlation:

predtimes <- unique(hill_predtimes$Time)

hill_bayes_fitted <- fitted(hill_bayes_fit, 
                            newdata=list(Time = predtimes)) %>% 
  as_tibble()

hill_bayes_pred <- predict(hill_bayes_fit,
                           newdata=list(Time = predtimes)) %>%  
  as_tibble()


hill_bayes_ribbons <- tibble(
  Time = predtimes,
  parentFraction=hill_bayes_fitted$Estimate,
  Estimate = hill_bayes_fitted$Estimate,
  pred_lower = hill_bayes_pred$Q2.5,
  pred_upper = hill_bayes_pred$Q97.5,
  fitted_lower = hill_bayes_fitted$Q2.5,
  fitted_upper = hill_bayes_fitted$Q97.5)

ggplot(pfdat, aes(x=Time, y=parentFraction)) +
  geom_point(size=3) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.2, fill=colourcodes[3]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.5, fill=colourcodes[3]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Estimate), colour=colourcodes[3], 
            size=1)

     
summary(hill_nls_fit)
summary(hill_bayes_fit)


# showing hill_stan code:

hillstan <- "
  real hill_stan(real Time, real a, real b, real c) {
  
    real pred;
    
    pred = 1 - ( ( (1-a) * Time^b) / ( 10^c + (Time)^b ) );
    
    return(pred);
  }
"

get_prior(bf(parentFraction ~ hill_stan(Time, a, b, c),
             # Nonlinear variables
             a + b + c ~ 1 + (1|k|PET),
             # Nonlinear fit
             nl = TRUE), data=pf_modeldata)


hillprior_multilevel <- c(
  set_prior("normal(0.2, 0.1)", nlpar = "a", lb=0, ub=1),
  set_prior("normal(2, 1)", nlpar = "b", lb=1),
  set_prior("normal(7, 3)", nlpar = "c", lb=0), 
  set_prior("normal(0.05, 0.02)", class="sigma"),
  set_prior("normal(0.03, 0.02)", class="sd", nlpar="a"),
  set_prior("normal(0.3, 0.1)", class="sd", nlpar="b"),
  set_prior("normal(0.7, 0.2)", class="sd", nlpar="c"),
  set_prior("lkj(2)", class = "cor"))

hill_multilevelbayes_formula <- bf(parentFraction ~ hill_stan(Time, a, b, c),
                                   # Nonlinear variables
                                   a + b + c ~ 1 + (1|k|PET),
                                   # Nonlinear fit
                                   nl = TRUE)


# check STAN code before fitting (showing what we DON'T have to write when using STAN)

make_stancode(hill_multilevelbayes_formula,
              family=gaussian(), 
              data = pf_modeldata,
              prior = hillprior_multilevel,
              stanvars = stanvar(scode = hillstan, 
                                 block="functions"))

# fit STAN

hill_multilevelbayes_fit <- brm(
  hill_multilevelbayes_formula,
  family=gaussian(), 
  data = pf_modeldata,
  prior = hillprior_multilevel,
  stanvars = stanvar(scode = hillstan, 
                     block="functions"), 
  cores = 4)


# check the fit:

summary(hill_multilevelbayes_fit)
plot(hill_multilevelbayes_fit, ask = FALSE)


# individual fits:

expose_functions(hill_multilevelbayes_fit, vectorize = TRUE)

hill_mlbayes_pred <- predict(hill_multilevelbayes_fit, 
                             newdata=hill_predtimes) %>% 
  as_tibble()
hill_mlbayes_fitted <- fitted(hill_multilevelbayes_fit, 
                              newdata=hill_predtimes) %>% 
  as_tibble()

hill_mlbayes_ribbons <- tibble(
  PET = hill_predtimes$PET,
  Time = hill_predtimes$Time,
  parentFraction = hill_mlbayes_fitted$Estimate,
  pred_lower = hill_mlbayes_pred$Q2.5,
  pred_upper = hill_mlbayes_pred$Q97.5,
  fitted_lower = hill_mlbayes_fitted$Q2.5,
  fitted_upper = hill_mlbayes_fitted$Q97.5)

ggplot(pf_modeldata, aes(x=Time, y=parentFraction)) +
  geom_point() +
  geom_line(data=hill_mlbayes_ribbons, aes(y=parentFraction), colour=colourcodes[3], 
            size=1)  +
  geom_ribbon(data=hill_mlbayes_ribbons, alpha=0.2, aes(ymin=pred_lower, 
                                                        ymax=pred_upper), 
              fill=colourcodes[3]) +
  geom_ribbon(data=hill_mlbayes_ribbons, alpha=0.5, aes(ymin=fitted_lower,
                                                        ymax=fitted_upper), 
              fill=colourcodes[3]) +
  facet_wrap(~PET, ncol=4)




#### AVERAGE curves ####
probnames <- c(20, 50, 80, 95)
probs <- c(40, 60, 25, 75, 10, 90, 2.5, 97.5)/100

probtitles <- probs[order(probs)]*100
probtitles <- paste("Q", probtitles, sep="")


hill_mlbayes_avgpred <- predict(hill_multilevelbayes_fit, 
                                newdata=list(Time=hill_predtimes$Time),
                                re_formula = NA,
                                probs=probs) %>% 
  as_tibble() %>% 
  mutate(Curve = "Prediction Intervals",
         Effects = "Fixed")

hill_mlbayes_avgfitted <- fitted(hill_multilevelbayes_fit, 
                                 newdata=list(Time=hill_predtimes$Time),
                                 re_formula = NA,
                                 probs=probs) %>% 
  as_tibble() %>% 
  mutate(Curve = "Credible Intervals",
         Effects = "Fixed")

hill_mlbayes_avgfitted_ns <- fitted(hill_multilevelbayes_fit, 
                                    newdata=list(Time=hill_predtimes$Time, 
                                                 PET=rep("new", nrow(hill_predtimes))),
                                    probs=probs, allow_new_levels=TRUE) %>% 
  as_tibble() %>% 
  mutate(Curve = "Credible Intervals",
         Effects = "Fixed + Random")

hill_mlbayes_avgpred_ns <- predict(hill_multilevelbayes_fit, 
                                   newdata=list(Time=hill_predtimes$Time, 
                                                PET=rep("new", nrow(hill_predtimes))),
                                   probs=probs, allow_new_levels=TRUE) %>% 
  as_tibble() %>% 
  mutate(Curve = "Prediction Intervals",
         Effects = "Fixed + Random")

hill_mlbayes_aribbons <- bind_rows(hill_mlbayes_avgpred, 
                                   hill_mlbayes_avgfitted,
                                   hill_mlbayes_avgpred_ns,
                                   hill_mlbayes_avgfitted_ns) %>% 
  mutate(Time = rep(hill_predtimes$Time, 4))


avg_pal <- viridis::plasma(n=4)
names(avg_pal) <- paste(probnames, "%", sep="")

ggplot(hill_mlbayes_aribbons, aes(x=Time, y=Estimate)) +
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

hill_mlbayes_arrays <- coef(hill_multilevelbayes_fit)

hill_mlbayes_outcomes <- rbind(a=hill_mlbayes_arrays$PET[, , 1],
                               b=hill_mlbayes_arrays$PET[, , 2],
                               c=hill_mlbayes_arrays$PET[, , 3]) %>% 
  as_tibble(rownames='PET') %>% 
  mutate(Model = 'MCMC',
         Parameter = rep(c('a', 'b', 'c'), each=nrow(pfdata))) %>% 
  select(PET, Parameter, Estimate, Model)

model_outcome_comparison <- bind_rows(nls_nlme_comparison, hill_mlbayes_outcomes)

ggplot(model_outcome_comparison, aes(x=Estimate, colour=Model, fill=Model)) +
  geom_density(alpha=0.3) +
  scale_fill_manual(values=colourpal) +
  scale_colour_manual(values=colourpal) +
  facet_wrap(~Parameter, scales="free")


preddata <- tibble(
  PET = hill_nlspreds$PET,
  Time = hill_nlspreds$Time,
  NLS = hill_nlspreds$.fitted,
  NLME = hill_nlmepreds$.fitted,
  MCMC = hill_mlbayes_fitted$Estimate
) %>% 
  gather(Model, parentFraction, -PET, -Time) %>% 
  mutate(Model = fct_inorder(Model)) %>% 
  arrange(PET, Time, Model)


ggplot(pf_modeldata, aes(x=Time, y=parentFraction)) +
  geom_point() +
  geom_line(data=preddata, aes(y=parentFraction, colour=Model), 
            size=0.7) +
  facet_wrap(~PET, ncol=2) +
  scale_colour_manual(values=colourpal)

