

# Maj Beldring Henningsen, majbh@sund.ku.dk

# wilmink for SCC curves
# nlme with nls starthelp

# tidy up:
## name each section: section 1 ####
## Mossas code: apply for: output, mean and median, plots
## for those possible: use nls output directly as start for nlme



# Packages and settings -------------------------------------

library(tidyverse)
library(ggpubr) # p values in boxplot
library(gridExtra) # gridarrange
ggplot2::theme_set(ggplot2::theme_bw())  # globally sets ggplot2 theme to theme_bw
library(GGally) # for ggpairs
library(nls.multstart)
library(nlme) # for nlslist

Sys.setlocale("LC_ALL","English") # date formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters



# Loading data and preparing data -------------------------------------

load("M:/PCR_data/curve6_SCC.RData") # from PCR crypted container

# load("M:/PCR_data/curve6_nlme_wilmink.RData")
save.image("M:/PCR_data/curve6_nlme_wilmink_copy.RData")


# wilmink function ----------------------------------------

f_wilmink <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}


# nls model - output for starting nlme ----------------------

# nls multistart (not grouped herd level):
# repeat for all 6 diff. datasets
nls_pos4 <- nls.multstart::nls_multstart(logSCC ~ f_wilmink(DIM, a, b, k, d),
                                         data = df4_pos,
                                         lower=c(a=0, b=0, k=-5, d=0),
                                         upper=c(a=9, b=1.5, k=0, d=5),
                                         start_lower = c(a=0, b=0, k=-5, d=0),
                                         start_upper = c(a=8, b=1, k=-0.01, d=4),
                                         iter = 500,
                                         supp_errors = "Y")

# repeat for all groups:
start_pos2 <- coef(nls_pos4) %>% 
  as_tibble() 


# nlme wilmink  ---------------------------------------

# repeat for all 6 different data sets
nlme_pos4 <- nlme(logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                  data=df4_pos,
                  fixed=a+b+k+d~1,
                  random=a+b+k+d~1,
                  groups=~BES_ID,
                  start=c(a = 4.5, b = 0.004, k = -1.5, d = 3.2),
                  # start = nlme_start_pos2$value,
                  na.action=na.exclude,
                  control = list(maxIter = 1200, msMaxIter = 1200))

out_nlme_pos4 <- coef(nlme_pos4)


# output stats
coef(nlme_pos4) %>% 
  as_tibble() %>% 
  summarise(across(everything(), mean))
  #summarise(across(everything(), median))
  #summarise(across(everything(), 
  #                 list(
  #                   "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1)))))



# intercept and quantiles  -------------------------------------

# all_quan_sum_curves: listing parity, PCR, intercept and DIM quantiles:
list(
  list(PCR ="NEG", PARITY = 2, output = out_nlme_neg2),
  list(PCR ="NEG", PARITY = 3, output = out_nlme_neg3),
  list(PCR ="NEG", PARITY = 4, output = out_nlme_neg4),
  list(PCR ="POS", PARITY = 2, output = out_nlme_pos2),
  list(PCR ="POS", PARITY = 3, output = out_nlme_pos3),
  list(PCR ="POS", PARITY = 4, output = out_nlme_pos4)
) %>%
  enframe() %>%
  select(-name) %>% 
  unnest_wider(value) %>% 
  # mutate(mean_output = output %>% map(colMeans)) %>% 
  # mutate(mean_output = output %>% map(~ apply(.x, MARGIN = 2, function(x) median(x)))) %>%
  # mutate(mean_output = output %>% map(~ lapply(.x, function(x) median(x)))) %>%
  mutate(mean_output = output %>% map(~ lapply(.x, function(x) 
    quantile(x, probs = seq(0, 1, length.out = 22))) %>% {
      # print(.);.
      append(., list(quantile_id = seq_len(length(.[[1]]))))
    }
    )) %>% 
  unnest_wider(mean_output) %>%
  
  select(-output) %>% 
  mutate(parity_pcr_id = row_number()) %>% 
  
  unnest(c(quantile_id, a,b,k,d)) %>%
  select(parity_pcr_id, quantile_id, everything()) %>% 
  # select(parity_pcr_id, everything()) %>% 
  
  expand_grid(DIM = seq(5, 305)) %>% 
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% 
  
  nest(data = -c(parity_pcr_id, quantile_id)) %>% 
  
  mutate(output = data %>% map(. %>% filter(which.min(logSCC) <= row_number()))) %>% 
  select(-data) %>%
  # print()
  
  # # VALIDATE: We got the curves for each paring
  # unnest(output) %>%
  # identity() %>% {
  #   ggplot(.) +
  #     # aes(DIM, logSCC, group = interaction(PCR, PARITY)) +
  #     # aes(DIM, logSCC, group = parity_pcr_id, color = factor(parity_pcr_id)) +
  #     aes(DIM, logSCC, group = interaction(parity_pcr_id, quantile_id), 
  #         color = interaction(parity_pcr_id, quantile_id)) +
  # 
  #     # labs(color = "PCR x Parity") +
  #     labs(color = "PCR x Parity x Quantile id") +
  #     geom_line(size = rel(1.5)) + # blue for NEG
  #     guides(color = "none") +
  #     ggpubr::theme_classic2() +
  #     NULL
  # }
  
  # post means post minima
  mutate(model_post = output %>% map(~lm(logSCC ~ DIM, data = .x))) %>% 
  select(-output) %>% 
  mutate(model_tidy = model_post %>% map(. %>% broom::tidy())) %>% 
  # mutate(model_glance = model_post %>% map(broom::glance)) %>% 
  # unnest(model_glance)
  unnest(model_tidy) %>% 
  
  select(-model_post) %>% 
  identity() -> all_quan_sum_curves


# following not useful
all_quan_sum_curves %>% 
  pivot_wider(names_from = term, values_from = estimate,id_cols = c(parity_pcr_id, quantile_id)) %>% 
  identity() %>% {
    ggplot(.) +
      aes(quantile_id, DIM) + 
      
      geom_col() + 
      facet_wrap(~parity_pcr_id)
  }




# nlme output dataframe --------------------------------------

# dataframe withh all output - repeat for all groups :
## instead insert in coef pipeline: out_nlme_neg2
all_out_nlme_pos2 <- out_nlme_pos2 %>%
  add_column(PARITY = "2", PCR = "POS")


all_out_nlme <- rbind(all_out_nlme_neg2, all_out_nlme_neg3, all_out_nlme_neg4, 
                      all_out_nlme_pos2, all_out_nlme_pos3, all_out_nlme_pos4) 
rm(all_out_nlme_neg2, all_out_nlme_neg3, all_out_nlme_neg4, 
   all_out_nlme_pos2, all_out_nlme_pos3, all_out_nlme_pos4); gc()

all_out_nlme$PARITY <- factor(all_out_nlme$PARITY)
all_out_nlme$PCR <- factor(all_out_nlme$PCR)

all_out_nlme <- all_out_nlme %>% rownames_to_column("BES_ID")



# Boxplot and histograms of parameters ---------------------------------

# visualize nlme parameters boxplot:
p_a <- ggplot(all_out_nlme, aes(x=PARITY, y=a, fill=PCR)) + 
  geom_boxplot()
  # + theme(legend.position = "none")
p_b <- ggplot(all_out_nlme, aes(x=PARITY, y=b, fill=PCR)) + 
  geom_boxplot()
p_k <- ggplot(all_out_nlme, aes(x=PARITY, y=k, fill=PCR)) + 
  geom_boxplot() 
p_d <- ggplot(all_out_nlme, aes(x=PARITY, y=d, fill=PCR)) + 
  geom_boxplot() 

# boxplot for all output
all_p <- ggarrange(p_a, p_b, p_k, p_d, ncol=2, nrow=2, common.legend = TRUE, legend="right")


# include t.test (p values) in plot (all are significant)
pa_ttest <- p_a + 
  stat_compare_means(method = "t.test") +
  labs(title="parameter a, PCR POS vs NEG",x="Parity", y = "a values") 
pb_ttest <- p_b + 
  stat_compare_means(method = "t.test") +
  labs(title="parameter b, PCR POS vs NEG",x="Parity", y = "b values") 
pk_ttest <- p_k + 
  stat_compare_means(method = "t.test") +
  labs(title="parameter k, PCR POS vs NEG",x="Parity", y = "k values") 
pd_ttest <- p_d + 
  stat_compare_means(method = "t.test") +
  labs(title="parameter d, PCR POS vs NEG",x="Parity", y = "d values")


# histogram of distribution example overlaying:
hist_out_p2a <- all_out_nlme %>%
  filter(PARITY == 2) %>%
  ggplot(aes(x=a, fill= PCR, color=PCR)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Wilmink parameter a distribution, Parity 2, PCR POS vs NEG")




# Overlaying histograms, all output: create 12, one for each parameter for each parity:
h_a4 <- all_out_nlme %>%
  filter(PARITY == 4) %>%
  ggplot(aes(x=a, fill=PCR, color=PCR)) +
  geom_histogram(position="identity", alpha=0.5)+
  labs(x="Parity 4, parameter a")

hist_all <- ggarrange(h_a2, h_b2, h_k2, h_d2, 
                      h_a3, h_b3, h_k3, h_d3, 
                      h_a4, h_b4, h_k4, h_d4, 
                      ncol=4, nrow=3, common.legend = TRUE, legend="right")




# grid arrange histogram all parities ----------------------------

# Create for POS and NEG!

#ha_pos2 <- ggplot(out_nlme_pos2, aes(x=a))+
ha_neg2 <- ggplot(out_nlme_neg2, aes(x=a)) +
  #geom_histogram(color="darkred", fill="coral") + # pos colous
  geom_histogram(color="darkblue", fill="lightblue") + # neg colous
  #labs(x="a: Parity 2 POS") +
  labs(x="a: Parity 2 NEG") +
  xlim(2.0, 6.5)

#hb_pos2 <- ggplot(out_nlme_pos2, aes(x=b))+
hb_neg2 <- ggplot(out_nlme_neg2, aes(x=b))+
  #geom_histogram(color="darkred", fill="coral") + # pos colous
  geom_histogram(color="darkblue", fill="lightblue") +
  #labs(x="b: Parity 2 POS") +
  labs(x="b: Parity 2 NEG") +
  xlim(0.00065, 0.007)

#hk_pos2 <- ggplot(out_nlme_pos2, aes(x=k))+
hk_neg2 <- ggplot(out_nlme_neg2, aes(x=k))+
  #geom_histogram(color="darkred", fill="coral") + # pos colous
  geom_histogram(color="darkblue", fill="lightblue") +
  #labs(x="k: Parity 2 POS") +
  labs(x="k: Parity 2 NEG") +
  xlim(-3.3, -0.7)

#hd_pos2 <- ggplot(out_nlme_pos2, aes(x=d))+
hd_neg2 <- ggplot(out_nlme_neg2, aes(x=d))+
  #geom_histogram(color="darkred", fill="coral") + # pos colous
  geom_histogram(color="darkblue", fill="lightblue") +
  #labs(x="d: Parity 2 POS") +
  labs(x="d: Parity 2 NEG") +
  xlim(0.8, 4.0)


# gridarrange all histograms POS (red) and NEG (blue) Parity 2
### NOTE: better x axis in above histograms needed
h_all_p2 <- ggarrange(ha_pos2, hb_pos2, hk_pos2, hd_pos2, 
                   ha_neg2, hb_neg2, hk_neg2, hd_neg2, 
                   ncol=4, nrow=2)



# example of intercorrelation between parameters, NEG parity 2
gg_neg2 <- out_nlme_neg2 %>% 
  ggpairs()





# percentile curves based on wilmink ----------------------------------------

# repeat for all groups:
pos3_pct <- out_nlme_pos3 %>% 
  as_tibble() %>% 
  summarise(across(everything(), 
                   list(
                     "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1))))) %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything())


# repeat for all groups:
curve_pct_pos4 <- pos4_pct %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      aes(group = pctile) +
      aes(color = pctile) +
      geom_line() +
      
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("SCC curve, pctile of parametres Wilmink: PCR POS parity 4") +
      #ylim(3.0, 11.0) +
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()




# plot mean curves based on wilmink ---------------------------------------

# repeat for all groups:
neg2_mean <- out_nlme_neg2 %>% 
  summarise(across(everything(), mean))

# repeat for all groups:
curve_mean_pos2 <- pos2_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      #aes(group = pctile) +
      #aes(color = median) +
      geom_line(size = rel(1.5), colour = "#E69F00") + # orange for POS
      #geom_line(size = rel(1.5), colour = "#56B4E9") + # blue for NEG
      #labs(caption = "PCR negative, parity 2") +
      ggtitle("SCC curve POS parity 2: Mean of Wilmink parametres") +
      #ylim(3.5, 10.5) +
      
      ggpubr::theme_classic2() +
      NULL
  } %>% 
  plotly::ggplotly()
  




# retrieve values from fitted curve -------------------------------------

all_out_nlme <- all_out_nlme %>% as_tibble() 

all_out_nlme %>% 
  select(-BES_ID) %>% 
  expand_grid(DIM = 5:305) %>% 
  # arrange(DIM, PARITY, PCR, .by_group = FALSE) %>% 
  arrange(PARITY, PCR, DIM, .by_group = FALSE) %>% 
  # group_by(PARITY, PCR) %>% 
  # nest()
  identity()


# repeat for all groups:
est_pos2_mean <- 
  pos2_mean %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% 
  
  # identity{} : min SCC, DIM at 100, DIM at 150
  identity() %>% {
    bind_rows(
      slice_min(., logSCC, n = 1),
      filter(., DIM == 100),
      filter(., DIM == 150)
    )
  }
  # # group_by(BES_ID) %>% 
  # arrange(DIM) %>% 
  # filter(which.min(logSCC) <= row_number()) %>% 
  # # identity() %>% {
  # #   ggplot(.) + 
  # #     aes(DIM, logSCC) + 
  # #     geom_line(size = rel(1.5), colour = "#56B4E9") + # blue for NEG
  # #     ggtitle("SCC curve NEG parity 2: Mean of Wilmink parametres") +
  # #     
  # #     ggpubr::theme_classic2() +
  # #     NULL
  # # }
  # nest(data = everything()) %>% 
  # 
  # mutate(model_post = data %>% map(~lm(logSCC ~ DIM, data = .x))) %>% 
  # select(-data) %>% 
  # mutate(model_tidy = model_post %>% map(broom::tidy)) %>% 
  # # mutate(model_glance = model_post %>% map(broom::glance)) %>% 
  # # unnest(model_glance)
  # unnest(model_tidy)



# repeat for all groups
## instead insert in est pipeline est_neg3_mean <-
est2neg <- est_neg2_mean %>%
  add_column(PARITY = "2", PCR = "NEG")
est2pos <- est_pos2_mean %>%
  add_column(PARITY = "2", PCR = "POS")
est3neg <- est_neg3_mean %>%
  add_column(PARITY = "3", PCR = "NEG")
est3pos <- est_pos3_mean %>%
  add_column(PARITY = "3", PCR = "POS")
est4neg <- est_neg4_mean %>%
  add_column(PARITY = "4", PCR = "NEG")
est4pos <- est_pos4_mean %>%
  add_column(PARITY = "4", PCR = "POS")

all_est <- rbind(est2neg, est2pos, est3neg, est3pos, est4neg, est4pos) 

all_est$PARITY <- factor(all_est$PARITY)
all_est$PCR <- factor(all_est$PCR)


#--------------------------------------------------------------------------
# display plots: 

# display parameter distribution:
all_p     # grid arrange, boxplot all parametres
gg_neg2   # example intercorrelation df2_neg
pa_ttest  # example of boxplot "a parameteres with p values"
hist_out_p2a # example parameter a distribution: Overlaying only one parity: Parity 2 POS and NEG
h_all_p2     # histogram parity 2
hist_all  # histogram all output all parities: Overlapping POS and NEG

# curves examples
curve_mean_neg2
curve_mean_pos2
curve_pct_neg2
curve_pct_pos4

#---------------------------------------------------------------------------

# THE END
