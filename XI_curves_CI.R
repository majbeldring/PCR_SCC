
# PCR project
# adding CI shades to curves


# settings and data ------------------------------
library("tidyverse")
library("nlme")
library("pbapply")

(load("K:/paperI/major4/Ouput_nlme.RData"))



# summarize --------------------------------------

# summarise mean for all 8 models
neg1_mean <- nlme_out_neg1 %>% 
  dplyr::summarise(across(everything(), mean))
pos1_mean <- nlme_out_pos1 %>% 
  dplyr::summarise(across(everything(), mean))



# wilmink function -----------------------------------------

f_wilmink <- function(DIM, a,b,k,d){
  a + b * DIM + exp(-(exp(k)) * DIM)*d
}



# Curves -----------------------------------------

# Trying to add a ribbon to the plot
# geom_ribbon(aes(ymin=low, ymax=high), alpha=0.3) +
# However, fails at retrieving lower and upper CI values from the NLME output
# and apply it to the plot, with both POS and NEG are presented at the same time:

# Parity 1, POS and NEG curve:
list(pos1_mean = pos1_mean, neg1_mean = neg1_mean) %>% 
  enframe() %>% 
  unnest(value) %>% 
  crossing(DIM = seq_len(305)) %>% 
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a, b, k, d), f_wilmink)) %>% {
    ggplot(., aes(DIM, logSCC, group = name, color=name)) + 
      geom_line(size = rel(1.5)) +
      ylim(3.5, 6.0) +
      ggpubr::theme_classic2() + 
      labs(title = "Parity 1", color= "PCR") + 
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5),
            #legend.position="none"
            #,legend.position = "top"
      ) +
      scale_color_discrete(breaks=c("neg1_mean", "pos1_mean"),
                           labels=c("NEG", "POS")) +
      geom_hline(yintercept=5.3, linetype="dashed") + #200.000 threshold
      NULL
  } 







# pct plot -  Alternative ? ----------------------------------

nlme_out_pos1 %>% 
  as_tibble() %>% 
  dplyr::summarise(across(everything(), 
                          list(
                            "percentile" = ~ quantile(.x, probs = seq(0.1, .9, 0.1))))) %>% 
  dplyr::rename_with(~ str_remove(.x, "_percentile"), ends_with("_percentile")) %>% 
  mutate(pctile = names(a)) %>% select(pctile, everything()) %>% 
  #' join parameters with x-axis (`DIM`)
  crossing(DIM = seq_len(305)) %>% 
  
  #' calculate the proper `logSCC`
  mutate(logSCC = pmap_dbl(select(., DIM, a,b,k,d), f_wilmink)) %>% {
    ggplot(., aes(DIM, logSCC)) + 
      # aes(group = BES_ID) +
      aes(group = pctile) +
      aes(color = pctile) +
      geom_line() +
      ggtitle("Parity 1, PCR POS. Percentiles of Wilmink Parameters") +
      #ylim(3.0, 11.0) +
      ggpubr::theme_classic2() +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5),
            #legend.position="none"
            #,legend.position = "top"
      ) +
      NULL
  }



