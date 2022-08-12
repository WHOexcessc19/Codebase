
### load in packages
library(cmdstanr)
library(posterior)
library(zoo)
library(tidyverse)
library(readxl)

### load in full country df
load("Generated_Data/mf.df.Rda")
### load in expecteds estimates
load("Generated_Data/acm_monthly_predictions_tier1.RData")
load("Generated_Data/acm_monthly_predictions_tier2.RData")
### bind together expecteds dfs and merge into covariate model df
new.expected <- rbind(acm_monthly_predictions_tier1,acm_monthly_predictions_tier2)
mf.df <- mf.df %>%
  left_join(new.expected,by=c("iso3","year","month")) %>%
  dplyr::select(-expected) %>%
  rename(expected = expected_acm)

### read in Indonesia national and subnational deaths data
Indoensia_deaths_national <- read_excel("Dropbox/WHO-COVID-TAG/SubNational Data/Indonesia/Indoensia_deaths_national.xlsx")
indonesia_subnational <- read_excel("Dropbox/WHO-COVID-TAG/SubNational Data/Indonesia/indonesia_subnational.xlsx")

### filter national data, use GBD source
indonesia_national <- Indoensia_deaths_national %>%
  filter(location=="Indoensia" & source=="GBD")

### format subnational Jakarta death data
indonesia_monthly <- indonesia_subnational %>% 
  rename(time = month_num) %>% 
  pivot_longer(cols = starts_with("deaths_"), names_to = "year", values_to = "deaths_cor") %>% 
  mutate(year = str_sub(year,-4) %>% as.numeric)

### sum subnational to annual historical, since we only have national annual historical data
indonesia.sum_toYear <- indonesia_monthly %>%
  group_by(year) %>%
  summarise(Jakarta_Deaths = sum(deaths_cor,na.rm=TRUE)) %>%
  left_join(indonesia_national,by="year") %>%
  rename(National_Deaths = deaths)


### Binomial Model for Subnational Estimates

#### Compile stan model ####
binomial_model  <- cmdstan_model('Excess/Stan_Code/Subnational_Binomial.stan')

### PC Prior Specification
pc.u <- 1
pc.alpha <- 0.01
### STAN Data
stan_data <- list(Tt=5, # number of historical time points
                  N=18, # number of prediction time points
                  Y1=as.integer(indonesia.sum_toYear$Jakarta_Deaths[1:5]), # subnational annual historical deaths
                  Y=as.integer(indonesia.sum_toYear$National_Deaths[1:5]), # national annual historical deaths
                  Y1T=c(indonesia_subnational$deaths_2020,indonesia_subnational$deaths_2021[1:6]), # subnational monthly predictions deaths
                  pc_U = pc.u,
                  pc_alpha = pc.alpha)

### Fit model
stan_fit_subnational <- binomial_model$sample(data = stan_data, seed = 60,
                                             chains = 4, parallel_chains = 2, 
                                             max_treedepth = 15, adapt_delta = 0.85,
                                             iter_warmup = 500, iter_sampling = 250)
### Draws from the posterior for National ACM
draws_Y2T <- stan_fit_subnational$draws() %>% subset_draws(variable = "Y2T", regex = TRUE) %>%
  merge_chains() %>% as.data.frame()

### DF of Indonesia samples where rows are 18 months of 2020 - Jun 2021 and columns are 1000 samples per month
Indonesia.est = draws_Y2T
save(Indonesia.est,file="Generated_Data/Indonesia.est.RData")

### Preliminary Summary DF (doesn't take into account expecteds uncertainty)
Ind.who = mf.df %>% filter(Country=="Indonesia") %>% filter(!(year==2021 & month>6))
Indonesia.preds <- data.frame(months=1:18,acm=apply(draws_Y2T,2,mean),lower=apply(draws_Y2T,2,quantile,probs=0.025),
                            upper=apply(draws_Y2T,2,quantile,probs=0.975),
                            excess=apply(draws_Y2T,2,mean)-Ind.who$expected,
                            excess.lower=apply(draws_Y2T,2,quantile,probs=0.025)-Ind.who$expected,
                            excess.upper=apply(draws_Y2T,2,quantile,probs=0.975)-Ind.who$expected)
