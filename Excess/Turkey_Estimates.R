
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

### load in and format subnational historical deaths
turkey_deaths_subnational <- read_excel("Imported_Data/turkey_deaths_subnational.xlsx")
turkey_subnational_monthly <- turkey_deaths_subnational %>% 
  rename(month = month_num) %>% 
  pivot_longer(cols = starts_with("deaths_"), names_to = "year", values_to = "deaths_cor") %>% 
  mutate(year = str_sub(year,-4) %>% as.numeric) %>%
  mutate(date = as.Date(as.yearmon(paste(year, month), "%Y %m"))) %>%
  arrange(date) %>%
  filter(year>2017) %>%
  filter(year<2020)

### filter to subnational prediction deaths
turkey.pred <- turkey_deaths_subnational %>%
  rename(month = month_num) %>%
  pivot_longer(cols = starts_with("deaths_"), names_to = "year", values_to = "deaths_cor") %>% 
  mutate(year = str_sub(year,-4) %>% as.numeric) %>%
  mutate(date = as.Date(as.yearmon(paste(year, month), "%Y %m"))) %>%
  arrange(date) %>%
  filter(year == 2020|year == 2021) %>%
  na.omit()

### load in and format Turkey national deaths
turkey_deaths_national <- read_excel("Imported_Data/turkey_deaths_national.xlsx")
turkey_national_monthly <- turkey_deaths_national %>% 
  rename(month = month_num) %>% 
  pivot_longer(cols = starts_with("deaths_"), names_to = "year", values_to = "deaths_cor") %>% 
  mutate(year = str_sub(year,-4) %>% as.numeric) %>%
  mutate(date = as.Date(as.yearmon(paste(year, month), "%Y %m"))) %>%
  arrange(date) %>%
  filter(year>2017) %>%
  filter(year<2021)


### Binomial Model for Subnational Estimates

#### Compile stan model ####
binomial_model  <- cmdstan_model('Excess/Stan_Code/Subnational_Binomial.stan')

### PC Prior Specification
pc.u <- 1
pc.alpha <- 0.01
### STAN Data
stan_data <- list(Tt=nrow(turkey_national_monthly), # number of historical time points
                  N=nrow(turkey.pred), # number of prediction time points
                  Y1=as.integer(turkey_subnational_monthly$deaths_cor), # subnational historical deaths
                  Y=as.integer(turkey_national_monthly$deaths_cor), # national historical deaths
                  Y1T=as.integer(turkey.pred$deaths_cor), # subnational predictions deaths
                  pc_U = pc.u,
                  pc_alpha = pc.alpha)

# Fit model
stan_fit_subnational <- binomial_model$sample(data = stan_data, seed = 60,
                                             chains = 4, parallel_chains = 2, 
                                             max_treedepth = 15, adapt_delta = 0.85,
                                             iter_warmup = 500, iter_sampling = 250)
### Draws from the posterior for National ACM
draws_Y2T <- stan_fit_subnational$draws() %>% subset_draws(variable = "Y2T", regex = TRUE) %>%
  merge_chains() %>% as.data.frame()

### DF of Turkey samples where rows are 24 months of 2020-2021 and columns are 1000 samples per month
Turkey.est = draws_Y2T
save(Turkey.est,file="Generated_Data/Turkey.est.RData")

### Preliminary Summary DF (doesn't take into account expecteds uncertainty)
Turk.who = mf.df %>% filter(Country=="Turkey")
Turkey.preds <- data.frame(months=1:24,acm=apply(draws_Y2T,2,mean),lower=apply(draws_Y2T,2,quantile,probs=0.025),
                            upper=apply(draws_Y2T,2,quantile,probs=0.975),
                            excess=apply(draws_Y2T,2,mean)-Turk.who$expected,
                            excess.lower=apply(draws_Y2T,2,quantile,probs=0.025)-Turk.who$expected,
                            excess.upper=apply(draws_Y2T,2,quantile,probs=0.975)-Turk.who$expected)

