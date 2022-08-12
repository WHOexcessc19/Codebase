
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


### load in and format Cordoba subnational deaths
argentina_cordoba_acm <- read_excel("Imported_Data/aregentina_cordoba_acm.xlsx")
cordoba_monthly <- argentina_cordoba_acm %>% 
  rename(time = month_num) %>% 
  pivot_longer(cols = starts_with("deaths_"), names_to = "year", values_to = "deaths_cor") %>% 
  mutate(year = str_sub(year,-4) %>% as.numeric)

### load in and format Argentina national deaths
argentina_monthly <- read_csv("https://raw.githubusercontent.com/akarlinsky/world_mortality/main/world_mortality.csv") %>% 
  filter(country_name == "Argentina") %>% 
  rename(deaths_arg = deaths)

### filter historical Argentina national data
argentina_monthly.mod <- argentina_monthly %>%
  filter(year>2018) %>%
  filter(year<2021)

### filter historical Cordoba subnational data
cordoba_monthly.mod <- cordoba_monthly %>%
  rename(month = time) %>%
  mutate(date = as.Date(as.yearmon(paste(year, month), "%Y %m"))) %>%
  arrange(date) %>%
  filter(year < 2021)

### filter Cordoba subnational predictions data
cordoba.pred <- cordoba_monthly %>%
  rename(month = time) %>%
  mutate(date = as.Date(as.yearmon(paste(year, month), "%Y %m"))) %>%
  arrange(date) %>%
  filter(year == 2021) %>%
  na.omit()


### Binomial Model for Subnational Estimates

#### Compile stan model ####
binomial_model  <- cmdstan_model('Excess/Stan_Code/Subnational_Binomial.stan')

### PC Prior Specification
pc.u <- 1
pc.alpha <- 0.01
### STAN Data
stan_data <- list(Tt=nrow(argentina_monthly.mod), # number of historical time points
                  N=nrow(cordoba.pred), # number of prediction time points
                  Y1=cordoba_monthly.mod$deaths_cor, # subnational historical deaths
                  Y=argentina_monthly.mod$deaths_arg, # national historical deaths
                  Y1T=cordoba.pred$deaths_cor, # subnational predictions deaths
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

### DF of Argentina samples where rows are 12 months of 2021 and columns are 1000 samples per month
Argentina.est = draws_Y2T
save(Argentina.est,file="Generated_Data/Argentina.est.RData")

### Preliminary Summary DF (doesn't take into account expecteds uncertainty)
Arg.who = mf.df %>% filter(Country=="Argentina") %>% filter(year==2021)
Argentina.preds <- data.frame(months=1:12,acm=apply(draws_Y2T,2,mean),lower=apply(draws_Y2T,2,quantile,probs=0.025),
                            upper=apply(draws_Y2T,2,quantile,probs=0.975),
                            excess=apply(draws_Y2T,2,mean)-Arg.who$expected,
                            excess.lower=apply(draws_Y2T,2,quantile,probs=0.025)-Arg.who$expected,
                            excess.upper=apply(draws_Y2T,2,quantile,probs=0.975)-Arg.who$expected)


