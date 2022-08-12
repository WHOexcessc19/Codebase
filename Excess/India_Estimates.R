
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

# load in subnational historical deaths
india_states_monthly <- read_csv("Imported_Data/india_states_monthly.csv")
# load in national historical expected deaths
load("Imported_Data/india_hist_monthly.RData")

# take out Gujarat state since no pre-pandemic data to fit the model on
india_states_monthly <- india_states_monthly %>%
  filter(state != "Gujarat State")

# national expecteds pre-pandemic
Ytplus = india_hist_monthly$gamma_E 

# totalling subnational pre-pandemic
Ytplusobs = india_states_monthly %>%
  group_by(time,year) %>%
  summarize(Ytplus=sum(deaths_scaled)) %>%
  arrange(year,time) %>%
  filter(year < 2020)
Ytplusobs = Ytplusobs$Ytplus

# making pre-pandemic grid of observed subnational deaths by state
Ytobs = india_states_monthly %>%
  filter(year < 2020)
Ytobs.grid = expand.grid(unique(india_states_monthly$state),1:12,2015:2019) %>% arrange(Var1) %>%
  rename(state = Var1,
         time = Var2,
         year = Var3)
Ytobs.grid = Ytobs.grid %>%
  left_join(Ytobs,by=c("time","year","state")) %>%
  mutate(Tt=rep(1:60,17)) %>%
  dplyr::select(-c(time,year,deaths,completeness_deaths)) %>%
  mutate(deaths_scaled = ifelse(is.na(deaths_scaled),0,deaths_scaled))
Yt = reshape(Ytobs.grid,direction = "wide",timevar = "state",idvar = "Tt") %>% dplyr::select(-Tt)
Yt = as.matrix(Yt)
Ytobs = ifelse(Yt==0,0,1) # making the grid of whether each subnational time point is observed

# totalling subnational for pandemic times
Ytplusobs_covid = india_states_monthly %>%
  mutate(time = ifelse(year==2020,time,time+12)) %>%
  group_by(time,year) %>%
  arrange(year,time) %>%
  summarize(Ytplus=sum(deaths_scaled)) %>%
  filter(year > 2019)
Ytplusobs_covid = Ytplusobs_covid$Ytplus

# making during pandemic grid of observed subnational deaths by state
Ytobs.2020 = india_states_monthly %>%
  filter(year > 2019) %>%
  mutate(time = ifelse(year==2020,time,time+12))
Ytobs.grid = expand.grid(unique(india_states_monthly$state),1:24) %>% arrange(Var1) %>%
  rename(state = Var1,
         time = Var2)
Ytobs.grid = Ytobs.grid %>%
  left_join(Ytobs.2020,by=c("time","state")) %>%
  mutate(Tt=rep(1:24,17)) %>%
  dplyr::select(-c(time,year,deaths,completeness_deaths)) %>%
  mutate(deaths_scaled = ifelse(is.na(deaths_scaled),0,deaths_scaled))
Yt.2020 = reshape(Ytobs.grid,direction = "wide",timevar = "state",idvar = "Tt") %>% dplyr::select(-Tt)
Yt.2020 = as.matrix(Yt.2020)
Ytobs.2020 = ifelse(Yt.2020==0,0,1) # making the grid of whether each subnational time point is observed

### Prior Specification
pc.u <- 1
pc.alpha <- 0.01

## Model Data list
stan_data <- list(Tt=60, # Number of non-COVID time points
                  K=17, # Number of Subnational Regions
                  Ytplus=Ytplus, # National Deaths in non-COVID times
                  Ytplusobs=Ytplusobs, # Sum of Subnational Deaths in non-COVID times
                  Yt=Yt, # Grid of Subnational Deaths in non-COVID times
                  Ytobs=Ytobs, # Grid of Subnational Missingness in non-COVID times
                  Tt_covid=24, # Number of COVID time points
                  Ytplusobs_covid=Ytplusobs_covid, # Sum of Subnational Deaths in COVID times
                  Yt_covid=Yt.2020, # Grid of Subnational Deaths in COVID times
                  Ytobs_covid=Ytobs.2020, # Grid of Subnational Missingness in COVID times
                  pc_U = pc.u, # Prior Parameters
                  pc_alpha = pc.alpha) # Prior Parameters

# Fit model
multin_model  <- cmdstan_model('Excess/Stan_Code/Multinomial_Subnational.stan')

stan_fit_subnational <- multin_model$sample(data = stan_data, seed = 60,
                                            chains = 4, parallel_chains = 2, 
                                            max_treedepth = 15, adapt_delta = 0.85,
                                            iter_warmup = 500, iter_sampling =250, save_warmup = FALSE)

# Draw from Posterior
draws_Ytplus_covid <- stan_fit_subnational$draws() %>% subset_draws(variable = "Ytplus_covid", regex = TRUE) %>%
  merge_chains() %>% as.data.frame()

### DF of India samples where rows are 24 months of 2020-2021 and columns are 1000 samples per month
India.est = draws_Y2T
save(India.est,file="Generated_Data/India.est.RData")

### Preliminary Summary DF (doesn't take into account expecteds uncertainty)
India.who = mf.df %>% filter(Country=="India")
India.preds <- data.frame(months=1:24,acm=apply(draws_Y2T,2,mean),lower=apply(draws_Y2T,2,quantile,probs=0.025),
                           upper=apply(draws_Y2T,2,quantile,probs=0.975),
                           excess=apply(draws_Y2T,2,mean)-India.who$expected,
                           excess.lower=apply(draws_Y2T,2,quantile,probs=0.025)-India.who$expected,
                           excess.upper=apply(draws_Y2T,2,quantile,probs=0.975)-India.who$expected)


