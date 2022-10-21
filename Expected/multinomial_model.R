library(INLA)
library(tidyverse)
library(matrixStats)

# The purpose of this file is to fit a multinomial temperature model
# to the countries with monthly mortality data to predict expected monthly
# mortality in 2020-2021 for countries with only annual mortality data

#### Set the seed for this script ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(42)

#### Load in data ####
load("../Generated_Data/multinom.expected.dfs.Rda")
load("../Generated_Data/labels.Rda")
load("../Generated_Data/acm_annual_predictions_tier2.RData")
load("../Generated_Data/expected_acm_samples_tier2.RData")
load("../Expected/temperature_data_clean.RData")

#### Create a data frame containing monthly deaths counts and temperature ####
# Don't use data from 2020 or 2021
inla.dat <- merge(exp_obs_bymonthyear, temp_data, 
                  by = c("iso3", "year", "month")) %>% 
  select(iso3, year, month, observed, temperature) %>%
  filter(year < 2020) %>% 
  mutate(country_year = paste(iso3, year, sep = "_"),
         country_year_num = as.numeric(as.factor(country_year)),
         observed = floor(observed),
         month_char = as.character(month))

#### Fit multinomial model using Poisson-multinomial trick ####
pois.model <- 
  inla(observed ~ -1 + f(country_year_num, model = "iid", 
                         hyper = list(prec = list(initial = -25, fixed = TRUE)))
       + temperature,
       family = "poisson",
       data = inla.dat,
       control.compute = list(config = TRUE))

# Save samples for regression parameter
num_samples <- 10000
pois.samples <- inla.posterior.sample(n = num_samples, pois.model)
beta_samples <- rep(0, num_samples)
sample_ind <- nrow(pois.samples[[1]]$latent)
for(i in 1:num_samples){
  beta_samples[i] <-  pois.samples[[i]]$latent[sample_ind]
}
save(beta_samples, file = "../Generated_Data/beta_samples.RData")

# Extract beta without uncertainty
# With more than 1 coefficient need to change this code
beta <- pois.model$summary.fixed$mean
beta_sd <- pois.model$summary.fixed$sd

#### Check fitted monthly values to observed monthly values for historical data ####
inla.dat$fitted <- NA
for(i in 1:max(inla.dat$country_year_num)){
  whichs <- which(inla.dat$country_year_num == i)
  temp <- inla.dat[whichs, ]
  N <- sum(temp$observed)
  g <- exp(temp$temperature * beta)
  probs <- g / sum(g)
  inla.dat[whichs, "fitted"] <- N * probs
}
inla.dat$obsminfit <- inla.dat$observed - inla.dat$fitted
plot(inla.dat$observed, inla.dat$fitted)
plot(inla.dat$country_year_num, inla.dat$obsminfit)
plot(inla.dat$country_year_num, inla.dat$obsminfit / inla.dat$observed)


#### Predict monthly acm for countries with only annual acm available ####
# Use 2020 temperatures for 2021
annual_data_2020 <- merge(acm_annual_predictions_tier2, temp_data,
                     by = c("iso3", "year")) %>% 
  mutate(expected_monthly_acm = NA) 
temp_data_2021 <- temp_data %>% filter(year == 2020) %>% 
  mutate(year = 2021)
annual_data_2021 <- merge(acm_annual_predictions_tier2, temp_data_2021,
                          by = c("iso3", "year")) %>% 
  mutate(expected_monthly_acm = NA) 
annual_data <- rbind(annual_data_2020, annual_data_2021) %>% 
  arrange(iso3, year, month) %>%
  mutate(country_year = paste(iso3, year, sep = "_"),
         country_year_num = as.numeric(as.factor(country_year)),
         expected_monthly_acm_samples = NA,
         expected_monthly_acm_se_samples = NA,
         expected_monthly_log_acm_samples = NA,
         expected_monthly_log_acm_se_samples = NA,
         gamma_E = NA,
         gamma_delta = NA,
         gamma_E_nb = NA,
         gamma_delta_nb = NA)

# Get point estimates for expecteds
for(i in 1:max(annual_data$country_year_num)){
  whichs <- which(annual_data$country_year_num == i)
  temp <- annual_data[whichs, ]
  N <- temp$expected_acm[1]
  g <- exp(temp$temperature * beta)
  probs <- g / sum(g)
  annual_data[whichs, "expected_monthly_acm"] <- N * probs
}

# Get uncertainty for expecteds
acm_annual_predictions_tier2 <-  acm_annual_predictions_tier2 %>% 
  mutate(country_year = paste(iso3, year, sep = "_"),
         country_year_num = as.numeric(as.factor(country_year)))

annual_data_samples <- matrix(0, nrow = nrow(annual_data), ncol = num_samples)

expected_monthly_acm_samples <- 
  matrix(0, nrow = nrow(annual_data), ncol = num_samples)
expected_monthly_acm_samples_nb <- 
  matrix(0, nrow = nrow(annual_data), ncol = num_samples)

for(i in 1:max(annual_data$country_year_num)){
  whichs <- which(annual_data$country_year_num == i)
  temp <- annual_data[whichs, ]
  
  # Old
  annual_acm_samples <- rnorm(num_samples, mean = temp$expected_log_acm[1],
                              sd = temp$expected_log_acm_se[1])
  for(j in 1:num_samples){
    N <-  exp(annual_acm_samples[j])
    g <-  exp(temp$temperature * beta_samples[j])
    probs <- g / sum(g)
    annual_data_samples[whichs, j] <- N * probs
  }
  annual_data[whichs, "expected_monthly_acm_samples"] <- 
    rowMeans(annual_data_samples[whichs, ])
  annual_data[whichs, "expected_monthly_acm_se_samples"] <- 
    rowSds(annual_data_samples[whichs, ])
  annual_data[whichs, "expected_monthly_log_acm_samples"] <- 
    rowMeans(log(annual_data_samples[whichs, ]))
  annual_data[whichs, "expected_monthly_log_acm_se_samples"] <- 
    rowSds(log(annual_data_samples[whichs, ]))
  
  # New
  whichs_samples <- 
    which(acm_annual_predictions_tier2$country_year_num == 
            temp$country_year_num[1])
  annual_acm_samples <- expected_acm_samples[whichs_samples, ]
  annual_acm_samples_nb <- expected_acm_samples_nb[whichs_samples, ]
  for(j in 1:num_samples){
    N <-  annual_acm_samples[j]
    N_nb <-  annual_acm_samples_nb[j]
    g <-  exp(temp$temperature * beta_samples[j])
    probs <- g / sum(g)
    expected_monthly_acm_samples[whichs, j] <- N * probs
    expected_monthly_acm_samples_nb[whichs, j] <- N_nb * probs
  }
  
  gamma_E <- rep(0, 12)
  gamma_delta <- rep(0, 12)
  gamma_E_nb <- rep(0, 12)
  gamma_delta_nb <- rep(0, 12)
  samples <- expected_monthly_acm_samples[whichs, ]
  samples_nb <- expected_monthly_acm_samples_nb[whichs, ]
  for(j in 1:12){
    
    gamma_E[j] <- mean(samples[j, ])
    gamma_delta[j] <- ((gamma_E[j]) ^ 2) / var(samples[j, ])
    
    gamma_E_nb[j] <- mean(samples_nb[j, ])
    gamma_delta_nb[j] <- ((gamma_E_nb[j]) ^ 2) / var(samples_nb[j, ])
  }
  annual_data[whichs, "gamma_E"] <- gamma_E
  annual_data[whichs, "gamma_delta"] <- gamma_delta
  annual_data[whichs, "gamma_E_nb"] <- gamma_E_nb
  annual_data[whichs, "gamma_delta_nb"] <- gamma_delta_nb
}

#### Save expecteds for 2020-2021 ####
acm_monthly_predictions_tier2 <- annual_data %>% arrange(iso3, year, month) %>%
  mutate(expected_annual_acm = expected_acm,
         expected_annual_acm_se = expected_acm_se,
         expected_annual_log_acm = expected_log_acm,
         expected_annual_log_acm_se = expected_log_acm_se,
         expected_acm = expected_monthly_acm,
         expected_acm_se = expected_monthly_acm_se_samples,
         expected_log_acm = log(expected_monthly_acm),
         expected_log_acm_se = expected_monthly_log_acm_se_samples,
         gamma_sd = sqrt((gamma_E ^ 2) / gamma_delta),
         gamma_sd_nb = sqrt((gamma_E_nb ^ 2) / gamma_delta_nb)) %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, gamma_E, gamma_delta, gamma_E_nb, gamma_delta_nb,
         gamma_sd, gamma_sd_nb)
save(acm_monthly_predictions_tier2, 
     file = "./Generated_Data/acm_monthly_predictions_tier2.RData")

