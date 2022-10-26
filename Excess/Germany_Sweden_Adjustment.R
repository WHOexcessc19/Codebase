
# Germany/Sweden Expected-Excess Adjustment Code

### load in packages
library(mgcv)
library(tidyverse)
library(readxl)

#### Set the seed for this script ####
set.seed(42)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load in data
load("../Generated_Data/multinom.expected.dfs.Rda")
expected <- read_csv("../Generated_Data/expected.csv")
excess <- read_csv("../Generated_Data/excess.csv")
acm <- read_csv("../Generated_Data/acm.csv")



#### GERMANY

#### Get expecteds with linear term for year ####

# Clean data for use in gam function
exp_obs_bymonthyear <- exp_obs_bymonthyear %>%
  filter(iso3 == "DEU") %>%
  select(iso3, year, month, observed) %>%
  filter(year < 2020) %>% 
  mutate(country_num = as.numeric(as.factor(iso3)),
         observed = floor(observed)) %>% 
  arrange(iso3, year, month) 

# Created data frame to store expected monthly mortality in 2020-2021
countries <- unique(exp_obs_bymonthyear$iso3)
acm_predictions <- data.frame(iso3 = rep(countries, each = 24),
                              year = rep(c(rep(2020, 12), 
                                           rep(2021, 12)), 
                                         times = length(countries)),
                              month = rep(c(1:12, 1:12), 
                                          times = length(countries)),
                              expected_acm = NA,
                              expected_acm_se = NA,
                              expected_log_acm = NA,
                              expected_log_acm_se = NA,
                              gamma_E = NA,
                              gamma_delta = NA,
                              gamma_E_nb = NA,
                              gamma_delta_nb = NA) %>%
  mutate(country_num = as.numeric(as.factor(iso3)))

exp_obs_bymonthyear$expected_acm <- NA
exp_obs_bymonthyear$expected_acm_se <- NA
exp_obs_bymonthyear$expected_log_acm <- NA
exp_obs_bymonthyear$expected_log_acm_se <- NA
exp_obs_bymonthyear$gamma_E <- NA
exp_obs_bymonthyear$gamma_delta <- NA

num_samples <- 10000

# Predict expected monthly mortality in 2020-2021 with country specific models
for(i in c(1)){
  whichs <- which(exp_obs_bymonthyear$country_num == i)
  temp <- exp_obs_bymonthyear[whichs, ]
  
  # Fit gam
  annual_model <- gam(observed ~ year +
                        s(month, bs = "cc", k = length(unique(temp$month))),
                      data = temp, family = nb(theta = NULL, link = "log"))
  
  overd <- exp(annual_model$family$getTheta())
  
  # Get predictions 
  pred <- predict(annual_model,
                  se.fit = TRUE,
                  type = "response",
                  newdata = data.frame(year = c(rep(2020, 12),
                                                rep(2021, 12)),
                                       month = c(1:12, 1:12)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  pred_log <- predict(annual_model,
                      se.fit = TRUE,
                      newdata = data.frame(year = c(rep(2020, 12),
                                                    rep(2021, 12)),
                                           month = c(1:12, 1:12)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  gamma_E <- rep(0, 24)
  gamma_delta <- rep(0, 24)
  gamma_E_nb <- rep(0, 24)
  gamma_delta_nb <- rep(0, 24)
  for(j in 1:24){
    samples <- exp(rnorm(num_samples, mean = pred_log$fit[j], 
                         sd = pred_log$se.fit[j]))
    
    gamma_E[j] <- mean(samples)
    gamma_delta[j] <- ((gamma_E[j]) ^ 2) / var(samples)
    
    samples_nb <- rnbinom(num_samples, size = overd, mu = samples)
    
    gamma_E_nb[j] <- mean(samples_nb)
    gamma_delta_nb[j] <- ((gamma_E_nb[j]) ^ 2) / var(samples_nb)
  }
  acm_predictions[whichs_pred, "gamma_E"] <- gamma_E
  acm_predictions[whichs_pred, "gamma_delta"] <- gamma_delta
  acm_predictions[whichs_pred, "gamma_E_nb"] <- gamma_E_nb
  acm_predictions[whichs_pred, "gamma_delta_nb"] <- gamma_delta_nb
  
  # Get fitted values
  pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
  exp_obs_bymonthyear[whichs, "expected_acm"] <- pred_hist$fit
  exp_obs_bymonthyear[whichs, "expected_acm_se"] <- pred_hist$se.fit
  pred_log_hist <- predict(annual_model, se.fit = TRUE)
  exp_obs_bymonthyear[whichs, "expected_log_acm"] <- pred_log_hist$fit
  exp_obs_bymonthyear[whichs, "expected_log_acm_se"] <- pred_log_hist$se.fit
  
  num_hist <- length(pred_log_hist$fit)
  gamma_E_hist <- rep(0, num_hist)
  gamma_delta_hist <- rep(0, num_hist)
  for(j in 1:num_hist){
    samples <- exp(rnorm(num_samples, mean = pred_log_hist$fit[j], 
                         sd = pred_log_hist$se.fit[j]))
    
    gamma_E_hist[j] <- mean(samples)
    gamma_delta_hist[j] <- ((gamma_E_hist[j]) ^ 2) / var(samples)
    
  }
  exp_obs_bymonthyear[whichs, "gamma_E"] <- gamma_E_hist
  exp_obs_bymonthyear[whichs, "gamma_delta"] <- gamma_delta_hist
}

acm_predictions_forbind <- acm_predictions %>% 
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta)
germany_linear_historical <- 
  exp_obs_bymonthyear %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta) %>%
  rbind(acm_predictions_forbind) %>% arrange(iso3, year, month)


acm_predictions <- data.frame(iso3 = rep(countries, each = 24),
                              year = rep(c(rep(2020, 12), 
                                           rep(2021, 12)), 
                                         times = length(countries)),
                              month = rep(c(1:12, 1:12), 
                                          times = length(countries)),
                              expected_acm = NA,
                              expected_acm_se = NA,
                              expected_log_acm = NA,
                              expected_log_acm_se = NA,
                              gamma_E = NA,
                              gamma_delta = NA,
                              gamma_E_nb = NA,
                              gamma_delta_nb = NA) %>%
  mutate(country_num = as.numeric(as.factor(iso3)))

exp_obs_bymonthyear$expected_acm <- NA
exp_obs_bymonthyear$expected_acm_se <- NA
exp_obs_bymonthyear$expected_log_acm <- NA
exp_obs_bymonthyear$expected_log_acm_se <- NA
exp_obs_bymonthyear$gamma_E <- NA
exp_obs_bymonthyear$gamma_delta <- NA

num_samples <- 10000

# Predict expected monthly mortality in 2020-2021 with country specific models
for(i in c(1)){
  whichs <- which(exp_obs_bymonthyear$country_num == i)
  temp <- exp_obs_bymonthyear[whichs, ]
  
  # Fit gam
  annual_model <- gam(observed ~ s(year, k = length(unique(temp$year))) +
                        s(month, bs = "cc", k = length(unique(temp$month))),
                      data = temp, family = nb(theta = NULL, link = "log"))
  overd <- exp(annual_model$family$getTheta())
  
  # Get predictions 
  pred <- predict(annual_model,
                  se.fit = TRUE,
                  type = "response",
                  newdata = data.frame(year = c(rep(2020, 12),
                                                rep(2021, 12)),
                                       month = c(1:12, 1:12)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  pred_log <- predict(annual_model,
                      se.fit = TRUE,
                      newdata = data.frame(year = c(rep(2020, 12),
                                                    rep(2021, 12)),
                                           month = c(1:12, 1:12)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  gamma_E <- rep(0, 24)
  gamma_delta <- rep(0, 24)
  gamma_E_nb <- rep(0, 24)
  gamma_delta_nb <- rep(0, 24)
  for(j in 1:24){
    samples <- exp(rnorm(num_samples, mean = pred_log$fit[j], 
                         sd = pred_log$se.fit[j]))
    
    gamma_E[j] <- mean(samples)
    gamma_delta[j] <- ((gamma_E[j]) ^ 2) / var(samples)
    
    samples_nb <- rnbinom(num_samples, size = overd, mu = samples)
    
    gamma_E_nb[j] <- mean(samples_nb)
    gamma_delta_nb[j] <- ((gamma_E_nb[j]) ^ 2) / var(samples_nb)
  }
  acm_predictions[whichs_pred, "gamma_E"] <- gamma_E
  acm_predictions[whichs_pred, "gamma_delta"] <- gamma_delta
  acm_predictions[whichs_pred, "gamma_E_nb"] <- gamma_E_nb
  acm_predictions[whichs_pred, "gamma_delta_nb"] <- gamma_delta_nb
  
  # Get fitted values
  pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
  exp_obs_bymonthyear[whichs, "expected_acm"] <- pred_hist$fit
  exp_obs_bymonthyear[whichs, "expected_acm_se"] <- pred_hist$se.fit
  pred_log_hist <- predict(annual_model, se.fit = TRUE)
  exp_obs_bymonthyear[whichs, "expected_log_acm"] <- pred_log_hist$fit
  exp_obs_bymonthyear[whichs, "expected_log_acm_se"] <- pred_log_hist$se.fit
  
  num_hist <- length(pred_log_hist$fit)
  gamma_E_hist <- rep(0, num_hist)
  gamma_delta_hist <- rep(0, num_hist)
  for(j in 1:num_hist){
    samples <- exp(rnorm(num_samples, mean = pred_log_hist$fit[j], 
                         sd = pred_log_hist$se.fit[j]))
    
    gamma_E_hist[j] <- mean(samples)
    gamma_delta_hist[j] <- ((gamma_E_hist[j]) ^ 2) / var(samples)
    
  }
  exp_obs_bymonthyear[whichs, "gamma_E"] <- gamma_E_hist
  exp_obs_bymonthyear[whichs, "gamma_delta"] <- gamma_delta_hist
}

acm_predictions_forbind <- acm_predictions %>% 
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta)
germany_historical_spline <- 
  exp_obs_bymonthyear %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta) %>%
  rbind(acm_predictions_forbind) %>% arrange(iso3, year, month)



#### Get results ####

load("../Generated_Data/multinom.expected.dfs.Rda")

all_dat <- merge(germany_historical_spline, exp_obs_bymonthyear %>%
                   filter(iso3 == "DEU"),
                 by = c("iso3", "year", "month"), all = TRUE) %>%
  mutate(observed_acm = observed) %>%
  dplyr::select(iso3, year, month, expected_acm, country_num,
         observed_acm, gamma_E, gamma_delta, expected_log_acm, 
         expected_log_acm_se) %>%
  arrange(iso3, year, month)


all_dat <- all_dat %>%
  mutate(gamma_E_spline = gamma_E,
         gamma_delta_spline = gamma_delta) %>%
  select(-c(gamma_E, gamma_delta)) %>%
  cbind(germany_linear_historical %>% select(gamma_E, gamma_delta)) %>%
  mutate(gamma_E_linear = gamma_E, gamma_delta_linear = gamma_delta) %>%
  select(iso3, year, month, observed_acm, gamma_E_spline, gamma_delta_spline,
         gamma_E_linear, gamma_delta_linear) 

num_samps <- 10000
excess_samples_spline <- germany_linear_historical %>% 
  select(iso3, year, month) %>%
  cbind(matrix(NA, nrow = nrow(germany_linear_historical), ncol = num_samps))
excess_samples_linear <- excess_samples_spline

expected_samples_linear <- excess_samples_spline

set.seed(42)
for(i in 1:nrow(all_dat)){
  excess_samples_spline[i, 4:10003] <- all_dat$observed_acm[i] -
    rgamma(num_samps, shape = all_dat$gamma_delta_spline[i], 
           rate = all_dat$gamma_delta_spline[i] / 
             all_dat$gamma_E_spline[i])
  temp <- rgamma(num_samps, shape = all_dat$gamma_delta_linear[i], 
                 rate = all_dat$gamma_delta_linear[i] / 
                   all_dat$gamma_E_linear[i])
  excess_samples_linear[i, 4:10003] <- all_dat$observed_acm[i] - temp
  expected_samples_linear[i, 4:10003] <- temp
  
}
excess_samples_linear_germany <- excess_samples_linear %>% 
  filter(year > 2019)
expected_samples_linear_germany <- expected_samples_linear
save(expected_samples_linear_germany, 
     file = "../Generated_Data/germany_linear_expected_samples.RData")

expected[which(expected$Country=="Germany"),4:1003] = expected_samples_linear[61:84,4:1003]
excess[which(excess$Country=="Germany"),4:1003] = acm[which(acm$Country=="Germany"),4:1003] - expected_samples_linear[61:84,4:1003]




#### SWEDEN

#### Get expecteds with linear term for year ####

load("../Generated_Data/multinom.expected.dfs.Rda")

# Updated data
swe_data <- readxl::read_excel("../Imported_Data/swe data.xlsx", 
                               sheet = "Consultation data")

#### Get expecteds with updated data ####

# Clean data for use in gam function
exp_obs_bymonthyear <- exp_obs_bymonthyear %>%
  filter(iso3 == "SWE") %>%
  select(iso3, year, month, observed) %>%
  filter(year < 2020) %>% 
  mutate(country_num = as.numeric(as.factor(iso3)),
         observed = floor(observed)) %>% 
  arrange(iso3, year, month) 


# Created data frame to store expected monthly mortality in 2020-2021
countries <- unique(exp_obs_bymonthyear$iso3)
acm_predictions <- data.frame(iso3 = rep(countries, each = 24),
                              year = rep(c(rep(2020, 12), 
                                           rep(2021, 12)), 
                                         times = length(countries)),
                              month = rep(c(1:12, 1:12), 
                                          times = length(countries)),
                              expected_acm = NA,
                              expected_acm_se = NA,
                              expected_log_acm = NA,
                              expected_log_acm_se = NA,
                              gamma_E = NA,
                              gamma_delta = NA,
                              gamma_E_nb = NA,
                              gamma_delta_nb = NA) %>%
  mutate(country_num = as.numeric(as.factor(iso3)))

exp_obs_bymonthyear$expected_acm <- NA
exp_obs_bymonthyear$expected_acm_se <- NA
exp_obs_bymonthyear$expected_log_acm <- NA
exp_obs_bymonthyear$expected_log_acm_se <- NA
exp_obs_bymonthyear$gamma_E <- NA
exp_obs_bymonthyear$gamma_delta <- NA

num_samples <- 10000

# Predict expected monthly mortality in 2020-2021 with country specific models
for(i in c(1)){
  whichs <- which(exp_obs_bymonthyear$country_num == i)
  temp <- exp_obs_bymonthyear[whichs, ]
  
  # Fit gam
  # annual_model <- gam(observed ~ s(year, k = length(unique(temp$year))) +
  #                       s(month, bs = "cc", k = length(unique(temp$month))),
  #                     data = temp, family = nb(theta = NULL, link = "log"))
  annual_model <- gam(observed ~ year +
                        s(month, bs = "cc", k = length(unique(temp$month))),
                      data = temp, family = nb(theta = NULL, link = "log"))
  overd <- exp(annual_model$family$getTheta())
  
  # Get predictions 
  pred <- predict(annual_model,
                  se.fit = TRUE,
                  type = "response",
                  newdata = data.frame(year = c(rep(2020, 12),
                                                rep(2021, 12)),
                                       month = c(1:12, 1:12)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  pred_log <- predict(annual_model,
                      se.fit = TRUE,
                      newdata = data.frame(year = c(rep(2020, 12),
                                                    rep(2021, 12)),
                                           month = c(1:12, 1:12)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  gamma_E <- rep(0, 24)
  gamma_delta <- rep(0, 24)
  gamma_E_nb <- rep(0, 24)
  gamma_delta_nb <- rep(0, 24)
  for(j in 1:24){
    samples <- exp(rnorm(num_samples, mean = pred_log$fit[j], 
                         sd = pred_log$se.fit[j]))
    
    gamma_E[j] <- mean(samples)
    gamma_delta[j] <- ((gamma_E[j]) ^ 2) / var(samples)
    
    samples_nb <- rnbinom(num_samples, size = overd, mu = samples)
    
    gamma_E_nb[j] <- mean(samples_nb)
    gamma_delta_nb[j] <- ((gamma_E_nb[j]) ^ 2) / var(samples_nb)
  }
  acm_predictions[whichs_pred, "gamma_E"] <- gamma_E
  acm_predictions[whichs_pred, "gamma_delta"] <- gamma_delta
  acm_predictions[whichs_pred, "gamma_E_nb"] <- gamma_E_nb
  acm_predictions[whichs_pred, "gamma_delta_nb"] <- gamma_delta_nb
  
  # Get fitted values
  pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
  exp_obs_bymonthyear[whichs, "expected_acm"] <- pred_hist$fit
  exp_obs_bymonthyear[whichs, "expected_acm_se"] <- pred_hist$se.fit
  pred_log_hist <- predict(annual_model, se.fit = TRUE)
  exp_obs_bymonthyear[whichs, "expected_log_acm"] <- pred_log_hist$fit
  exp_obs_bymonthyear[whichs, "expected_log_acm_se"] <- pred_log_hist$se.fit
  
  num_hist <- length(pred_log_hist$fit)
  gamma_E_hist <- rep(0, num_hist)
  gamma_delta_hist <- rep(0, num_hist)
  for(j in 1:num_hist){
    samples <- exp(rnorm(num_samples, mean = pred_log_hist$fit[j], 
                         sd = pred_log_hist$se.fit[j]))
    
    gamma_E_hist[j] <- mean(samples)
    gamma_delta_hist[j] <- ((gamma_E_hist[j]) ^ 2) / var(samples)
    
  }
  exp_obs_bymonthyear[whichs, "gamma_E"] <- gamma_E_hist
  exp_obs_bymonthyear[whichs, "gamma_delta"] <- gamma_delta_hist
}

acm_predictions_forbind <- acm_predictions %>% 
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta)
sweden_historical_linear <- 
  exp_obs_bymonthyear %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta) %>%
  rbind(acm_predictions_forbind) %>% arrange(iso3, year, month)


load("../Generated_Data/multinom.expected.dfs.Rda")
# Clean data for use in gam function
exp_obs_bymonthyear <- exp_obs_bymonthyear %>%
  filter(iso3 == "SWE") %>%
  select(iso3, year, month, observed) %>%
  filter(year < 2020) %>% 
  mutate(country_num = as.numeric(as.factor(iso3)),
         observed = floor(observed)) %>% 
  arrange(iso3, year, month) 

# Created data frame to store expected monthly mortality in 2020-2021
countries <- unique(exp_obs_bymonthyear$iso3)
acm_predictions <- data.frame(iso3 = rep(countries, each = 24),
                              year = rep(c(rep(2020, 12), 
                                           rep(2021, 12)), 
                                         times = length(countries)),
                              month = rep(c(1:12, 1:12), 
                                          times = length(countries)),
                              expected_acm = NA,
                              expected_acm_se = NA,
                              expected_log_acm = NA,
                              expected_log_acm_se = NA,
                              gamma_E = NA,
                              gamma_delta = NA,
                              gamma_E_nb = NA,
                              gamma_delta_nb = NA) %>%
  mutate(country_num = as.numeric(as.factor(iso3)))

exp_obs_bymonthyear$expected_acm <- NA
exp_obs_bymonthyear$expected_acm_se <- NA
exp_obs_bymonthyear$expected_log_acm <- NA
exp_obs_bymonthyear$expected_log_acm_se <- NA
exp_obs_bymonthyear$gamma_E <- NA
exp_obs_bymonthyear$gamma_delta <- NA

num_samples <- 10000

# Predict expected monthly mortality in 2020-2021 with country specific models
for(i in c(1)){
  whichs <- which(exp_obs_bymonthyear$country_num == i)
  temp <- exp_obs_bymonthyear[whichs, ]
  
  # Fit gam
  annual_model <- gam(observed ~ s(year, k = length(unique(temp$year))) +
                        s(month, bs = "cc", k = length(unique(temp$month))),
                      data = temp, family = nb(theta = NULL, link = "log"))
  # annual_model <- gam(observed ~ year +
  #                       s(month, bs = "cc", k = length(unique(temp$month))),
  #                     data = temp, family = nb(theta = NULL, link = "log"))
  overd <- exp(annual_model$family$getTheta())
  
  # Get predictions 
  pred <- predict(annual_model,
                  se.fit = TRUE,
                  type = "response",
                  newdata = data.frame(year = c(rep(2020, 12),
                                                rep(2021, 12)),
                                       month = c(1:12, 1:12)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  pred_log <- predict(annual_model,
                      se.fit = TRUE,
                      newdata = data.frame(year = c(rep(2020, 12),
                                                    rep(2021, 12)),
                                           month = c(1:12, 1:12)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  gamma_E <- rep(0, 24)
  gamma_delta <- rep(0, 24)
  gamma_E_nb <- rep(0, 24)
  gamma_delta_nb <- rep(0, 24)
  for(j in 1:24){
    samples <- exp(rnorm(num_samples, mean = pred_log$fit[j], 
                         sd = pred_log$se.fit[j]))
    
    gamma_E[j] <- mean(samples)
    gamma_delta[j] <- ((gamma_E[j]) ^ 2) / var(samples)
    
    samples_nb <- rnbinom(num_samples, size = overd, mu = samples)
    
    gamma_E_nb[j] <- mean(samples_nb)
    gamma_delta_nb[j] <- ((gamma_E_nb[j]) ^ 2) / var(samples_nb)
  }
  acm_predictions[whichs_pred, "gamma_E"] <- gamma_E
  acm_predictions[whichs_pred, "gamma_delta"] <- gamma_delta
  acm_predictions[whichs_pred, "gamma_E_nb"] <- gamma_E_nb
  acm_predictions[whichs_pred, "gamma_delta_nb"] <- gamma_delta_nb
  
  # Get fitted values
  pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
  exp_obs_bymonthyear[whichs, "expected_acm"] <- pred_hist$fit
  exp_obs_bymonthyear[whichs, "expected_acm_se"] <- pred_hist$se.fit
  pred_log_hist <- predict(annual_model, se.fit = TRUE)
  exp_obs_bymonthyear[whichs, "expected_log_acm"] <- pred_log_hist$fit
  exp_obs_bymonthyear[whichs, "expected_log_acm_se"] <- pred_log_hist$se.fit
  
  num_hist <- length(pred_log_hist$fit)
  gamma_E_hist <- rep(0, num_hist)
  gamma_delta_hist <- rep(0, num_hist)
  for(j in 1:num_hist){
    samples <- exp(rnorm(num_samples, mean = pred_log_hist$fit[j], 
                         sd = pred_log_hist$se.fit[j]))
    
    gamma_E_hist[j] <- mean(samples)
    gamma_delta_hist[j] <- ((gamma_E_hist[j]) ^ 2) / var(samples)
    
  }
  exp_obs_bymonthyear[whichs, "gamma_E"] <- gamma_E_hist
  exp_obs_bymonthyear[whichs, "gamma_delta"] <- gamma_delta_hist
}

acm_predictions_forbind <- acm_predictions %>% 
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta)
sweden_historical_spline <- 
  exp_obs_bymonthyear %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta) %>%
  rbind(acm_predictions_forbind) %>% arrange(iso3, year, month)


#### Get results ####

load("../Generated_Data/multinom.expected.dfs.Rda")

all_dat <- merge(sweden_historical_spline, exp_obs_bymonthyear %>%
                   filter(iso3 == "SWE"),
                 by = c("iso3", "year", "month"), all = TRUE) %>%
  mutate(observed_acm = observed) %>%
  select(iso3, year, month, expected_acm, country_num,
         observed_acm, gamma_E, gamma_delta, expected_log_acm, 
         expected_log_acm_se) %>%
  arrange(iso3, year, month)



all_dat <- all_dat %>%
  mutate(gamma_E_spline = gamma_E,
         gamma_delta_spline = gamma_delta) %>%
  select(-c(gamma_E, gamma_delta)) %>%
  cbind(sweden_historical_linear %>% select(gamma_E, gamma_delta)) %>%
  mutate(gamma_E_linear = gamma_E, gamma_delta_linear = gamma_delta) %>%
  select(iso3, year, month, observed_acm, gamma_E_spline, gamma_delta_spline,
         gamma_E_linear, gamma_delta_linear) 

num_samps <- 10000
excess_samples_spline <- sweden_historical_linear %>% 
  select(iso3, year, month) %>%
  cbind(matrix(NA, nrow = nrow(sweden_historical_linear), ncol = num_samps))
excess_samples_linear <- excess_samples_spline

expected_samples_linear <- excess_samples_spline

set.seed(42)
for(i in 1:nrow(all_dat)){
  excess_samples_spline[i, 4:10003] <- all_dat$observed_acm[i] -
    rgamma(num_samps, shape = all_dat$gamma_delta_spline[i], 
           rate = all_dat$gamma_delta_spline[i] / 
             all_dat$gamma_E_spline[i])
  temp <- rgamma(num_samps, shape = all_dat$gamma_delta_linear[i], 
                 rate = all_dat$gamma_delta_linear[i] / 
                   all_dat$gamma_E_linear[i])
  excess_samples_linear[i, 4:10003] <- all_dat$observed_acm[i] - temp
  expected_samples_linear[i, 4:10003] <- temp
}
excess_samples_linear_sweden <- excess_samples_linear %>% 
  filter(year > 2019)
expected_samples_linear_sweden <- expected_samples_linear
save(expected_samples_linear_sweden, 
     file = "../Generated_Data/sweden_linear_expected_samples.RData")

expected[which(expected$Country=="Sweden"),4:1003] = expected_samples_linear[61:84,4:1003]
excess[which(excess$Country=="Sweden"),4:1003] = acm[which(acm$Country=="Sweden"),4:1003] - expected_samples_linear[61:84,4:1003]





### Update the excess, expected, excess_summarized, expected_summarized csv within Generated_Data folder

excess_summarized = cbind(excess[,1:3],est=rowMeans(excess[,4:(num_inla_samps+3)]),
                          lwr=apply(excess[,4:(num_inla_samps+3)],1,quantile,0.025),
                          uppr=apply(excess[,4:(num_inla_samps+3)],1,quantile,0.975))
expected_summarized = cbind(expected[,1:3],est=rowMeans(expected[,4:(num_inla_samps+3)]),
                            lwr=apply(expected[,4:(num_inla_samps+3)],1,quantile,0.025),
                            uppr=apply(expected[,4:(num_inla_samps+3)],1,quantile,0.975))

### Write Estimates Full Sample DFs and Summarized Estimate DFs to csv
write.csv(excess, file="../Generated_Data/excess.csv", row.names=FALSE)
write.csv(expected, file="../Generated_Data/expected.csv", row.names=FALSE)

write.csv(excess_summarized, file="../Generated_Data/excess_summarized.csv", row.names=FALSE)
write.csv(expected_summarized, file="../Generated_Data/expected_summarized.csv", row.names=FALSE)


