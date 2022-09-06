library(mgcv)
library(tidyverse)

# The purpose of this file is to fit a negative binomial spline model
# to the countries with monthly mortality data to get expected monthly
# mortality in 2020-2021

#### Set the seed for this script ####
set.seed(42)

#### Load in data ####
load("../Generated_Data/multinom.expected.dfs.Rda")
load("../Generated_Data/labels.Rda")

#### Clean data for use in gam function ####
exp_obs_bymonthyear <- exp_obs_bymonthyear %>% 
  select(iso3, year, month, observed) %>%
  filter(year < 2020) %>% 
  mutate(country_num = as.numeric(as.factor(iso3)),
         observed = floor(observed)) %>% 
  arrange(iso3, year, month)

#### Create data frame to store expected monthly mortality in 2020-2021 ####
# gamma_E and gamma_delta will contain, for each country time period, the 
# parameters of the gamma distribution for the expected mortality
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

#### Predict expected monthly mortality in 2020-2021 ####
for(i in 1:max(exp_obs_bymonthyear$country_num)){
  whichs <- which(exp_obs_bymonthyear$country_num == i)
  temp <- exp_obs_bymonthyear[whichs, ]
  
  # Fit gam
  # If there are less than 3 years of historical data, use a linear annual 
  # trend, rather than a spline trend
  if(length(unique(temp$year)) < 3){
    print(temp$iso3[1])
    annual_model <- gam(observed ~ year +
                          s(month, bs = "cc", k = length(unique(temp$month))),
                        data = temp, family = nb(theta = NULL, link = "log"))
  } else{
    annual_model <- gam(observed ~ s(year, k = length(unique(temp$year))) +
                          s(month, bs = "cc", k = length(unique(temp$month))),
                        data = temp, family = nb(theta = NULL, link = "log"))
  }
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
  
  # Get gamma parameters
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
  
  # Get fitted values for historical time periods
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

#### Save all expecteds, i.e. for pre 2020 and 2020-2021 ####
acm_predictions_forbind <- acm_predictions %>% 
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta)
acm_monthly_predictions_tier1_hist <- 
  exp_obs_bymonthyear %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, country_num, gamma_E, gamma_delta) %>%
  rbind(acm_predictions_forbind) %>% arrange(iso3, year, month)
save(acm_monthly_predictions_tier1_hist, 
     file = "../Generated_Data/acm_monthly_predictions_tier1_hist.RData")

#### Save expecteds for 2020-2021 ####
acm_monthly_predictions_tier1 <- acm_predictions %>%
  select(iso3, year, month, expected_acm, expected_acm_se, expected_log_acm,
         expected_log_acm_se, gamma_E, gamma_delta, gamma_E_nb, 
         gamma_delta_nb) %>%
  mutate(gamma_sd = sqrt((gamma_E ^ 2) / gamma_delta),
         gamma_sd_nb = sqrt((gamma_E_nb ^ 2) / gamma_delta_nb))
save(acm_monthly_predictions_tier1, 
     file = "../Generated_Data/acm_monthly_predictions_tier1.RData")
