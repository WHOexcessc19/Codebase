library(mgcv)
library(tidyverse)

# The purpose of this file is to fit a negative binomial spline model
# to the countries with only annual mortality data to get expected annual
# mortality in 2020-2021

#### Set the seed for this script ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(42)

#### Load in data ####
load("../Generated_Data/multinom.expected.dfs.Rda")
load("../Generated_Data/labels.Rda")
AnnualMortality <- read_csv("../Imported_Data/AnnualMortality.csv")

#### Clean data for use in gam function ####
countries_monthly <- unique(exp_obs_bymonthyear$iso3)
countries <- setdiff(unique(acm.byyear$iso3), countries_monthly)

acm.byyear <- acm.byyear %>% filter(iso3 %in% countries) %>%
  mutate(country_num = as.numeric(as.factor(iso3)),
         deaths = floor(deaths)) 

acm.byyear$expected_acm <- NA
acm.byyear$expected_acm_se <- NA
acm.byyear$expected_log_acm <- NA
acm.byyear$expected_log_acm_se <- NA

#### Create data frame to store expected annual mortality in 2020-2021 ####
acm_predictions <- data.frame(iso3 = rep(countries, each = 2),
                              year = rep(c(2020, 2021), 
                                         times = length(countries)),
                              expected_acm = NA,
                              expected_acm_se = NA,
                              expected_log_acm = NA,
                              expected_log_acm_se = NA) %>%
  mutate(country_num = as.numeric(as.factor(iso3)))

num_samples <- 10000
expected_acm_samples <- matrix(0, nrow = nrow(acm_predictions), 
                               ncol = num_samples)
expected_acm_samples_nb <- matrix(0, nrow = nrow(acm_predictions), 
                                  ncol = num_samples)

#### Predict expected annual mortality in 2020-2021 ####
for(i in 1:max(acm.byyear$country_num)){
  whichs <- which(acm.byyear$country_num == i)
  temp <- acm.byyear[whichs, ]
  
  # For China, use data from AnnualMortality file
  # For Indonesia, do not use data from 2004
  if(temp$iso3[1] == "CHN"){
    whichs <- which(acm.byyear$country_num == i & acm.byyear$year > 2014)
    temp <- temp %>% filter(year > 2014)
    temp_chn <- AnnualMortality %>% filter(country == "China" & year < 2020)
    temp$deaths <- temp_chn$deaths
    
    # Fit gam
    # China only has 5 years of historical data from WMD, need to change df
    annual_model <- gam(deaths ~ s(year, k = nrow(temp)), data = temp,
                        family = nb(theta = NULL, link = "log"))
  } else if(temp$iso3[1] == "IDN"){
    temp <- temp %>% filter(year != 2004)
    annual_model <- gam(deaths ~ s(year), data = temp,
                        family = nb(theta = NULL, link = "log"))
  }
  else{
    # Fit gam
    annual_model <- gam(deaths ~ s(year), data = temp,
                        family = nb(theta = NULL, link = "log"))
  }
  overd <- exp(annual_model$family$getTheta())
  
  # Get predictions 
  pred <- predict(annual_model, se.fit = TRUE, type = "response",
                  newdata = data.frame(year = c(2020, 2021)))
  whichs_pred <- which(acm_predictions$country_num == i)
  acm_predictions[whichs_pred, "expected_acm"] <- pred$fit
  acm_predictions[whichs_pred, "expected_acm_se"] <- pred$se.fit
  
  pred_log <- predict(annual_model, se.fit = TRUE, 
                      newdata = data.frame(year = c(2020, 2021)))
  acm_predictions[whichs_pred, "expected_log_acm"] <- pred_log$fit
  acm_predictions[whichs_pred, "expected_log_acm_se"] <- pred_log$se.fit
  
  for(j in 1:2){
    samples <- exp(rnorm(num_samples, mean = pred_log$fit[j], 
                         sd = pred_log$se.fit[j]))
    expected_acm_samples[whichs_pred[j], ] <- samples
    
    samples_nb <- rnbinom(num_samples, size = overd, mu = samples)
    expected_acm_samples_nb[whichs_pred[j], ] <- samples_nb
  }
  
  # Get fitted values for historical time periods
  years <- unique(temp$year)
  range_years <- min(years):max(years)
  pred_hist <- predict(annual_model, se.fit = TRUE, type = "response",
                       newdata = data.frame(year = range_years))
  acm.byyear[whichs, "expected_acm"] <- pred_hist$fit
  acm.byyear[whichs, "expected_acm_se"] <- pred_hist$se.fit
  pred_log_hist <- predict(annual_model, se.fit = TRUE,
                           newdata = data.frame(year = range_years))
  acm.byyear[whichs, "expected_log_acm"] <- pred_log_hist$fit
  acm.byyear[whichs, "expected_log_acm_se"] <- pred_log_hist$se.fit
}

#### Save all annual expecteds, i.e. for pre 2020 and 2020-2021 ####
acm_annual_predictions_tier2_hist <- acm.byyear
save(acm_annual_predictions_tier2_hist, 
     file = "../Generated_Data/acm_annual_predictions_tier2_hist.RData")

#### Save annual expecteds for 2020-2021 ####
acm_annual_predictions_tier2 <- acm_predictions
save(acm_annual_predictions_tier2, 
     file = "../Generated_Data/acm_annual_predictions_tier2.RData")

#### Save samples of annual expecteds for 2020-2021 ####
save(expected_acm_samples, expected_acm_samples_nb, 
     file = "../Generated_Data/expected_acm_samples_tier2.RData")
