#' @param model      : INLA model object 
#' @param df.a       : the data frame used to fit the INLA object
#' @param model.lab  : A label for the model used 
#' @param obsd       : Boolean "yes" or "no" to replace fit with observed values
#' @param ci         : The interval to use in summarising outputs, e.g. 95% or 80%
#' @return           : A data frame containing all relevant data for visualisation or summary
#' 

rm(list = ls(all = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pacman::p_load(data.table, lubridate, tidyr, dplyr, INLA, xlsx)

load("mf.df.Rda")

get.intervals <- function(model.lab){
  
  mf.df <- mf.df %>% mutate(months = month + 12*(year - 2020)) %>%
    arrange(Country, months)
  
#########################################################################
# Function to get mean, lower and upper CI
  get.dist <- function(x){
    fit  <- mean(x, na.rm = T)
    lwr  <- quantile(x, 0.025, na.rm = T)
    uppr <- quantile(x, 0.975, na.rm = T)
    c(fit, lwr, uppr)
  }
  
#########################################################################
  mk.mat <- function(df){
    df %>% mutate(Country = ifelse(Country == "CÃ´te d'Ivoire","Côte d'Ivoire", Country)) %>% 
      arrange(Country, months) %>% select(-c("Country", "WHO_region", "months")) %>% as.matrix()
  }
  
exp.samps        <- fread("expected.csv") %>% mk.mat()
ac.samps         <- fread("acm.csv") %>% mk.mat()
excess.samps     <- ac.samps - exp.samps
pscore.samps     <- 100*excess.samps/exp.samps

rep.covid    <- mf.df %>% select(covid_deaths) %>% as.matrix()
observed.df  <- mf.df %>% select(observed) %>% as.matrix()
mx.m         <- 24
########################################################################
  
#########################################################################
# Now get distributions
#########################################################################
print("Distributions of monthly mortality predictions")
  
  # All cause deaths, expected and covid
  # 1) Get distribution by country
  print(".................................By Country")
  
  country.ac            <- t(apply(ac.samps, 1, get.dist)) # Predicted all cause deaths
  country.exp           <- t(apply(exp.samps, 1, get.dist)) # Expected all cause deaths
  country.exc           <- t(apply(excess.samps, 1, get.dist)) # Excess deaths
  country.pscore        <- t(apply(pscore.samps, 1, get.dist)) # pscore
  country.covid         <- rep.covid[,1]
  
  # 2 Get distribution for global 
  print(".................................By Global")
  global.pscore <- global.exc <- global.exp <- global.ac  <- array(dim = c(mx.m, 3)) # empty matrix to populate with distribution
  global.covid          <- rep(NA, mx.m)
  for (m in 1:mx.m){
    # get row indices for month m for entire matrix of posterior draws 
    rows  <- mf.df %>% filter(months == m) %>% pull(rowid) 
    # first sum by draw for the selected rows and then get distribution
    global.ac[m,]      <- get.dist(apply(ac.samps[rows,], 2, sum, na.rm = T)) 
    global.exp[m,]     <- get.dist(apply(exp.samps[rows,], 2, sum, na.rm = T)) 
    global.exc[m,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T)) 
    global.pscore[m,]  <- get.dist(100*apply(excess.samps[rows,], 2, sum, na.rm = T)/
                                     apply(exp.samps[rows,], 2, sum, na.rm = T)) 
    global.covid[m]    <- sum(rep.covid[rows,1]) 
  }
  
  # 3) by region
  print(".................................By Region")
  region.pscore <- region.exc <- region.exp <- region.ac <- array(dim = c(mx.m*6, 3))
  region.covid         <- rep(NA, mx.m*6)
  r = 1
  for (reg in c("AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO")){
    for (m in 1:mx.m){
      # get row indices for month m and region reg for entire matrix of posterior draws 
      rows   <- mf.df %>% filter(WHO_region == reg & months == m) %>% pull(rowid)
      # first sum by draw for he selected rows and then get distribution
      region.ac[r,]      <- get.dist(apply(ac.samps[rows,], 2, sum, na.rm = T))
      region.exp[r,]     <- get.dist(apply(exp.samps[rows,], 2, sum, na.rm = T))
      region.exc[r,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
      region.pscore[r,]  <- get.dist(100*apply(excess.samps[rows,], 2, sum, na.rm = T)/
                                       apply(exp.samps[rows,], 2, sum, na.rm = T))
      region.covid[r]    <- sum(rep.covid[rows,1]) 
      r = r + 1
    }
  }

  # 4) by income group
  print(".................................By Income Group")
  income.pscore <- income.exc <- income.exp <- income.ac  <- array(dim = c(mx.m*4, 3))
  income.covid         <- rep(NA, mx.m*4)
  r = 1
  for (reg in c("LIC", "LMIC", "UMIC", "HIC")){
    for (m in 1:mx.m){
      # get row indices for month m and income reg for entire matrix of posterior draws 
      rows   <- mf.df %>% filter(income_group == reg & months == m) %>% pull(rowid)
      # first sum by draw for he selected rows and then get distribution
      income.ac[r,]      <- get.dist(apply(ac.samps[rows,], 2, sum, na.rm = T))
      income.exp[r,]     <- get.dist(apply(exp.samps[rows,], 2, sum, na.rm = T))
      income.exc[r,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
      income.pscore[r,]     <- get.dist(100*apply(excess.samps[rows,], 2, sum, na.rm = T)/
                                       apply(exp.samps[rows,], 2, sum, na.rm = T))
      income.covid[r]    <- sum(rep.covid[rows,1]) 
      r = r + 1
    }
  }
  
  ####################################################################################
  
  
  colnames(income.ac)<-colnames(region.ac)<-colnames(global.ac)<-colnames(country.ac)<- c("ac.fit", "ac.lwr", "ac.uppr")
  colnames(income.exp)<-colnames(region.exp)<-colnames(global.exp)<-colnames(country.exp)<- c("exp.fit", "exp.lwr", "exp.uppr")
  colnames(income.exc)<-colnames(region.exc)<-colnames(global.exc)<-colnames(country.exc)<-  c("fit", "lwr", "uppr")
  colnames(income.pscore)<-colnames(region.pscore)<-colnames(global.pscore)<-colnames(country.pscore)<-  c("psc.fit", "psc.lwr", "psc.uppr")
  
  ####################################################################################
  all.covid    <- c(country.covid, global.covid, region.covid, income.covid)
  all.expected <- rbind(country.exp, global.exp, region.exp, income.exp)  
  all.allcause <- rbind(country.ac, global.ac, region.ac, income.ac)
  all.excess   <- rbind(country.exc, global.exc, region.exc, income.exc)  
  all.pscore   <- rbind(country.pscore, global.pscore, region.pscore, income.pscore)  
  
  ####################################################################################
  
  
  ####################################################################################
  # CUMULATIVE EXCESS and COVID just for year 2021
  # Similar approach but instead of pulling rows for month m, it is month m and earlier 
  # 1) by country
  print("Cumulative excess estimates for year 2021")
  print(".................................By Country")
  isos                    <- sort(unique(mf.df$iso3))
  country.ex.c21          <- array(dim = c(mx.m*length(isos), 3))
  country.covid.c21       <- rep(NA, mx.m*length(isos))
  r = 1
  for (reg in isos){
    for (m in 1:mx.m){
      if (m <= 12){
        country.ex.c21[r,]    <- 0
        country.covid.c21[r]  <- 0
      } 
      if (m == 13){
        rows                  <- mf.df %>% filter(year == 2021 & months <= m & iso3 == reg) %>% pull(rowid)
        country.ex.c21[r,]    <- get.dist(excess.samps[rows,])
        country.covid.c21[r]  <- sum(rep.covid[rows,1])
      }
      
      if (m >= 14){
        rows                  <- mf.df %>% filter(year == 2021 & months <= m & iso3 == reg) %>% pull(rowid)
        country.ex.c21[r,]    <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
        country.covid.c21[r]  <- sum(rep.covid[rows,1])
      }
      r = r + 1
    }
  }
  
  # 2 Global
  print(".................................By Global")
  global.ex.c21             <- array(dim = c(mx.m, 3))
  global.covid.c21          <- rep(NA, mx.m)
  for (m in 1:mx.m){
    if (m <= 12){
      global.ex.c21[m,]       <- 0
      global.covid.c21[m]     <- 0
    } else {
    rows          <- mf.df %>% filter(months <= m & year == 2021) %>% pull(rowid)
    global.ex.c21[m,]       <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
    global.covid.c21[m]     <- sum(rep.covid[rows,1])
    }
  }
  
  # 3) by region
  print(".................................By Region")
  region.ex.c21             <- array(dim = c(mx.m*6, 3))
  region.covid.c21          <- rep(NA, mx.m*6)
  r = 1
  for (reg in c("AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO")){
    for (m in 1:mx.m){
      if (m <= 12){
        region.ex.c21[r,]     <- 0
        region.covid.c21[r]   <- 0
      } else  {
        rows                  <- mf.df %>% filter(year == 2021 & months <= m & WHO_region == reg) %>% pull(rowid)
        region.ex.c21[r,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
        region.covid.c21[r]   <- sum(rep.covid[rows,1])
      }     
      r = r + 1
    }
  }
  
  # 4) by income
  print(".................................By income group")
  income.ex.c21             <- array(dim = c(mx.m*4, 3))
  income.covid.c21          <- rep(NA, mx.m*4)
  r = 1
  for (reg in c("LIC", "LMIC", "UMIC", "HIC")){
    for (m in 1:mx.m){
      if (m <= 12){
      income.ex.c21[r,]     <- 0
      income.covid.c21[r]   <- 0
      } else {
      rows                <- mf.df %>% filter(year == 2021 & months <= m & income_group == reg) %>% pull(rowid)
      income.ex.c21[r,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
      income.covid.c21[r]   <- sum(rep.covid[rows,1])  
      }
      r = r + 1
    }
  }
  
  colnames(income.ex.c21)<-colnames(region.ex.c21)<-colnames(global.ex.c21)<-colnames(country.ex.c21)  <- c("c.fit21", "c.lwr21", "c.uppr21")
  all.excess.c21         <- rbind(country.ex.c21, global.ex.c21, region.ex.c21, income.ex.c21) 
  all.covid.c21          <- c(country.covid.c21, global.covid.c21, region.covid.c21, income.covid.c21)
  
  
  ####################################################################################
  # CUMULATIVE EXCESS and COVID
  # Similar approach but instead of pulling rows for month m, it is month m and earlier 
  # 1) by country
  print("Cumulative excess estimates for entire period")
  print(".................................By Country")
  isos                    <- sort(unique(mf.df$iso3))
  country.ex.c            <- array(dim = c(mx.m*length(isos), 3))
  country.covid.c         <- rep(NA, mx.m*length(isos))
  r = 1
  for (reg in isos){
    for (m in 1:mx.m){
      rows                <- mf.df %>% filter(months <= m & iso3 == reg) %>% pull(rowid)
      if (m == 1){
        country.ex.c[r,]    <- get.dist(excess.samps[rows,])
      } else {
        country.ex.c[r,]    <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
      }
      country.covid.c[r]  <- sum(rep.covid[rows,1])
      r = r + 1
    }
  }
  
  # 2 Global
  print(".................................By Global")
  global.ex.c             <- array(dim = c(mx.m, 3))
  global.covid.c          <- rep(NA, mx.m)
  for (m in 1:mx.m){
    rows          <- mf.df %>% filter(months <= m) %>% pull(rowid)
    global.ex.c[m,]       <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
    global.covid.c[m]     <- sum(rep.covid[rows,1])
  }
  
  # 3) by region
  print(".................................By Region")
  region.ex.c             <- array(dim = c(mx.m*6, 3))
  region.covid.c          <- rep(NA, mx.m*6)
  r = 1
  for (reg in c("AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO")){
    for (m in 1:mx.m){
      rows                <- mf.df %>% filter(months <= m & WHO_region == reg) %>% pull(rowid)
      region.ex.c[r,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
      region.covid.c[r]   <- sum(rep.covid[rows,1])
      r = r + 1
    }
  }
  
  # 4) by income
  print(".................................By income group")
  income.ex.c             <- array(dim = c(mx.m*4, 3))
  income.covid.c          <- rep(NA, mx.m*4)
  r = 1
  for (reg in c("LIC", "LMIC", "UMIC", "HIC")){
    for (m in 1:mx.m){
      rows                <- mf.df %>% filter(months <= m & income_group == reg) %>% pull(rowid)
      income.ex.c[r,]     <- get.dist(apply(excess.samps[rows,], 2, sum, na.rm = T))
      income.covid.c[r]   <- sum(rep.covid[rows,1])
      r = r + 1
    }
  }
  
  colnames(income.ex.c)<-colnames(region.ex.c)<-colnames(global.ex.c)<-colnames(country.ex.c)  <- c("c.fit", "c.lwr", "c.uppr")
  all.excess.c         <- rbind(country.ex.c, global.ex.c, region.ex.c, income.ex.c) 
  all.covid.c          <- c(country.covid.c, global.covid.c, region.covid.c, income.covid.c)
  
  ############################################################################################################
  # Labels for the final data frame
  full.labs <- mf.df %>%
    select(iso3, months) %>% rename(location = iso3) %>%
    rbind(data.table(location = "Global", months = 1:mx.m)) %>%
    rbind(data.table(location = rep(c("AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO"), each = mx.m), 
                     months = rep(1:mx.m, 6))) %>%
    rbind(data.table(location = rep(c("LIC", "LMIC", "UMIC", "HIC"), each = mx.m), 
                     months = rep(1:mx.m, 4)))
  
output.table <- data.table(full.labs, 
             use.observed  = "yes",
             covid         = all.covid,
             c.covid       = all.covid.c,
             all.expected,
             all.allcause,
             all.excess,
             all.excess.c,
             all.pscore, 
             c.covid21    = all.covid.c21,
             all.excess.c21,
             model.name    = model.lab) %>% 
    mutate(year  = ifelse(months <= 12, 2020, 
                          ifelse(months <= 24 ,2021, 2022)),
           month = months - 12*(year - 2020)) %>%
    rename(expected = exp.fit)

return(output.table)
}

loc.labs <- mf.df %>% select(Country, iso3) %>% unique() %>%
  rename(location = iso3) %>%
  rbind(data.table(Country = c("Global", "AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO", 
                               "LIC", "LMIC", "UMIC", "HIC"),
                   location = c("Global", "AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO",
                                "LIC", "LMIC", "UMIC", "HIC")))

who.mod.ests <- get.intervals("Model All") %>% left_join(loc.labs, by = "location")

##################################################################################################
# 

write.xlsx(who.mod.ests,row.names = FALSE, 
           file="Estimate.Database.xlsx",
           sheetName="Estimates", append=FALSE)

# Compare models with IHME and economist

ihme2411 <- fread("https://www.dropbox.com/s/v9twekf68eecf8s/ihme1811.csv?dl=1") %>%
  filter(! WHO_region %in% c("China", "India"))

econ2411 <- rbind(fread("https://raw.github.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/output-by-alternative-regions/export_regions_WHO_region_cumulative.csv")  %>% 
                    select(region, date, cumulative_estimated_daily_excess_deaths,
                           cumulative_estimated_daily_excess_deaths_ci_95_bot,
                           cumulative_estimated_daily_excess_deaths_ci_95_top),
                  fread("https://raw.github.com/TheEconomist/covid-19-the-economist-global-excess-deaths-model/main/output-data/export_world_cumulative.csv") %>% rename(region = world) %>%
                    select(region, date, cumulative_estimated_daily_excess_deaths,
                           cumulative_estimated_daily_excess_deaths_ci_95_bot,
                           cumulative_estimated_daily_excess_deaths_ci_95_top))  %>%
  rename(Country = region, 
         c.fit  = cumulative_estimated_daily_excess_deaths, 
         c.lwr  = cumulative_estimated_daily_excess_deaths_ci_95_bot,
         c.uppr = cumulative_estimated_daily_excess_deaths_ci_95_top) %>%
  mutate(date  = as.Date(date), 
         year   =  lubridate::year(lubridate::ymd(date)),
         month  =  lubridate::month(lubridate::ymd(date)), 
         months = ifelse(year == 2020, month, month + 12), 
         Country = ifelse(Country == "World", "Global", Country)) %>%
  filter(Country %in% ihme2411$WHO_region & year < 2022) %>%
  group_by(Country, year) %>% filter(date == max(date)) %>% ungroup() %>% 
  select(Country, months, c.fit, c.lwr, c.uppr) %>% 
  rename(WHO_region = Country, period = months, mean = c.fit, lwr = c.lwr, uppr = c.uppr) %>% 
  mutate(source = "Economist")


who2411 <- who.mod.ests %>%
  filter(model.name == "Model All" &
           use.observed  == "yes") %>% 
  filter(Country %in% ihme2411$WHO_region & months %in% c(12,24)) %>% 
  select(Country, months, c.fit, c.lwr, c.uppr) %>% 
  rename(WHO_region = Country, period = months, mean = c.fit, lwr = c.lwr, uppr = c.uppr) %>% 
  mutate(source = "WHO")

reg.ests <- rbind(ihme2411, econ2411, who2411) %>% 
  mutate(WHO_region = factor(WHO_region, levels = c("Global", "AFRO", "AMRO", "EMRO", "EURO", "SEARO", "WPRO"))) %>%
  arrange(period, WHO_region, source) 

write.xlsx(reg.ests,row.names = FALSE, 
           file="Model.Comparison.xlsx",
           sheetName="Estimates", append=FALSE)
