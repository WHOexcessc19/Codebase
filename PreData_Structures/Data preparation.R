#####################################################################################################################################
# House keeping
rm(list = ls(all = TRUE))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(data.table, lubridate, tidyr, dplyr, countrycode, ggplot2, tm, mgcv, httr,
               prophet, tsibble, fable.prophet, fable, wktmo)

select    <- dplyr::select

#####################################################################################################################################
# Additional files for analysis
#####################################################################################################################################

# 1) Codes for mapping iso2 to iso3 for WHO public data
# 2) Total population numbers by country, age and sex for year 20001 to 2021 (WPP2019 for 194 Member States)
# 3) GHE data by broad cause group, age, sex and year 2000 to 2019 (183/194 member states)
# 4) GBD data by broad cause group, age, sex and year 2000 to 2019 (11/194 member states)
# *3 and 4 are combined in "ghe.all.causes".

ccodes         <- fread("ccodes.csv")
population     <- fread("population.csv")
ghe.all.causes <- fread("ghe.all.causes.csv")

#####################################################################################################################################
# 5) WHO COVID case and death data by date and country

# Read in the WHO data from https://covid19.who.int/WHO-COVID-19-global-data.csv and clean up some of the names
# The version used in the analysis is from the 24th of March 2022. Commented out cleaning, final version loaded

# who.covid      <- fread("https://covid19.who.int/WHO-COVID-19-global-data.csv") %>%
#   filter(Country != "Other") %>%
#   mutate(Country_code = ifelse(Country == "Namibia", "NA", Country_code)) %>%
#   left_join(ccodes, by = "Country_code") %>%
#   mutate(Country = case_when(Country == "Curaçao" ~ "Curacao",
#                              Country == "Saint Barthélemy" ~ "Saint Barthelemy",
#                              Country == "Réunion" ~ "Reunion",
#                              Country == "Kosovo[1]" ~ "Kosovo",
#                              Country %in% c("Bonaire","Saba","Sint Eustatius") ~ "Bonaire, Sint Eustatius and Saba",
#                              TRUE ~ as.character(Country)),
#          iso3    = case_when(Country == "Kosovo" ~ "XKX",
#                              Country == "Saint Martin"~ "MAF",
#                              Country == "Sint Maarten"~ "SXM",
#                              Country == "Bonaire, Sint Eustatius and Saba" ~ "BES",
#                              Country == "Curacao" ~ "CUW",
#                              Country == "Saint Barthelemy" ~ "BLM",
#                              Country == "South Sudan" ~ "SSD",
#                              Country == "Namibia" ~ "NAM",
#                              TRUE ~ as.character(iso3)),
#          date    = as.Date(Date_reported,"%y-%m-%d"),
#          WHO_region = factor(WHO_region, levels = c("AFRO","AMRO","EMRO","EURO","SEARO","WPRO"))) %>%
#   mutate(year   =  lubridate::year(lubridate::ymd(date)),
#          month  =  lubridate::month(lubridate::ymd(date))) %>%
#   group_by(Country, iso3, WHO_region, year, month) %>%
#   summarise(covid_cases  = sum(New_cases, na.rm = T),
#             covid_deaths = sum(New_deaths, na.rm = T), .groups = "drop") %>% ungroup() %>%
#   mutate(covid_cases = ifelse(covid_cases < 0, 0, round(covid_cases)),
#          covid_deaths = ifelse(covid_deaths < 0, 0, round(covid_deaths))) %>%
#   filter(year < 2022)
#   fwrite(who.covid, file = "clean.who.covid.csv")

who.covid         <- fread("clean.who.covid.csv")

#####################################################################################################################################
#####################################################################################################################################
# 6)  WHO all cause mortality death data
#####################################################################################################################################
#####################################################################################################################################

who.allcause <- fread("WHO_Allcause_Mortality_Data.csv") %>%
  mutate(sex = ifelse(sex == "", "Unknown", sex),
         time = ifelse(time_unit == "Week" & time == 53, 12,
                       ifelse(time_unit == "Week" & time < 53,
                              lubridate::month(as.Date(paste0(year,"-",time,"-",1),format="%Y-%U-%u")), time)),
         time_unit = ifelse(time_unit== "Week", "Month", time_unit),
         age_cat_s = ifelse(age_cat_s %in% c("0", "1-4", "0-4", "< 5"), "0-4",
                            ifelse(age_cat_s == "total age", "total",
                                   ifelse(age_cat_s %in% c("85-90", "90+", "85-89", "90-94", "95-100", "95+", "100+"),
                                          "85+", age_cat_s)))) %>%
  group_by(country, sex, year, time, time_unit, source, format_age, age_cat_s, coverage) %>%
  summarise(deaths = sum(deaths, na.rm = T), .groups = "drop") %>% ungroup() %>%
  right_join(data.frame(age_cat_s = c(paste0(seq(0,80,5),"-",seq(0,80,5)+4), "85+", "total"),
                        age = c(seq(0,85,5), 999)), by = "age_cat_s") %>%
  spread(sex, deaths) %>%  rename(iso3 = country) %>%
  rowwise() %>% mutate(Both = sum(c(Female, Male, Unknown), na.rm = T)) %>% ungroup() %>%
  mutate(source = "WHO Data Call") %>% rename(age_cat = age_cat_s) %>%
  select(source, coverage, iso3, year, time, age, age_cat, Female, Male, Both) %>%
  gather(sex, deaths, -source, -coverage, -iso3, -year, -time, -age, -age_cat) %>%
  arrange(iso3, sex, year, time, age) %>% rename(month = time)

#####################################################################################################################################
#####################################################################################################################################
# 7) STMF Short term mortality fluctations

# Read in the STMF and clean up some of the names
# The version used in the analysis is from the 24th of March 2022. Commented out cleaning, final version loaded
#####################################################################################################################################
#####################################################################################################################################

# stmf      <- fread("https://www.mortality.org/File/GetDocument/Public/STMF/Outputs/stmf.csv", skip = 2) %>%
#   select(CountryCode, Sex, Year, Week, D0_14, D15_64, D65_74, D75_84, D85p, DTotal) %>%
#   rename(year = Year, iso3 = CountryCode, sex = Sex,
#          `0-14`=D0_14,`15-64`=D15_64,`65-74`=D65_74,`75-84`=D75_84,`85+`=D85p, `total`=DTotal) %>%
#   mutate(sex = ifelse(sex=="m","Male", ifelse(sex=="f","Female","Both")),
#          time = ifelse(Week == 53, 12,
#                        lubridate::month(as.Date(paste0(year,"-",Week,"-",1),format="%Y-%U-%u"))),
#          iso3 = ifelse(iso3 %in% c("GBR_NIR","GBR_SCO","GBRTENW"),"GBR",
#                        ifelse(iso3 == "DEUTNP", "DEU",
#                               ifelse(iso3 == "FRATNP", "FRA",
#                                      ifelse(iso3 == "AUS2", "AUS",
#                                             ifelse(iso3 == "NZL_NP", "NZL", iso3)))))) %>%
#   filter(year > 2010) %>%
#   gather(age_cat, deaths, -iso3, -sex, -year, -time, -Week)
# 
# stmf <- stmf %>% mutate(time = 999) %>% rbind(stmf) %>%
#   group_by(age_cat, iso3, sex, year, time) %>%
#   summarise(deaths = sum(deaths, na.rm = T), .groups = "drop") %>% ungroup() %>%
#   right_join(data.frame(age_cat = c("0-14", "15-64", "65-74", "75-84","85+", "total"),
#                         age = c(0,15,65,75,85,999)), by = "age_cat") %>% rename(month = time) %>%
#   mutate(source = "STMF", coverage = "National") %>% select(names(who.allcause))
# 
# fwrite(stmf, file = "stmf.clean.csv")

stmf   <- fread("stmf.clean.csv")

#########################################################################################################################
#########################################################################################################################
# 8) EUROSTAT data by age and sex from 
# https://ec.europa.eu/eurostat/estat-navtree-portlet-prod/BulkDownloadListing?file=data/demo_r_mwk_05.tsv.gz
#########################################################################################################################

eurostat.df <- fread("eurostat.csv")

eurostat    <- eurostat.df %>%
  separate(week,c("year","week"), "W") %>%
  right_join(data.frame(ageg = c("Less than 5 years",
                                 paste0("From ", seq(5,85,5)," to ",seq(5,85,5)+4, " years"), "90 years or over", "Total"),
                        age = c(seq(0,85,5), 85, 999),
                        age_cat = c(paste0(seq(0,80,5),"-",seq(0,80,5)+4), "85+", "85+", "total")), by = "ageg") %>%
  filter(!is.na(sex) & val != ":") %>% spread(sex, val) %>%
  rename(Female = Females, Male = Males) %>%
  mutate(Female = as.numeric(Female), Male = as.numeric(Male)) %>%
  filter(!(is.na(Female)|is.na(Male))) %>%
  mutate(Both = Female + Male,
         iso3 = countrycode::countrycode(country, "country.name", "iso3c"),
         year = as.numeric(year), week = as.numeric(week),
         time = ifelse(week >= 53, 12,
                       lubridate::month(as.Date(paste0(year,"-",week,"-",1),format="%Y-%U-%u")))) %>%
  group_by(iso3, year, time, age, age_cat) %>%
  summarize(Female = sum(Female, na.rm = T), Male = sum(Male, na.rm = T),
            Both   = sum(Both, na.rm = T), .groups = "drop") %>% ungroup() %>%
  filter(year < 2022)

eurostat  <- eurostat %>% mutate(time = 999) %>%
  group_by(iso3, year, time, age, age_cat) %>%
  summarize(Female = sum(Female, na.rm = T), Male = sum(Male, na.rm = T),
            Both   = sum(Both, na.rm = T), .groups = "drop") %>% ungroup() %>%
  rbind(eurostat) %>%
  mutate(source = "Eurostat", coverage = "National") %>%
  select(source, coverage, iso3, year, time, age, age_cat, Female, Male, Both) %>%
  gather(sex, deaths, -source, -coverage, -iso3, -year, -time, -age, -age_cat) %>%
  arrange(iso3, sex, year, time, age) %>% rename(month = time)

#####################################################################################################################################
#####################################################################################################################################
# 9) WMD: Observed and expected data by week from Karlinsky et al and multiple sources (world mortality data)

# baselines.csv file from "https://raw.github.com/dkobak/excess-mortality/main/baselines.csv"
# world_mortality.csv file from "https://raw.github.com/akarlinsky/world_mortality/main/world_mortality.csv"
# Loading clean and monthly aggregated as of the 24th March 2022

wmd.clean <- fread("wmd.clean.csv")

#####################################################################################################################################
# 10) Expected deaths

load("../Generated_Data/acm_monthly_predictions_tier1.RData")
load("../Generated_Data/acm_monthly_predictions_tier2.RData")

expected.out     <- rbind(acm_monthly_predictions_tier1 %>% select(iso3, year, month, expected_acm, expected_acm_se),
                          acm_monthly_predictions_tier2 %>% select(iso3, year, month, expected_acm, expected_acm_se) %>% 
                            filter(!iso3 %in% acm_monthly_predictions_tier1$iso3)) %>%
  rename(expected = expected_acm, se = expected_acm_se) %>%
  mutate(expected.lwr = expected - 2*se, expected.uppr = expected + 2*se, se = NULL)

rm(list = c("acm_monthly_predictions_tier1", "acm_monthly_predictions_tier2"))

wmd.df <- wmd.clean %>%
  right_join(expected.out, by = c("iso3", "month", "year"))

#####################################################################################################################################
#####################################################################################################################################
# COVARIATES
#####################################################################################################################################
#####################################################################################################################################

#####################################################################################################################################
# 11) Oxford tracker government policy ("https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_nat_latest.csv")
# Using version from 24 March 2022

# covid_pol    <- fread("https://github.com/OxCGRT/covid-policy-tracker/raw/master/data/OxCGRT_nat_latest.csv")  %>%
#   rename(iso3 = CountryCode) %>%
#   mutate(year   =  lubridate::year(lubridate::ymd(Date)),
#          month  =  lubridate::month(lubridate::ymd(Date))) %>%
#   group_by(iso3, month, year) %>%
#   summarize(Stringency = median(StringencyIndex_Average, na.rm = T),
#             Government = median(GovernmentResponseIndex_Average, na.rm = T),
#             Containment = median(ContainmentHealthIndex_Average, na.rm = T),
#             Economic = mean(EconomicSupportIndex, na.rm = T), .groups = "drop") %>%
#   ungroup() %>% filter(!is.na(Stringency))
# fwrite(covid_pol, file = "covid_pol.clean.csv")

covid_pol <- fread("covid_pol.clean.csv")

#####################################################################################################################################
# 12) SDI, filter for country locations (leave out region or subnational; & Georgia in the USA)

sdi.loc <- fread("sdi.country.csv")

sdi.gbd <- readxl::read_excel("IHME_GBD_2019_SDI_1990_2019_Y2020M10D15.XLSX", skip = 1)[-c(105), ] %>%
  filter(Location %in% sdi.loc$Location & Location != "Virgin Islands" ) %>%
  mutate(iso3 = countrycode(Location, "country.name", "iso3c")) %>%
  gather(year, sdi, -Location, -iso3) %>%
  mutate(sdi=as.numeric(gsub( "?", ".",sdi)), year = as.numeric(year), Location = NULL) %>%
  filter(year == 2019 & iso3 %in% unique(population$iso3)) %>%
  select(-c(year))

#####################################################################################################################################
# 13) Income group from World Bank, assume Cook Islands and Niue are upper middle income like Tonga and Tuvalu

wbincome <- fread("WBG.csv") %>%
  select(Code, `Income group`) %>%
  rename(iso3 = Code, group = "Income group") %>%
  filter(iso3 %in% unique(population$iso3))

wbincome <- wbincome  %>% filter(iso3 %in% c("TON","TUV")) %>% mutate(iso3 = c("COK", "NIU")) %>%
  rbind(wbincome) %>%
  mutate(income_group = ifelse(group == "Low income", "LIC",
                               ifelse(group == "Lower middle income", "LMIC",
                                      ifelse(group == "Upper middle income", "UMIC","HIC"))))

inc.groups <- wbincome %>%
  mutate(high.income  = ifelse(group == "High income", 1, 0),
         upper.income = ifelse(group %in% c("Upper middle income","High income"), 1, 0),
         group = NULL)
#####################################################################################################################################
# 14) Our world in Data (test positivity rate, population density and life-expectancy)
# Taken from https://covid.ourworldindata.org/data/owid-covid-data.csv
# Version date 24 March 2022

# owid <- fread("https://covid.ourworldindata.org/data/owid-covid-data.csv") %>%
#   rename(iso3 = iso_code) %>%
#   mutate(date = as.Date(date, "%Y-%m-%d")) %>%
#   select(iso3, date, gdp_per_capita, population_density, life_expectancy, positive_rate, new_tests, new_cases) %>%
#   mutate(year   =  lubridate::year(lubridate::ymd(date)),
#          month  =  lubridate::month(lubridate::ymd(date)),
#          date = NULL) %>%
#   group_by(iso3, year, month) %>%
#   summarise(population_density      = median(population_density, na.rm = T),
#             life_expectancy         = median(life_expectancy, na.rm = T),
#             positive_rate           = median(positive_rate, na.rm = T),
#        .groups = "drop") %>% ungroup()
# fwrite(owid, file = "owid.clean.csv")

owid <- fread("owid.clean.csv")

#####################################################################################################################################
# 15) GBD2019 - by age and sex for 2019 - diabetes prevalence, cardiovascular, ncd and hiv mortality

gbd_in       <- fread("gbd_input.csv") %>%
  rename(ages = age) %>%
  mutate(iso3 = countrycode::countrycode(location, "country.name", "iso3c"))

gbd_in_diab <- gbd_in %>%
  filter(ages == "All Ages" & cause == "Diabetes mellitus" & measure == "Prevalence") %>%
  group_by(iso3) %>% summarise(val = mean(val)) %>% ungroup() %>%
  rename(diabetes_allages = val)

gbd_in_oth <- gbd_in %>%
  filter(ages != "All Ages" & cause != "Diabetes mellitus" & measure != "Prevalence") %>%
  select(iso3, sex, ages, cause, val) %>% spread(cause, val) %>%
  rename(cardiovascular = "Cardiovascular diseases", hiv = "HIV/AIDS", ncds = "Non-communicable diseases") %>%
  left_join(data.frame(ages = c("Under 5", paste0(seq(5,80,5), " to ",seq(5,80,5)+4 ), "85 plus"),
                       age  = seq(0,85,5)), by = "ages") %>%
  mutate(ages = NULL) %>% arrange(iso3, sex, age) %>%
  right_join(population %>% filter(year == 2019) %>% select(iso3, age, sex, Nx), by = c("iso3","age","sex")) %>%
  mutate(cardiovascular = cardiovascular*Nx/1e5, hiv = hiv*Nx/1e5, ncds = ncds*Nx/1e5) %>%
  gather(cause, val, -iso3, -sex, -age) %>%
  spread(sex, val) %>%
  rowwise() %>% mutate(Both = sum(c(Female, Male), na.rm = T)) %>% ungroup() %>%
  gather(sex, val, -iso3, -cause, -age)

gbd_in_oth <- gbd_in_oth %>%
  group_by(iso3, sex, cause) %>%
  summarise(val = sum(val, na.rm = T), .groups = "drop") %>%
  ungroup() %>% mutate(age = 999) %>%
  rbind(gbd_in_oth) %>%
  spread(cause, val) %>%
  mutate(cardiovascular = cardiovascular/Nx*1e5,
         hiv = hiv/Nx*1e5, ncds = ncds/Nx*1e5, Nx = NULL)

#####################################################################################################################################
# 16) WPP2019 - proportion of the population over 65 and proportion of the population under 15

pop.sum <- population %>% filter(year == 2019) %>%
  group_by(iso3, age) %>% summarise(Nx = sum(Nx), .groups = "drop") %>%
  mutate(group = ifelse(age <= 15, "under15",
                        ifelse(age > 65, "over65", "other"))) %>%
  group_by(iso3) %>% mutate(totp = sum(Nx)) %>% ungroup() %>%
  group_by(iso3, group, totp) %>% summarise(totg = sum(Nx), .groups = "drop") %>% ungroup() %>%
  mutate(prop = totg/totp, totp = NULL, totg = NULL) %>%
  spread(group, prop) %>% select(-c(other))

#####################################################################################################################################
# 17) Temparature by country and month

monthsl         <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

temp_celsius    <- fread("temparature.csv") %>%
  filter(Year >= 2015) %>% 
  rename(year = Year, temperature = "Temperature - (Celsius)") %>%
  mutate(iso3 = countrycode::countrycode(Country, "country.name", "iso3c")) %>%
  mutate(month = removeWords(Statistics, " Average")) %>%
  mutate(month = as.numeric(factor(month, levels = monthsl))) %>%
  select(iso3, year, month, temperature) %>% arrange(iso3, year, month)

temperatures    <- temp_celsius %>%
  group_by(iso3, month) %>%
  summarise(temperature = mean(temperature), .groups = "drop") %>%
  ungroup() %>%  arrange(iso3, month) %>% select(iso3, month, temperature)  %>%
  group_by(iso3) %>%
  mutate(ztemp = scale(temperature)) %>% ungroup()

# Also the HDI covariate
hdi           <- fread("hdi.csv") %>%
  mutate(iso3 = countrycode::countrycode(country, "country.name", "iso3c"), hdi = 100*hdi) %>%
  select(-c(country))

#####################################################################################################################################
#####################################################################################################################################
# 18) Bring the covariates data all together
# Function to replace with medians by group
#####################################################################################################################################
#####################################################################################################################################

rep_with_median <- function(df, inds){
  for (j in inds){
    df <- df %>% rename(cvar = j) %>%
      mutate(cvar = ifelse(is.na(cvar), median(cvar, na.rm = T), cvar))
    data.table::setnames(df, "cvar", j)
  }
  df
}

###########################################################################################
# Also work out things like log_u5mr, 45q15, 40q30 and e60 to include in mf.df
# function to estimate lifetable elements

est.lt.elem <- function(mx){
  nx  <- c(rep(5, 17), 1)

  qx  <- c(1 - exp(-nx * mx)); qx[18] <- 1
  ax  <- (nx + 1/mx - nx/qx)
  a   <- rep(0.5, 17)
  lx  <- c(1, cumprod(1 - qx)[1:17])
  dx  <- c(rev(diff(rev(lx))), lx[1] - sum(rev(diff(rev(lx)))))
  Lx  <- nx * lx - (nx - ax) * dx
  Tx  <- rev(cumsum(rev(Lx)))
  ex  <- Tx/lx

  q5      <- 1000*(lx[2] - lx[1])/lx[1]
  q15     <- 1000*(lx[4] - lx[13])/lx[4]
  q30     <- 1000*(lx[7] - lx[15])/lx[7]
  e60     <- ex[13]
  e0      <- ex[1]

  data.table(q5=q5, q15=q15, q30=q30, e60=e60, e0=e0)
}

other.ltl <- list()
for (is in unique(population$iso3)){
  print(is)
  pop.is   <- population %>% filter(iso3 == is & year == 2019) %>%
    group_by(age) %>% summarise(Nx = sum(Nx)) %>% ungroup()
  dea.is <- ghe.all.causes %>% filter(iso3 == is & year == 2019)  %>%
    group_by(age) %>% summarise(val = sum(val)) %>% ungroup()
  mx.is <- left_join(pop.is, dea.is, by = "age") %>% arrange(age) %>%
    mutate(mx = val/Nx) %>% pull(mx)
  mx.is[mx.is > 1] <- 0.999
  mx.is[mx.is < 0|is.na(mx.is)] <- 0
  other.ltl[[is]] <- data.table(iso3 = is, est.lt.elem(mx.is))
}
other.lt <- rbindlist(other.ltl)

###############################################################################################
# Final raw covariates data base

covariates.raw1 <- who.covid %>%
  right_join(population %>% select(iso3) %>% unique(), by = "iso3") %>%
  right_join(other.lt, by = "iso3") %>%
  left_join(owid, by = c("iso3", "month", "year")) %>%
  left_join(covid_pol, by = c("iso3", "month", "year")) %>%
  left_join(temperatures, by = c("iso3", "month")) %>%
  left_join(pop.sum, by = "iso3") %>%
  left_join(sdi.gbd, by = "iso3") %>%
  left_join(gbd_in_diab , by = "iso3") %>%
  select(-c(life_expectancy)) %>%
  left_join(hdi , by = "iso3")

covariates.raw2 <- covariates.raw1 %>%
  mutate(Containment = ifelse(is.na(Containment) & year == 2020 & month < 4,
                              0, Containment),
         positive_rate = ifelse(is.na(positive_rate) & year == 2020 & month < 4,
                              0, positive_rate))

covariates_sexage <- covariates.raw2 %>% 
  group_by(WHO_region, month, year) %>%
  rep_with_median(., c("population_density", "life_expectancy", "positive_rate",
                       "Stringency", "Government", "Containment", "Economic",
                       "temperature", "ztemp",
                       "over65", "under15", "sdi", "diabetes_allages",
                       "hdi", "myeduc", "eyeduc", "gnipc")) %>%
  ungroup() %>%
  left_join(inc.groups, by = "iso3") %>%  full_join(gbd_in_oth, by = "iso3")

covariates <- covariates_sexage %>%
  filter(sex == "Both" & age == 999) %>% select(-c(sex, age)) %>%
  mutate(Containment = ifelse(is.na(Containment), 0, Containment),
         positive_rate = ifelse(is.na(positive_rate), 0, positive_rate))

#############################################################################################
# Check covariates

pacman::p_load(tidyverse, gghighlight)

covp <- covariates %>%
  mutate(months = month + 12*(year - 2020)) %>%
  select(Country, iso3, WHO_region, months, Containment, positive_rate)

ggplot(covp, aes(x = months, y = positive_rate, col = Country)) +
  geom_point() + geom_line() +
  facet_wrap(~WHO_region) +
  gghighlight(Country %in% c("Australia", "New Zealand", "Republic of Korea",
                             "United States of America", "Italy", "Egypt",
                             "South Africa", "India"),
              use_direct_label = FALSE, calculate_per_facet = TRUE)

ggplot(covp, aes(x = months, y = Containment, col = Country)) +
  geom_point() + geom_line() +
  facet_wrap(~WHO_region) +
  gghighlight(Country %in% c("Australia", "New Zealand", "Republic of Korea",
                             "United States of America", "Italy", "Egypt",
                             "South Africa", "India"),
              use_direct_label = FALSE, calculate_per_facet = TRUE)

########################################################################################################################
# Final datasets
#####################################################################################################################################

mf.df.in <- covariates %>%
  select(Country, iso3, WHO_region, year, income_group, high.income, upper.income,
         month, covid_cases, covid_deaths, temperature,ztemp,
         life_expectancy, q5, q15, q30, e60, e0,
         positive_rate, sdi, hdi, eyeduc, myeduc, gnipc,
         ncds, diabetes_allages, population_density,
         over65, under15, cardiovascular, hiv, Containment, Economic, Government, Stringency) %>%
  filter(year %in% c(2020, 2021)) %>% 
  left_join(wmd.df, by = c("iso3","year", "month")) %>%
  left_join(population %>% filter(year > 2019) %>% group_by(iso3, year) %>%
              summarise(Nx = sum(Nx), .groups = "drop"), by = c("iso3","year")) %>%
  left_join(ghe.all.causes %>% filter(year == 2019) %>% group_by(iso3) %>%
              summarise(ghe = sum(val), .groups = "drop"), by = c("iso3")) %>%
  left_join(fread("consulted.data.csv"), by = c("iso3", "year", "month")) %>% 
  mutate(observed = adj.ac) %>% 
  arrange(iso3, year, month)  

###########################################################################################################################
###########################################################################################################################
# Final dataset used in model

mf.df <- mf.df.in %>%
  mutate(months = month + 12*(year - 2020), 
         expectedr = 1e5*expected/Nx, 
         observedr = 1e5*observed/Nx,
         excess = observed - expected, excessr = 1e5*excess/Nx, 
         excessr = 1e5*excess/Nx,
         covidr = 1e5*covid_deaths/Nx, cmr = 1e5*covid_cases/Nx, expluscov = expectedr + covidr,
         wdate = format(date_decimal(year + (month-0.1)/12), "%Y-%m-%d"),
         wdate = as.Date(wdate), 
         positive_rate = positive_rate * 100, sdi = sdi*100) %>%
  rename(expectedl="expected.lwr", expectedu="expected.uppr") %>%
  select("Country", "iso3", "WHO_region", "income_group",
         "high.income", "upper.income", "Nx",
         "year", "month","wdate","temperature",
         "covid_cases","cmr", "covid_deaths", "covidr",
         "population_density","sdi","hdi", "eyeduc", "myeduc", "gnipc",
         "Stringency", "Government", "Containment", "Economic",
         "ncds","cardiovascular", "hiv", "over65", "under15","diabetes_allages", "positive_rate",
         "life_expectancy", q5, q15, q30, e60, e0, "observed", "observedr", 
         "expected","expectedl","expectedu", "expectedr", "excess", "excessr",
         "expluscov") %>%
  arrange(iso3, year, month) %>% mutate(rowid = 1:n()) %>%
  mutate(Country = ifelse(iso3 == "CIV", "Côte d'Ivoire", Country))

# Save input file
save(mf.df, file = "../Generated_Data/mf.df.Rda")

