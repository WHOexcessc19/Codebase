### load in packages
library(tidyverse)
library(INLA) # make sure you have installed the latest testing branch of INLA package
library(posterior)

### load in general country and covariate data
load("Generated_Data/mf.df.Rda")
### load in expecteds estimates
load("Generated_Data/acm_monthly_predictions_tier1.RData")
load("Generated_Data/acm_monthly_predictions_tier2.RData")
### load in annual mortality data
AnnualMortality <- read_csv("Imported_Data/AnnualMortality.csv") %>%
  filter(country == "China")

### mutate country, covariate dataframe for covariate model, and standardize relevant variables
df.inla <- mf.df %>% 
  mutate(months = ifelse(year == 2020, month, month + 12)) %>%
  mutate(months.pr = months,
         months.sq_covidr = months,
         months.contain = months,
         months.temp = months,
         observed = round(observed),
         expected = round(expected),
         covidr = ifelse(covidr<0,0,covidr),
         sqrt_covidr = sqrt(covidr),
         months.income_pr = months,
         months.income_covidr = months,
         pr_high.income = high.income*positive_rate,
         covidr_high.income  = high.income*sqrt_covidr,
         months.income_contain = months,
         AMRO.reg = ifelse(WHO_region=="AMRO",1,0),
         AFRO.reg = ifelse(WHO_region=="AFRO",1,0),
         EURO.reg = ifelse(WHO_region=="EURO",1,0),
         EMRO.reg = ifelse(WHO_region=="EMRO",1,0),
         WPRO.reg = ifelse(WHO_region=="WPRO",1,0),
         SEARO.reg = ifelse(WHO_region=="SEARO",1,0),
         contain_high.income = Containment*high.income,
         country.time = as.factor(paste(Country,months,sep = ".")),
         country.pred = Country) %>%
  mutate(observed = ifelse(Country=="India",NA,observed)) %>%
  mutate(observed = ifelse(Country=="Suriname"&months==24,NA,observed)) %>%
  mutate_at(c("sqrt_covidr","sdi","positive_rate","Containment","temperature","population_density","ncds","cardiovascular","hiv","over65","diabetes_allages","life_expectancy"), ~(scale(.) %>% as.vector))

### bind together expecteds dfs and merge into covariate model df
new.expected <- rbind(acm_monthly_predictions_tier1,acm_monthly_predictions_tier2)
df.inla <- df.inla %>%
  left_join(new.expected,by=c("iso3","year","month")) %>%
  dplyr::select(-expected) %>%
  rename(expected = expected_acm)

### Panamas ACM in April 2020 was counted in May 2020
### So adjust Panama months 3,4 to have same sum observed ACM, but slope from surrounding points
Panama.ind = which(df.inla$Country=="Panama")
Panama.slope = (df.inla$observed[Panama.ind[6]]-df.inla$observed[Panama.ind[3]])/(6-3)
new.Panama.month3and4 = c(sum(df.inla$observed[Panama.ind[4:5]])/2-Panama.slope/2,sum(df.inla$observed[Panama.ind[4:5]])/2+Panama.slope/2)
df.inla$observed[Panama.ind[4:5]]= new.Panama.month3and4

### Set Up INLA Covariate Model

### Model parameters
pc.u <- 1
pc.alpha <- 0.01
hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u, pc.alpha)))
control.family1 = list(control.link = list(model = "log"))

### Model formula specification with RW2 and sum to 0 constraint on time-varying covariates
model.formula = formula(observed ~ high.income  + positive_rate 
                        + sqrt_covidr + pr_high.income 
                        + covidr_high.income + Containment
                        + contain_high.income + temperature
                        + cardiovascular + diabetes_allages
                        + f(months.pr, positive_rate, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(months.sq_covidr, sqrt_covidr, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(months.income_pr, pr_high.income, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(months.income_covidr, covidr_high.income, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(months.contain, Containment, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(months.income_contain, contain_high.income, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(months.temp,temperature, model = "rw2", constr = TRUE, hyper = hyperpc1)
                        + f(country.time, model = "iid", hyper = hyperpc1))

### Run INLA Model
### Uses expecteds as offset, with gamma_delta expecteds uncertainty using "variant = 2" from INLA testing branch
pois.pred.INLA <- INLA::inla(model.formula,
                             data = df.inla, offset = log(gamma_E), scale = gamma_delta, family = "nbinomial",
                             control.predictor= list(compute = TRUE, link = 1),
                             control.compute = list(config = TRUE, cpo = TRUE),
                             control.family = list(variant = 2,
                                                   control.link = list(model = "log"),
                                                   hyper =  list(theta = list(initial = 0, fixed = TRUE))))

### Number of samples of estimates
num_inla_samps = 1000

### sample from posterior of INLA model for theta's
sample.df <- INLA::inla.posterior.sample(num_inla_samps,pois.pred.INLA)

ind = which(df.inla$Country=="China")
p = matrix(NA,nrow = 24,ncol = 1000)
for(k in 1:1000){
  ### calculate expecteds
  E = rgamma(24,shape=df.inla$gamma_delta[ind],rate=df.inla$gamma_delta[ind]/df.inla$gamma_E[ind])
  theta = sample.df[[k]]$latent[ind]
  # proportion of annual deaths in each month, defined by covariate model
  p[,k] = exp((E*theta)/sum(E*theta))
}
China.est = rbind(p[1:12,]*AnnualMortality$deaths[which(AnnualMortality$country=="China"&AnnualMortality$year==2020)],
                  p[13:24,]*AnnualMortality$deaths[which(AnnualMortality$country=="China"&AnnualMortality$year==2021)])

save(China.est,file = "Generated_Data/China.est.RData")
