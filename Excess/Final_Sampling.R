
# Combined Final Sampling Code

### load in packages
library(tidyverse)
library(INLA) # make sure you have installed the latest testing branch of INLA package
library(posterior)

#### Set the seed for this script ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### load in general country and covariate data
load("../Generated_Data/mf.df.Rda")
### load in expecteds estimates
load("../Generated_Data/acm_monthly_predictions_tier1.RData")
load("../Generated_Data/acm_monthly_predictions_tier2.RData")
### load in subnational and mixed data estimates
load("../Generated_Data/Argentina.est.RData")
load("../Generated_Data/India.est.RData")
load("../Generated_Data/Indonesia.est.RData")
load("../Generated_Data/Turkey.est.RData")
load("../Generated_Data/China.est.RData")
load("../Generated_Data/AnnualData.Ests.RData")

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
df.inla$observed[Panama.ind[4:5]]= as.integer(new.Panama.month3and4)



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
### How many country time points
country.time.n = nrow(df.inla)
### Define models and model names to consider
models <- c("Final Model")
pred.in.pois <- pois.pred.INLA

### Function for Getting ACM / Expecteds  / Excess Estimates

INLA.estimates.sampling <- function(poisson.pred,mod,
                                        country.time.n,num_inla_samps){
  
  final.df <- vector(mode = "list", length = 3)
  
  nb.draw <- function(expected.delta,theta.i){
    rnbinom(1, size = expected.delta, mu = exp(theta.i))
  }
  
  ### find indices of transition country-time points from observed to modeled
  ### this is for the purpose of benchmarking
  summarise.df <- df.inla %>%
    mutate(observe.na.ind = ifelse(is.na(observed),1,0)) %>%
    group_by(Country) %>%
    summarise(sum.na = sum(observe.na.ind)) %>%
    filter(sum.na>0 & sum.na<24)
  trans.first.df <- df.inla %>%
    mutate(ind = 1:nrow(df.inla)) %>%
    filter(Country %in% summarise.df$Country) %>%
    mutate(observe.na.ind = ifelse(is.na(observed),1,0)) %>%
    group_by(Country) %>%
    filter(observe.na.ind==1) %>%
    slice(1) %>%
    dplyr::select(Country,months,observe.na.ind,ind)
  partial.country.ind <- trans.first.df$ind
  
    
    # sample of num_inla_samps from posterior
    sample.df <- INLA::inla.posterior.sample(num_inla_samps,poisson.pred,seed = 25)
    
    # create dfs of Country, WHO Region, and months to merge in estimates with later
    excess.df.expec <- df.inla %>%
      dplyr::select(Country,WHO_region,months)
    excess.df.acm <- df.inla %>%
      dplyr::select(Country,WHO_region,months)
    excess.df.excess <- df.inla %>%
      dplyr::select(Country,WHO_region,months)
    
    # loop over samples to replace estimates from subnational and mixed models, and fix benchmarking
    for(i in 1:num_inla_samps){
      
      print(paste("Starting Iteration: ",i))
      
      # extracting posterior estimates of parameters
      theta.i <- sample.df[[i]]$latent[1:country.time.n]
      
      # for the time points for countries with partial data
      set.seed(42 + i)
      
      # benchmarking for countries with partial data
      for(m in 1:length(partial.country.ind)){
        
        country.ind <- partial.country.ind[m]
        prev.ind <- country.ind-1
        
        month.m <- df.inla$months[country.ind]
        change.ind <- country.ind:(country.ind+(24-month.m))
        
        theta.i[change.ind] <- theta.i[change.ind] + log(df.inla$observed[prev.ind]) - theta.i[prev.ind]
      }
      
      # benchmarking for countries with annual data estimates
      for(m in 1:nrow(AnnualData.Ests.df)){
        
        country <- AnnualData.Ests.df$Country[m]
        month.m <- AnnualData.Ests.df$months[m]
        
        if(month.m==12){
          
          country.ind = which(df.inla$Country==country & df.inla$months==13)
          prev.ind = country.ind-1
          
          change.ind <- country.ind:(country.ind+(24-month.m))
          
          theta.i[change.ind] <- theta.i[change.ind] + log(unname(rowMeans(AnnualData.Ests.df[,3:(num_inla_samps+2)]))[m]) - theta.i[prev.ind]
        }
      }
      
      # Create iteration i covariate model estimates df
      set.seed(42 + i)
      excess.df.i <- data.frame(theta.i=theta.i,expected.delta=df.inla$gamma_delta,
                                expected.E=df.inla$gamma_E,Country=df.inla$Country,
                                months=df.inla$months,
                                observed.ind=ifelse(!is.na(df.inla$observed),TRUE,FALSE),
                                observed=df.inla$observed,
                                expected.sim=rgamma(country.time.n,shape=df.inla$gamma_delta,rate=df.inla$gamma_delta/df.inla$gamma_E))
      
      # iteration i vectors for expecteds, acm, and excess
      # estimated acm are replaced by observed when available
      expec.i <- excess.df.i$expected.sim
      acm.i <- ifelse(excess.df.i$observed.ind,
                      df.inla$observed,
                      rnbinom(country.time.n, size = excess.df.i$expected.delta, mu = exp(theta.i)))
      excess.i <- acm.i- excess.df.i$expected.sim 
      
      # loop through and replace relevant acm/excess vector elements with subnational and mixed model estimates
      for(h in 1:country.time.n){
        
        if(excess.df.i$Country[h]=="Argentina" & excess.df.i$months[h]>12){
          month.i <- excess.df.i$months[h]-12
          Arg.est <- Argentina.est[month.i,i]
          acm.i[h] <- Arg.est
          excess.i[h] <- Arg.est - rgamma(1,shape = excess.df.i$expected.delta[h],rate = excess.df.i$expected.delta[h]/excess.df.i$expected.E[h])
        }
        
        if(excess.df.i$Country[h]=="Turkey"){
          month.i <- excess.df.i$months[h]
          Turk.est <- Turkey.est[month.i,i]
          acm.i[h] <- Turk.est
          excess.i[h] <- Turk.est - rgamma(1,shape = excess.df.i$expected.delta[h],rate = excess.df.i$expected.delta[h]/excess.df.i$expected.E[h])
        }
        
        if(excess.df.i$Country[h]=="India"){
          month.i <- excess.df.i$months[h]
          Ind.est <- India.est[month.i,i]
          acm.i[h] <- Ind.est
          excess.i[h] <- Ind.est - rgamma(1,shape = excess.df.i$expected.delta[h],rate = excess.df.i$expected.delta[h]/excess.df.i$expected.E[h])
        }
        
        if(excess.df.i$Country[h]=="China"){
          month.i <- excess.df.i$months[h]
          Ch.est <- China.est[month.i,i]
          acm.i[h] <- Ch.est
          excess.i[h] <- Ch.est - rgamma(1,shape = excess.df.i$expected.delta[h],rate = excess.df.i$expected.delta[h]/excess.df.i$expected.E[h])
        }
        
        if((excess.df.i$Country[h] %in% unique(AnnualData.Ests.df$Country)) & excess.df.i$months[h]<13){
          month.i <- excess.df.i$months[h]
          country.i <- excess.df.i$Country[h]
          AnnualData.est <- AnnualData.Ests.df %>%
            filter(Country == country.i) %>%
            filter(months == month.i) %>%
            dplyr::select(-c(Country,months))
          AnnualData.est <- AnnualData.est[,i]
          acm.i[h] <- AnnualData.est
          excess.i[h] <- AnnualData.est - rgamma(1,shape = excess.df.i$expected.delta[h],rate = excess.df.i$expected.delta[h]/excess.df.i$expected.E[h])
        }
        
        if(excess.df.i$Country[h]=="Indonesia" & excess.df.i$months[h]<19){
          month.i <- excess.df.i$months[h]
          Ind.est <- Indonesia.est[month.i,i]
          acm.i[h] <- Ind.est
          excess.i[h] <- Ind.est - rgamma(1,shape = excess.df.i$expected.delta[h],rate = excess.df.i$expected.delta[h]/excess.df.i$expected.E[h])
        }
      }
      
      # Add Sample i excess for each country, time point to df
      excess.df.expec <- cbind(excess.df.expec,as.data.frame(expec.i))
      excess.df.acm <- cbind(excess.df.acm,as.data.frame(unlist(acm.i)))
      excess.df.excess <- cbind(excess.df.excess,as.data.frame(unlist(excess.i)))
    }
    
  final.df <- list(excess.df.acm, excess.df.expec, excess.df.excess)
  return(final.df)
  
} #INLA.estimates.sampling function



## Setup and Run the Country Estimates Function

country.df <- INLA.estimates.sampling(pred.in.pois,models,country.time.n,
                                                    num_inla_samps)

### Format Estimates

expected = country.df[[2]]
acm = country.df[[1]]
excess = country.df[[3]]
colnames(expected) <- c("Country", "WHO_region", "months", 
                        paste0("expec", 1:num_inla_samps))
colnames(acm) <- c("Country", "WHO_region", "months", 
                   paste0("acm", 1:num_inla_samps))
colnames(excess) <- c("Country", "WHO_region", "months", 
                      paste0("excess", 1:num_inla_samps))

excess_summarized = cbind(excess[,1:3],est=rowMeans(excess[,4:(num_inla_samps+3)]),
                          lwr=apply(excess[,4:(num_inla_samps+3)],1,quantile,0.025),
                          uppr=apply(excess[,4:(num_inla_samps+3)],1,quantile,0.975))
expected_summarized = cbind(expected[,1:3],est=rowMeans(expected[,4:(num_inla_samps+3)]),
                            lwr=apply(expected[,4:(num_inla_samps+3)],1,quantile,0.025),
                            uppr=apply(expected[,4:(num_inla_samps+3)],1,quantile,0.975))
acm_summarized = cbind(acm[,1:3],est=rowMeans(acm[,4:(num_inla_samps+3)]),
                       lwr=apply(acm[,4:(num_inla_samps+3)],1,quantile,0.025),
                       uppr=apply(acm[,4:(num_inla_samps+3)],1,quantile,0.975))

### Write Estimates Full Sample DFs and Summarized Estimate DFs to csv
write.csv(excess, file="Generated_Data/excess.csv", row.names=FALSE)
write.csv(expected, file="Generated_Data/expected.csv", row.names=FALSE)
write.csv(acm, file="Generated_Data/acm.csv", row.names=FALSE)

write.csv(excess_summarized, file="Generated_Data/excess_summarized.csv", row.names=FALSE)
write.csv(expected_summarized, file="Generated_Data/expected_summarized.csv", row.names=FALSE)
write.csv(acm_summarized, file="Generated_Data/acm_summarized.csv", row.names=FALSE)

