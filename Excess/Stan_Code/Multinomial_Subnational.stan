functions {
    /*
    * The log density for the penalized complexity prior for the precision of a 
    * normal random effect.
    * U and alpha are chosen such that P(1 / sqrt(tau) > U) = alpha
    *
    * @param U 
    * @param alpha 
    * @return The log density of tau for the given U, alpha
    */
    real pc_normal_precision_lpdf(real tau, real U, real alpha) {
        real par = -log(alpha) / U;
        return log(par) - (1.5) * log(tau) - par * (tau) ^ (-0.5);
    }
    
}

data {
  int<lower=0> Tt; 
  int<lower=0> K; 
  vector[Tt] Ytplus;
  vector[Tt] Ytplusobs;
  matrix[Tt,K] Yt;
  matrix[Tt,K] Ytobs;
  int<lower=0> Tt_covid; 
  vector[Tt_covid] Ytplusobs_covid;
  matrix[Tt_covid,K] Yt_covid;
  matrix[Tt_covid,K] Ytobs_covid;
  real<lower=0> pc_U; // 
  real<lower=0, upper=1> pc_alpha; //
}

parameters {
  real alpha[K];
  real epsilon[Tt];
  real<lower=0> tau_eps;
}

transformed parameters{
  matrix[Tt,K + 1] logitp;
  vector[Tt] norm;
  vector[Tt] p_unobs;
  matrix[Tt,K+1] p;
  
  for(i in 1:Tt){
    for(k in 1:K){
      
      logitp[i,k] = alpha[k] + epsilon[i];
      
    }
    logitp[i,K+1] = 0;
    norm[i] = sum(exp(logitp[i, ]));
    p[i, ] = exp(logitp[i, ]) / norm[i];
    
    p_unobs[i] = p[i, K+1];
    for(k in 1:K){
      if(Ytobs[i, k] == 0){
        p_unobs[i] += p[i, k];
      }
    }
  }
  
}


model {
  
  // priors for parametters
  alpha ~ normal(0, 10);
  tau_eps ~ pc_normal_precision(pc_U, pc_alpha);
  epsilon ~ normal(0, inv_sqrt(tau_eps));
  
  for(i in 1:Tt){
    for(k in 1:K){
      if(Ytobs[i, k] == 1){
        target += Yt[i, k] * log(p[i, k]);
      }
    }
    target += (Ytplus[i] - Ytplusobs[i]) * log(p_unobs[i]);
  }
  
}

generated quantities{
  real epsilon_covid[Tt_covid];
  matrix[Tt_covid, K + 1] logitp_covid;
  vector[Tt_covid] norm_covid;
  vector[Tt_covid] p_unobs_covid;
  matrix[Tt_covid,K+1] p_covid;
  real Ytplus_covid[Tt_covid];
  
  for(i in 1:Tt_covid){
    epsilon_covid[i] = normal_rng(0, inv_sqrt(tau_eps));
  }
  
  for(i in 1:Tt_covid){
    for(k in 1:K){
      
      logitp_covid[i,k] = alpha[k]+epsilon_covid[i];
      
    }
    logitp_covid[i,K+1] = 0;
    norm_covid[i] = sum(exp(logitp_covid[i, ]));
    p_covid[i, ] = exp(logitp_covid[i, ]) / norm_covid[i];
    
    p_unobs_covid[i] = p_covid[i, K+1];
    for(k in 1:K){
      if(Ytobs_covid[i, k] == 0){
        p_unobs_covid[i] += p_covid[i, k];
      }
    }
  }
  
  
  for(i in 1:Tt_covid){
    Ytplus_covid[i] = Ytplusobs_covid[i] + 
    neg_binomial_rng(Ytplusobs_covid[i],  (1 - p_unobs_covid[i]) / (p_unobs_covid[i]));
  }

}








