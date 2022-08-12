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
  int<lower=0> N;
  int Y1 [Tt];
  int Y [Tt];
  int Y1T [N];
  real<lower=0> pc_U; // 
  real<lower=0, upper=1> pc_alpha; //
}

parameters {
  real mu;
  vector[Tt] eps_raw;
  real<lower = 0> eps_tau;
}

transformed parameters{
  vector[Tt] logp;
  vector[Tt] eps;
  real<lower = 0> eps_sd;
  
  eps_sd = inv_sqrt(eps_tau);
  eps = eps_sd * eps_raw;
  
  for(i in 1:Tt){
    logp[i]  = mu + eps[i];
  }
}  

model {
  //p ~  beta(1,1);
  eps_tau ~ pc_normal_precision(pc_U, pc_alpha);
  eps_raw ~ std_normal();
  
  for(i in 1:Tt){
    Y1[i] ~  binomial(Y[i],(1/(1+exp(-logp[i]))));
  }
  
}

generated  quantities {
  
  int Y2T [N];
  real<lower=0,upper=1> p2T [N];
  
  for(k in 1:N){
  p2T[k] = 1/(1+exp(-(mu + normal_rng(0,eps_sd))));
  Y2T[k] = neg_binomial_rng(Y1T[k],(p2T[k]/(1-p2T[k])));
  }
}
