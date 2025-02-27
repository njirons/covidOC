data {
  int<lower=1> n_compartments; // number of compartments in epi model
  
  int n_d; // number of days of death observations/days to model    
  array[n_d] int d; // deaths over time
  array[n_d] int c; // cases over time
  int d0; // initial deaths
  int c0; // initial cases
  int n; // state populations
  array[n_d] int<lower=1> error_term_inds;
  int<lower=1> n_error_terms;
  
  real<lower=0, upper=1> init_props;
  
  real<lower=0> r0_upper;
  real<lower=0, upper=1> ifr_upper;
  real<lower=0> ccd_upper;
  
  real<lower=0> r0_lower;
  real<lower=0, upper=1> ifr_lower;
  real<lower=0> ccd_lower;
  
  real<lower=0> r0_mean;
  real<lower=0> r0_median;
  // real<lower=0, upper=1> ifr;
  real<lower=0, upper=1> ifr_mean;
  real<lower=1> ccd_mean;
  
  real<lower=0> r0_sd;
  real<lower=0> ifr_sd;
  real<lower=0> ccd_sd;
  
  // fixed epi parameters
  real<lower=0, upper=1> gamma;
  real<lower=0, upper=1> delta;
  real<lower=0, upper=1> mu;
  real<lower=0, upper=1> rr;
  
  // handle zero-inflation in data, if present
  int<lower=0, upper=n_d> n0_d;
  int<lower=0, upper=n_d> nz_d;
  array[nz_d] int<lower=1> d_nz_inds;
  int<lower=1> n0_d_inds;
  array[n0_d_inds] int d_z_inds;
  
  int<lower=0, upper=n_d> n0_c;
  int<lower=0, upper=n_d> nz_c;
  array[nz_c] int<lower=1> c_nz_inds;
  int<lower=1> n0_c_inds;
  array[n0_c_inds] int c_z_inds;
  
  real<lower=0> sigma_upper;
  real<lower=0> alpha_sigma;
  real<lower=0> beta_sigma;
  
  real<lower=0> kappa_d_upper;
}
parameters {
  vector<lower=0, upper=r0_upper>[n_error_terms] r0; // time-varying R0(t)
  // real<lower=r0_lower,upper=r0_upper> r0_max;
  // vector<lower=0, upper=1>[n_error_terms] r0_prop; // time-varying R0(t)
  real log_sigma;
  
  real<lower=ccd_lower> ccd; // case confirmation delay
  // real<lower=ccd_lower, upper=ccd_upper> ccd; // case confirmation delay
  // real log_c0;
  real<lower=0,upper=1> c0_prop;
  
  vector<lower=0, upper=1>[n_error_terms] car; // time-varying case ascertainment rate
  real log_car_sigma;
  
  simplex[n_compartments] x0; // initial SIR states   
  
  // real<lower=0,upper=kappa_d_upper> kappa_c;
  real log_kappa_c;
  
  real<lower=0, upper=1> d_hurdle; // death hurdle model zero probability
  
  real<lower=0, upper=1> c_hurdle; // case hurdle model zero probability
  real<lower=0, upper=1> c_disp; // case overdispersion mixture parameter
  
  // real<lower=ifr_lower,upper=ifr_upper> ifr; // infection fatality rate
  real<lower=0,upper=1> ifr;
}
transformed parameters {
  // real death_par0;
  // real case_par0;
  // real c_disp0;  
  
  vector[n_d] death_par; // death mean
  vector[n_d] death_par_obs; // death mean on observation days
  
  vector[n_d] case_par; // cases mean
  vector[n_d] case_par_obs; // cases mean on observation days
  vector[n_d] c_disp_par; // case dispersion
  vector[n_d] c_disp_obs; // case dispersion on observation days
  
  vector[n_d] sus; // proportion susceptible
  vector[n_d] exps; // proportion exposed  
  vector[n_d] inf; // proportion infected 
  vector[n_d] rem_s; // proportion removed, will survive
  vector[n_d] rem_d; // proportion removed, will die   
  
  vector[n_d] nu; // new infections
  vector[n_d] inf_c; // infections waiting to be confirmed
  
  // vector[n_error_terms] r0;
  // r0[1] = r0_max * r0_prop[1];
  // r0[1] = r0_max;
  
  sus[1] = init_props * x0[1] + (1 - init_props);
  exps[1] = init_props * x0[2];
  inf[1] = init_props * x0[3];
  rem_s[1] = init_props * (x0[4] + x0[n_compartments]) * (1-ifr);
  rem_d[1] = init_props * x0[4] * ifr;
  // rem_s[1] = 0;
  // rem_d[1] = 0;
  
  nu[1] = fmin(1,fmax(0,gamma * r0[1] * sus[1] * inf[1]));
  
  death_par[1] = n * fmax(0, rem_d[1] * mu);
  death_par_obs[1] = death_par[1];
  
  inf_c[1] = n * (1-c0_prop) * (1-sus[1]) * (1 - 1 / ccd) + n * car[1] * nu[1];
  case_par[1] =  n * (1-c0_prop) * (1-sus[1]) / ccd;  
  // inf_c[1] = exp(log_c0) * (1 - 1 / ccd) + n * car[1] * nu[1];
  // case_par[1] = exp(log_c0) / ccd;
  case_par_obs[1] = case_par[1];
  c_disp_par[1] = c_disp * case_par[1] * exp(-log_kappa_c)
                  + (1 - c_disp) * exp(-log_kappa_c);
  c_disp_obs[1] = c_disp_par[1]; 
  
  // death_par0 = n * init_props * x0[n_compartments] * ifr;
  // case_par0 = n * c0_prop * (1-sus[1]); 
  // c_disp0 = c_disp * case_par0 * exp(-log_kappa_c)
  //                 + (1 - c_disp) * exp(-log_kappa_c);  
  
  for (day in 1 : (n_d - 1)) {
    sus[day + 1] = fmin(1, fmax(0, sus[day] - nu[day] + rr * rem_s[day]));
    exps[day + 1] = fmin(1, fmax(0, exps[day] * (1 - delta) + nu[day]));
    inf[day + 1] = fmin(1,
                        fmax(0, inf[day] * (1 - gamma) + delta * exps[day]));
    rem_s[day + 1] = fmin(1,
                          fmax(0,
                               rem_s[day] * (1 - rr)
                               + inf[day] * gamma * (1 - ifr)));
    rem_d[day + 1] = fmin(1,
                          fmax(0,
                               rem_d[day] * (1 - mu) + inf[day] * gamma * ifr));
    
    // r0[error_term_inds[day + 1]] = r0_max * r0_prop[error_term_inds[day + 1]];
    nu[day + 1] = fmin(1,
                       fmax(0,gamma
                            * r0[error_term_inds[day + 1]]
                            * sus[day + 1] * inf[day + 1]));
    
    death_par[day + 1] = n * fmax(0, rem_d[day + 1] * mu);
    death_par_obs[day+1] = death_par[day+1];
    if(d[day] == 0){ // account for days with no reporting
        death_par_obs[day + 1] += d_hurdle * death_par_obs[day];      
    }
    
    inf_c[day + 1] = inf_c[day] * (1 - 1 / ccd)
                     + n * car[error_term_inds[day + 1]] * nu[day + 1];
    case_par[day + 1] = inf_c[day] / ccd;
    case_par_obs[day+1] = case_par[day+1];
    if(c[day] == 0){ // account for days with no reporting
        case_par_obs[day + 1] += c_hurdle * case_par_obs[day];      
    }    
    c_disp_par[day + 1] = c_disp * case_par[day + 1] * exp(-log_kappa_c)
                          + (1 - c_disp) * exp(-log_kappa_c);    
    c_disp_obs[day + 1] = c_disp * case_par_obs[day + 1] * exp(-log_kappa_c)
                          + (1 - c_disp) * exp(-log_kappa_c);
  }
}
model {
  // kappa_c ~ uniform(0,kappa_d_upper);
  // log_kappa_c ~ uniform(-10,kappa_d_upper);
  
  // ccd ~ normal(ccd_mean, ccd_sd) T[ccd_lower, ccd_upper];
  ccd ~ normal(ccd_mean, ccd_sd) T[ccd_lower,];  
  
  d_hurdle ~ uniform(0, 1);
  c_hurdle ~ uniform(0, 1);
  c_disp ~ uniform(0, 1);
  
  // d[d_nz_inds] ~ poisson(death_par_obs[d_nz_inds]);
  // c[c_nz_inds] ~ neg_binomial_2(case_par_obs[c_nz_inds], c_disp_obs[c_nz_inds]);
  if (n0_d > 0) {
    // ZIP model on deaths
    target += log_sum_exp(bernoulli_lpmf(1 | d_hurdle),
                          bernoulli_lpmf(0 | d_hurdle)
                          + poisson_lpmf(0 | death_par_obs[d_z_inds]));
    target += nz_d * bernoulli_lpmf(0 | d_hurdle);
    target += poisson_lpmf(d[d_nz_inds] | death_par_obs[d_nz_inds]);
  }
  if (n0_d == 0) {
    // Poisson model on deaths
    d ~ poisson(death_par_obs);
  }
  if (n0_c > 0) {
    // ZINB on cases
    target += log_sum_exp(bernoulli_lpmf(1 | c_hurdle),
                          bernoulli_lpmf(0 | c_hurdle)
                          + neg_binomial_2_lpmf(0 | case_par_obs[c_z_inds], c_disp_obs[c_z_inds]));
    target += nz_c * bernoulli_lpmf(0 | c_hurdle);
    target += neg_binomial_2_lpmf(c[c_nz_inds] | case_par_obs[c_nz_inds], c_disp_obs[c_nz_inds]);
  }
  if (n0_c == 0) {
    // NB model on cases
    c ~ neg_binomial_2(case_par_obs, c_disp_obs);
  }  
  
  // d0 ~ poisson(death_par0);
  // c0 ~ neg_binomial_2(case_par0,c_disp0);  
  d0 ~ poisson(n * init_props * (x0[4] + x0[n_compartments]) * ifr);
  // d0 ~ poisson(n * init_props * x0[n_compartments] * ifr);
               // + sum(death_par[1 : (min(d_nz_inds) - 1)]));
  // c0 ~ poisson(n * c0_prop * (1-sus[1]));            
  
  // independent error terms 
  r0[1] ~ uniform(0, r0_upper);
  // for (j in 2 : n_error_terms) {
  //   r0[j] ~ lognormal(r0[j-1],exp(log_sigma)) T[0,r0_upper];
  // }
  // r0[1] ~ normal(r0_median, r0_sd)T[r0_lower,r0_upper];
  for (j in 2 : n_error_terms) {
    target += beta_lpdf(r0[j]/r0_upper |
    exp(log_sigma)*r0[j-1]/r0_upper,
    exp(log_sigma)*(1-r0[j-1]/r0_upper));
  }
  // r0_max ~ normal(r0_median, r0_sd)T[r0_lower,r0_upper];
  // r0_prop[1] ~ uniform(0,1);
  // for (j in 2 : n_error_terms) {
  //   r0_prop[j] ~ beta(r0_prop[j-1]*exp(log_sigma),(1-r0_prop[j-1])*exp(log_sigma));
  // }

  car[1] ~ uniform(0, 1);
  car[2 : n_error_terms] ~ beta(exp(log_car_sigma)
                                * car[1 : (n_error_terms - 1)],
                                exp(log_car_sigma)
                                * (1 - car[1 : (n_error_terms - 1)]));
  
  // uniform prior on initial SEIRD compartment values
  x0 ~ dirichlet(rep_vector(1, n_compartments));
  
  // prior on IFR from Irons and Raftery (2021)
  // ifr ~ normal(ifr_mean,ifr_sd)T[ifr_lower,ifr_upper];
  ifr ~ normal(ifr_mean,ifr_sd)T[0,1];
}
generated quantities {
  // array[n_d] int dhat;
  // array[n_d] int chat;
  // array[n_d] int dhat_obs;
  // array[n_d] int chat_obs;
  // 
  // dhat = poisson_rng(death_par);
  // chat = neg_binomial_2_rng(case_par, c_disp_par);
  // dhat_obs = poisson_rng(death_par_obs);
  // chat_obs = neg_binomial_2_rng(case_par_obs, c_disp_obs); 
  
  array[n_d] int dhat_nz;
  array[n_d] int chat_nz;
  array[n_d] int dhat_nz_obs;
  array[n_d] int chat_nz_obs;
  
  dhat_nz = poisson_rng(death_par);
  chat_nz = neg_binomial_2_rng(case_par, c_disp_par);
  dhat_nz_obs = poisson_rng(death_par_obs);
  chat_nz_obs = neg_binomial_2_rng(case_par_obs, c_disp_obs);   
}


