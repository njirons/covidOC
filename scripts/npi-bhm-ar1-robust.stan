functions {
  /* compute correlated group-level effects
   * Args:
   *   z: matrix of unscaled group-level effects
   *   SD: vector of standard deviation parameters
   *   L: cholesky factor correlation matrix
   * Returns:
   *   matrix of scaled group-level effects
   */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  
  array[] vector tMVN(vector mu, matrix L, vector b, vector s, vector u) {
    int K = rows(mu);
    vector[K] d;
    vector[K] z;
    array[2] vector[K] out;
    for (k in 1 : K) {
      int km1 = k - 1;
      if (s[k] != 0) {
        real z_star = (b[k]
                       - (mu[k]
                          + ((k > 1) ? L[k, 1 : km1] * head(z, km1) : 0)))
                      / L[k, k];
        real v;
        real u_star = Phi(z_star);
        if (s[k] == -1) {
          v = u_star * u[k];
          d[k] = u_star;
        } else {
          d[k] = 1 - u_star;
          v = u_star + d[k] * u[k];
        }
        z[k] = inv_Phi(v);
      } else {
        z[k] = inv_Phi(u[k]);
        d[k] = 1;
      }
    }
    out[1] = z;
    out[2] = d;
    return out;
  }
  
  array[] real normal_ub_rng(vector mu, real sigma, real ub) {
    int n = rows(mu);
    array[n] real y;
    for (j in 1 : n) {
      real p_ub = normal_cdf(ub | mu[j], sigma);
      real u = uniform_rng(0, p_ub);
      y[j] = mu[j] + sigma * inv_Phi(u);
    }
    return y;
  }
}
data {
  int<lower=1> p; // number of NPIs
  int<lower=1> nX; // number of non-NPI covariates (deaths, infections, etc.)
  int<lower=p + nX + 1> K; // total number of regressors (intercept plus NPIs plus other covariates)
  vector[K] b; // bounds
  // number of dimensions
  // lower or upper bound
  // b_type[k] ==  0 implies no constraint; otherwise
  // b_type[k] == -1 -> b[k] is an upper bound
  // b_type[k] == +1 -> b[k] is a lower bound
  vector<lower=-1, upper=1>[K] b_type;
  
  int<lower=1> s; // number of states
  array[s + 1] int t; // days of observations
  matrix[t[s + 1], p] u; // NPIs
  row_vector[p] u_mean; // NPI mean
  matrix[t[s + 1], nX] X; // non-NPI covariates
  
  vector[t[s + 1]] log_rt; // log(R(t))
  
  real log_r0_upper;
  real r0_upper;
  real r0_lower;
  real r0_median;
  real r0_mean;
  real r0_sd;
}
parameters {
  array[s] vector<lower=0, upper=1>[K] u_s;
  
  real<lower=-1, upper=1> phi_ar; // ar(1)
  real<lower=0> sigma; // dispersion parameter
  real log_nu; // degrees of freedom
  
  vector<lower=0>[K] sd_re; // group-level standard deviations
  cholesky_factor_corr[K] L_corr; // cholesky factor of random effect correlation matrix
  vector[K] mu; // global effects
}
transformed parameters {
  cholesky_factor_cov[K] L_cov = diag_pre_multiply(sd_re, L_corr);
  
  real R0 = exp(mu[1] - dot_product(mu[2 : (p + 1)], u_mean));
  vector[s] R0_s;
  
  array[s] vector[K] mu_s; // state-specific random effects
  
  vector[t[s + 1]] lin_pred; //linear predictor
  vector[t[s + 1]] reg_mean; //regression mean
  
  real lprior = 0; // prior contributions to the log posterior
  vector[t[s + 1]] log_lik; // pointwise log likelihood   
  
  for (state in 1 : s) {
    array[2] vector[K] tMVNs = tMVN(mu, L_cov, b, b_type, u_s[state]);
    mu_s[state] = mu + L_cov * tMVNs[1];
    lprior += sum(log(tMVNs[2])); // Jacobian adjustments
    
    R0_s[state] = exp(mu_s[state, 1]
                      - dot_product(mu_s[state, 2 : (p + 1)], u_mean));
    
    for (day in 1 : (t[state + 1] - t[state])) {
      lin_pred[day + t[state]] = mu_s[state, 1]
                                 + dot_product(mu_s[state, 2 : (p + 1)],
                                               u[day + t[state],  : ]
                                               - u_mean)
                                 + dot_product(mu_s[state, (p + 2) : K],
                                               X[day + t[state],  : ]);
      
      reg_mean[day + t[state]] = lin_pred[day + t[state]];
      if (day > 1) {
        reg_mean[day + t[state]] += phi_ar
                                    * (log_rt[day + t[state] - 1]
                                       - lin_pred[day + t[state] - 1]);
      }
      
      log_lik[day + t[state]] = student_t_lpdf(log_rt[day + t[state]] | exp(log_nu), reg_mean[day
                                                                    + t[state]], sigma);
    }
  }
  
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_re | 3, 0, 2.5)
            - K * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_corr | 1);
}
model {
  // PRIOR
  target += lprior;
  
  // LIKELIHOOD
  target += sum(log_lik);
}
generated quantities {
  array[t[s + 1]] real log_rt_hat = student_t_rng(exp(log_nu), reg_mean,
                                                  sigma);
  vector[t[s + 1]] eps = log_rt - reg_mean;
}


