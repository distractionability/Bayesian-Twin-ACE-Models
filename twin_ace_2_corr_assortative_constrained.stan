data {
  // Meta-data
  int<lower=0> mz_N;// Number of mz-pairs total
  int<lower=0> dz_N;
  int<lower=0> cohort_N; // Number of distinct birth years
  int<lower = 1> outcomes_N;//Number of observed outcome types (ordinal categorical - written for education)
  int<lower = 0> factor_levels_N[2];//Number of levels within the hierarchical variable
  int factor_levels[cohort_N, 2];
  int sub_factor_levels_N[factor_levels_N[1]]; // Number of levels within each level of the top factor

  //Individual observations
  int<lower = 1, upper = outcomes_N> y_mz_1[mz_N]; // Sorted by twin-pair id
  int<lower = 1, upper = outcomes_N> y_mz_2[mz_N]; 
  int<lower = 1, upper = cohort_N> cohort_mz[mz_N]; 
  
  int<lower = 1, upper = outcomes_N> y_dz_1[dz_N]; 
  int<lower = 1, upper = outcomes_N> y_dz_2[dz_N]; 
  int<lower = 1, upper = cohort_N> cohort_dz[dz_N]; 
}
transformed data {
  int cutpoints_N = outcomes_N - 1;
}
parameters {
  // For cutpoints
  simplex[outcomes_N] outcome_shares[cohort_N];
  
  real share_ac_mu;
  real share_a_of_ac_mu;
  
  // Level 2
  real<lower = 0> share_ac_sigma;
  real<lower = 0> share_a_of_ac_sigma;
  
  
  // Level 3
  real<lower = 0> share_ac_2_sigma;
  real<lower = 0> share_a_of_ac_2_sigma;
  

  // Shared within twin-pair (random effect - unscaled)
  vector[mz_N + dz_N] twinpair_std;
  vector[factor_levels_N[1]] share_ac_1_std;
  vector[factor_levels_N[1]] share_a_of_ac_1_std;
  vector[factor_levels_N[2]] share_ac_2_std;
  vector[factor_levels_N[2]] share_a_of_ac_2_std;

}
transformed parameters {
  matrix[cutpoints_N, cohort_N] cohort_cutpoints;
  vector[cohort_N] MZ_shared;
  vector[cohort_N] DZ_shared;
  vector[cohort_N] MZ_coef;
  vector[cohort_N] DZ_coef;
  // vector[cohort_N] mz_re;
  // vector[cohort_N] dz_re;
  vector[cohort_N] MZ_sigma; // Unexplained variance
  vector[cohort_N] DZ_sigma;
  vector[cohort_N] share_ac_re;
  vector[cohort_N] share_a_of_ac_re;
  vector[cohort_N] A_share; // Set equal to the assortative so as not to break later code with missing
  vector[cohort_N] C_share;// Set equal to the assortative so as not to break later code with missing
  vector[cohort_N] E_share;// Set equal to the assortative so as not to break later code with missing
  vector[cohort_N] A_share_assortative;
  vector[cohort_N] C_share_assortative;
  vector[cohort_N] E_share_assortative;

{
    vector[factor_levels_N[1]] factor_1_share_ac_re;
    vector[factor_levels_N[2]] factor_2_share_ac_re;
    vector[factor_levels_N[1]] factor_1_share_a_of_ac_re;
    vector[factor_levels_N[2]] factor_2_share_a_of_ac_re;
    int temp_passed = 0;

    factor_1_share_ac_re = share_ac_sigma * (share_ac_1_std - mean(share_ac_1_std));
    factor_1_share_a_of_ac_re = share_a_of_ac_sigma * (share_a_of_ac_1_std - mean(share_a_of_ac_1_std));

    //Separate sum-to-zero constraints for factor 2 within each level of factor 1
    for (level_i in 1:factor_levels_N[1]){
      factor_2_share_ac_re[(temp_passed + 1):(temp_passed + sub_factor_levels_N[level_i])] =
        share_ac_2_sigma * (share_ac_2_std[(temp_passed + 1):(temp_passed + sub_factor_levels_N[level_i])] -
                      mean(share_ac_2_std[(temp_passed + 1):(temp_passed + sub_factor_levels_N[level_i])]));
      factor_2_share_a_of_ac_re[(temp_passed + 1):(temp_passed + sub_factor_levels_N[level_i])] =
        share_a_of_ac_2_sigma * (share_a_of_ac_2_std[(temp_passed + 1):(temp_passed + sub_factor_levels_N[level_i])] -
                      mean(share_a_of_ac_2_std[(temp_passed + 1):(temp_passed + sub_factor_levels_N[level_i])]));
      temp_passed += sub_factor_levels_N[level_i];
    }
    for (cohort_i in 1:cohort_N){
        share_ac_re[cohort_i] = share_ac_mu  + factor_1_share_ac_re[factor_levels[cohort_i, 1]] +
            factor_2_share_ac_re[factor_levels[cohort_i, 2]];
        share_a_of_ac_re[cohort_i] = share_a_of_ac_mu  + factor_1_share_a_of_ac_re[factor_levels[cohort_i, 1]] +
            factor_2_share_a_of_ac_re[factor_levels[cohort_i, 2]];
    }

  }  
  

  for (cohort_index in 1:cohort_N){
  
    // With an assumption that DZ twins share 0.68 of their genes
    A_share_assortative[cohort_index] = inv_logit(share_a_of_ac_re[cohort_index]) * inv_logit(share_ac_re[cohort_index]);
    C_share_assortative[cohort_index] = (1 - inv_logit(share_a_of_ac_re[cohort_index])) * inv_logit(share_ac_re[cohort_index]);
    E_share_assortative[cohort_index] = 1 - inv_logit(share_ac_re[cohort_index]);
    // NOTE: THE STANDARD SHARES ARE SET EQUAL JUST TO AVOID MISSING - THEY WILL NOT BE USED
    A_share[cohort_index] =  A_share_assortative[cohort_index];
    C_share[cohort_index] =   C_share_assortative[cohort_index];
    E_share[cohort_index] =  E_share_assortative[cohort_index];
    
  
  }
  

  
  for (cohort_index in 1:cohort_N){
    cohort_cutpoints[, cohort_index] = inv_Phi(cumulative_sum(outcome_shares[cohort_index][1:cutpoints_N]));
    
    MZ_shared[cohort_index] = A_share_assortative[cohort_index] + C_share_assortative[cohort_index];
    DZ_shared[cohort_index] = 0.68 * A_share_assortative[cohort_index] + C_share_assortative[cohort_index];
    
  }
  

    MZ_coef = sqrt(MZ_shared);
    DZ_coef = sqrt(DZ_shared);
    MZ_sigma = sqrt(1 - MZ_shared);
    DZ_sigma = sqrt(1 - DZ_shared);

}
model {
  real temp_sum = 0;
  real upper_limit = 0;
  real lower_limit = 0;
  int twinpair_counter = 1;
  int twin_same = 0;
  
  // Priors

  for (cohort_index in 1:cohort_N){
    outcome_shares[cohort_index] ~ dirichlet(rep_vector(1,outcomes_N));
  }
  
  share_ac_mu ~ normal(0, 1.5); 
  share_a_of_ac_mu ~ normal(0, 1.5);
  share_ac_sigma ~ normal(0, 1);
  share_a_of_ac_sigma ~ normal(0, 1);
  share_ac_2_sigma ~ normal(0, 0.5);
  share_a_of_ac_2_sigma ~ normal(0, 0.5);
  
  
  share_ac_1_std ~ std_normal();
  share_a_of_ac_1_std ~ std_normal();
  share_ac_2_std ~ std_normal();
  share_a_of_ac_2_std ~ std_normal();
  twinpair_std ~ std_normal();
  
  // Model
  
  // MZ twins
  
  for (mz_index in 1:mz_N){
    twin_same = (y_mz_1[mz_index] == y_mz_2[mz_index])? 1 : 0; //If same outcome for both twins
    
    upper_limit = (y_mz_1[mz_index] == outcomes_N) ? 1 : 
                    exp(normal_lcdf(cohort_cutpoints[y_mz_1[mz_index], cohort_mz[mz_index]] | 
                                      MZ_coef[cohort_mz[mz_index]] * twinpair_std[mz_index],
                                      MZ_sigma[cohort_mz[mz_index]]));
    
    lower_limit = (y_mz_1[mz_index] == 1) ? 0 : 
                    exp(normal_lcdf(cohort_cutpoints[y_mz_1[mz_index] - 1, cohort_mz[mz_index]] | 
                                      MZ_coef[cohort_mz[mz_index]] * twinpair_std[mz_index],
                                      MZ_sigma[cohort_mz[mz_index]]));
    
    temp_sum += log(upper_limit - lower_limit) * ((twin_same == 1) ? 2 : 1);
    
    if (twin_same == 0){
      upper_limit = (y_mz_2[mz_index] == outcomes_N) ? 1 : 
                      exp(normal_lcdf(cohort_cutpoints[y_mz_2[mz_index], cohort_mz[mz_index]] | 
                                        MZ_coef[cohort_mz[mz_index]] * twinpair_std[mz_index],
                                        MZ_sigma[cohort_mz[mz_index]]));
      
      lower_limit = (y_mz_2[mz_index] == 1) ? 0 : 
                      exp(normal_lcdf(cohort_cutpoints[y_mz_2[mz_index] - 1, cohort_mz[mz_index]] | 
                                        MZ_coef[cohort_mz[mz_index]] * twinpair_std[mz_index],
                                        MZ_sigma[cohort_mz[mz_index]]));
      
      temp_sum += log(upper_limit - lower_limit);
    }
    
  }

  for (dz_index in 1:dz_N){
    twin_same = (y_dz_1[dz_index] == y_dz_2[dz_index])? 1 : 0; //If same outcome for both twins
    
    upper_limit = (y_dz_1[dz_index] == outcomes_N) ? 1 : 
                    exp(normal_lcdf(cohort_cutpoints[y_dz_1[dz_index], cohort_dz[dz_index]] | 
                                      DZ_coef[cohort_dz[dz_index]] * twinpair_std[dz_index + mz_N],
                                      DZ_sigma[cohort_dz[dz_index]]));
    
    lower_limit = (y_dz_1[dz_index] == 1) ? 0 : 
                    exp(normal_lcdf(cohort_cutpoints[y_dz_1[dz_index] - 1, cohort_dz[dz_index]] | 
                                      DZ_coef[cohort_dz[dz_index]] * twinpair_std[dz_index + mz_N],
                                      DZ_sigma[cohort_dz[dz_index]]));
    
    temp_sum += log(upper_limit - lower_limit) * ((twin_same == 1) ? 2 : 1);
    
    if (twin_same == 0){
      upper_limit = (y_dz_2[dz_index] == outcomes_N) ? 1 : 
                      exp(normal_lcdf(cohort_cutpoints[y_dz_2[dz_index], cohort_dz[dz_index]] | 
                                        DZ_coef[cohort_dz[dz_index]] * twinpair_std[dz_index + mz_N],
                                        DZ_sigma[cohort_dz[dz_index]]));
      
      lower_limit = (y_dz_2[dz_index] == 1) ? 0 : 
                      exp(normal_lcdf(cohort_cutpoints[y_dz_2[dz_index] - 1, cohort_dz[dz_index]] | 
                                        DZ_coef[cohort_dz[dz_index]] * twinpair_std[dz_index + mz_N],
                                        DZ_sigma[cohort_dz[dz_index]]));
      
      temp_sum += log(upper_limit - lower_limit);
    }
    
  }
  target += temp_sum;
  
  
}
