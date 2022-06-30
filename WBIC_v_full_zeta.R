t <- DS_C$TObs
x <- DS_C$EntAge
bL <- DS_C$BL
bH <- DS_C$BH
c0 <- DS_C$C0
c1 <- DS_C$C1
c2 <- DS_C$C2
d <- DS_C$Status

t1 <- DS1$TObs
x1 <- DS1$EntAge
bL <- DS1$BL
bH <- DS1$BH
d1 <- DS1$Status

t2 <- DS2$TObs
x2 <- DS2$EntAge
c0 <- DS2$C0
c1 <- DS2$C1
c2 <- DS2$C2
d2 <- DS2$Status


######################################################## - J=0  - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00004, 0.1)

init_tau2 <- c(0.04, 0.06)

# - Set prior distributions
phi_0_prior_shape <- 0.09

phi_scale <- 1


# - Parameter storage

par_postJ0_gibbs_norm1 <- matrix(NA, 100000, 2)
#par_postJ0_gibbs_norm2 <- matrix(NA, 100000, 2)

#par_postJ0_metr_norm1 <- matrix(NA, 100000, 2)
#par_postJ0_metr_norm2 <- matrix(NA, 100000, 2)

par_postJ0_gibbs_norm1[1,] <- init_tau1
#par_postJ0_gibbs_norm2[1,] <- init_tau2

#par_postJ0_metr_norm1[1,] <- init_tau1
#par_postJ0_metr_norm2[1,] <- init_tau2

colnames(par_postJ0_gibbs_norm1) <- c("exp_alpha", "beta")
#colnames(par_postJ0_gibbs_norm2) <- c("exp_alpha", "beta")

#colnames(par_postJ0_metr_norm1) <- c("exp_alpha", "beta")
#colnames(par_postJ0_metr_norm2) <- c("exp_alpha", "beta")


##### - MCMC with Gibbs steps

for(i in 2:100000){
  
  phi0_1 <- par_postJ0_gibbs_norm1[i-1,1]
  beta_1 <- par_postJ0_gibbs_norm1[i-1,2]
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t) - 1) / beta_1) * exp(beta_1 * (x - 77.5)))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d / log(nrow(DS_C))), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ0_gibbs_norm1[i, 1] <- phi0
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t) - 1) / beta_star) * phi0 * exp(beta_star * (x - 77.5)) + d * (log(phi0) + beta_star * (x + t - 77.5))) / log(nrow(DS_C)) - 
           sum(- ((exp(beta_1    * t) - 1) / beta_1   ) * phi0 * exp(beta_1    * (x - 77.5)) + d * (log(phi0) + beta_1    * (x + t - 77.5))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ0_gibbs_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ0_gibbs_norm2[i-1,1]
  beta_1 <- par_postJ0_gibbs_norm2[i-1,2]
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t) - 1) / beta_1) * exp(beta_1 * (x - 77.5)))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d / log(nrow(DS_C))), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ0_gibbs_norm2[i, 1] <- phi0
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t) - 1) / beta_star) * phi0 * exp(beta_star * (x - 77.5)) + d * (log(phi0) + beta_star * (x + t - 77.5))) / log(nrow(DS_C)) - 
           sum(- ((exp(beta_1    * t) - 1) / beta_1   ) * phi0 * exp(beta_1    * (x - 77.5)) + d * (log(phi0) + beta_1    * (x + t - 77.5))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ0_gibbs_norm2[i, 2] <- beta
  
}


##### - MCMC with Metropolis steps

#for(i in 2:100000){
  
  phi0_1 <- par_postJ0_metr_norm1[i-1,1]
  beta_1 <- par_postJ0_metr_norm1[i-1,2]
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t) - 1) / beta_1) * phi0_star * exp(beta_1 * (x - 77.5)) + d * (log(phi0_star) + beta_1 * (x + t - 77.5))) / log(nrow(DS_C)) - 
           sum(- ((exp(beta_1 * t) - 1) / beta_1) * phi0_1    * exp(beta_1 * (x - 77.5)) + d * (log(phi0_1   ) + beta_1 * (x + t - 77.5))) / log(nrow(DS_C))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.008) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.008))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ0_metr_norm1[i, 1] <- phi0
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t) - 1) / beta_star) * phi0 * exp(beta_star * (x - 77.5)) + d * (log(phi0) + beta_star * (x + t - 77.5))) / log(nrow(DS_C)) - 
           sum(- ((exp(beta_1    * t) - 1) / beta_1   ) * phi0 * exp(beta_1    * (x - 77.5)) + d * (log(phi0) + beta_1    * (x + t - 77.5))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ0_metr_norm1[i, 2] <- beta
  
}

#

#for(i in 2:100000){
  
  phi0_1 <- par_postJ0_metr_norm2[i-1,1]
  beta_1 <- par_postJ0_metr_norm2[i-1,2]
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t) - 1) / beta_1) * phi0_star * exp(beta_1 * (x - 77.5)) + d * (log(phi0_star) + beta_1 * (x + t - 77.5))) / log(nrow(DS_C)) - 
           sum(- ((exp(beta_1 * t) - 1) / beta_1) * phi0_1    * exp(beta_1 * (x - 77.5)) + d * (log(phi0_1   ) + beta_1 * (x + t - 77.5))) / log(nrow(DS_C))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.008) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.008))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ0_metr_norm2[i, 1] <- phi0
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t) - 1) / beta_star) * phi0 * exp(beta_star * (x - 77.5)) + d * (log(phi0) + beta_star * (x + t - 77.5))) / log(nrow(DS_C)) - 
           sum(- ((exp(beta_1    * t) - 1) / beta_1   ) * phi0 * exp(beta_1    * (x - 77.5)) + d * (log(phi0) + beta_1    * (x + t - 77.5))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ0_metr_norm2[i, 2] <- beta
  
}

##################################################### - J=1 - #################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00003, 0.095, 0.07)
init_zeta1 <- c(0.2, 0.23, 0.01, 0.01, 0.19, 0.12, 0.01, 0.01, 0.01, 0.01, 0.11, 0.09)

#init_tau2 <- c(0.0003, 0.08, 0.9)
#init_zeta2 <- c(0.10, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.1)

# - Set prior distributions
phi_0_prior_shape <- 0.0005
phi_scale <- 1
phi_1_prior_shape <- 0.05

zeta_prior1 <- c(15, 20, 3, 4, 23, 10, 1, 3, 2, 2, 10, 7)

# - Parameter storage

par_postJ1_gibbs_norm1 <- matrix(NA, 100000, 15)
#par_postJ1_gibbs_norm2 <- matrix(NA, 100000, 15)

#par_postJ1_metr_norm1 <- matrix(NA, 100000, 15)
#par_postJ1_metr_norm2 <- matrix(NA, 100000, 15)

par_postJ1_gibbs_norm1[1,c(1:3)] <- init_tau1
par_postJ1_gibbs_norm1[1,c(4:15)] <- init_zeta1

#par_postJ1_gibbs_norm2[1,c(1:3)] <- init_tau2
#par_postJ1_gibbs_norm2[1,c(4:15)] <- init_zeta2

#par_postJ1_metr_norm1[1,c(1:3)] <- init_tau1
#par_postJ1_metr_norm1[1,c(4:15)] <- init_zeta1

#par_postJ1_metr_norm2[1,c(1:3)] <- init_tau2
#par_postJ1_metr_norm2[1,c(4:15)] <- init_zeta2

colnames(par_postJ1_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")
#colnames(par_postJ1_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")
#colnames(par_postJ1_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")
#colnames(par_postJ1_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")

##### - MCMC with Gibbs steps

for(i in 2:100000){
  
  phi0_1 <- par_postJ1_gibbs_norm1[i-1,1]
  beta_1 <- par_postJ1_gibbs_norm1[i-1,2]
  phi1_1 <- par_postJ1_gibbs_norm1[i-1,3]
  
  zeta0_1 <- par_postJ1_gibbs_norm1[i-1,4]
  zeta1_1 <- par_postJ1_gibbs_norm1[i-1,5]
  zeta2_1 <- par_postJ1_gibbs_norm1[i-1,6]
  zeta3_1 <- par_postJ1_gibbs_norm1[i-1,7]
  zeta4_1 <- par_postJ1_gibbs_norm1[i-1,8]
  zeta5_1 <- par_postJ1_gibbs_norm1[i-1,9]
  zeta6_1 <- par_postJ1_gibbs_norm1[i-1,10]
  zeta7_1 <- par_postJ1_gibbs_norm1[i-1,11]
  zeta8_1 <- par_postJ1_gibbs_norm1[i-1,12]
  zeta9_1 <- par_postJ1_gibbs_norm1[i-1,13]
  zeta10_1 <- par_postJ1_gibbs_norm1[i-1,14]
  zeta11_1 <- par_postJ1_gibbs_norm1[i-1,15]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  ## - for dataset with missing geo-demographic profile
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1
  
  ## - for dataset with benefit missing
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  
  c10 <- DS1$G0C0 + DS1$G1C0
  c11 <- DS1$G0C1 + DS1$G1C1
  c12 <- DS1$G0C2 + DS1$G1C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  
  b2L <- DS2$G0BL + DS2$G1BL
  b2H <- DS2$G0BH + DS2$G1BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 12)
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  
  par_postJ1_gibbs_norm1[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ g11)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ g21))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ1_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum(g11 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0) + sum(g21 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0)
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * g11) + sum(d2 * g21)) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1)
  
  par_postJ1_gibbs_norm1[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_gibbs_norm1[i, 2] <- beta
  
}

#for(i in 93914:100000){
  
  phi0_1 <- par_postJ1_gibbs_norm2[i-1,1]
  beta_1 <- par_postJ1_gibbs_norm2[i-1,2]
  phi1_1 <- par_postJ1_gibbs_norm2[i-1,3]
  
  zeta0_1 <- par_postJ1_gibbs_norm2[i-1,4]
  zeta1_1 <- par_postJ1_gibbs_norm2[i-1,5]
  zeta2_1 <- par_postJ1_gibbs_norm2[i-1,6]
  zeta3_1 <- par_postJ1_gibbs_norm2[i-1,7]
  zeta4_1 <- par_postJ1_gibbs_norm2[i-1,8]
  zeta5_1 <- par_postJ1_gibbs_norm2[i-1,9]
  zeta6_1 <- par_postJ1_gibbs_norm2[i-1,10]
  zeta7_1 <- par_postJ1_gibbs_norm2[i-1,11]
  zeta8_1 <- par_postJ1_gibbs_norm2[i-1,12]
  zeta9_1 <- par_postJ1_gibbs_norm2[i-1,13]
  zeta10_1 <- par_postJ1_gibbs_norm2[i-1,14]
  zeta11_1 <- par_postJ1_gibbs_norm2[i-1,15]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  ## - for dataset with missing geo-demographic profile
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1
  
  ## - for dataset with benefit missing
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  
  c10 <- DS1$G0C0 + DS1$G1C0
  c11 <- DS1$G0C1 + DS1$G1C1
  c12 <- DS1$G0C2 + DS1$G1C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  
  b2L <- DS2$G0BL + DS2$G1BL
  b2H <- DS2$G0BH + DS2$G1BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 12)
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  
  par_postJ1_gibbs_norm2[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ g11)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ g21))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ1_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum(g11 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0) + sum(g21 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0)
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * g11) + sum(d2 * g21)) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1)
  
  par_postJ1_gibbs_norm2[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_gibbs_norm2[i, 2] <- beta
  
}

##### - Full metropolis

#for(i in 2:100000){
  
  phi0_1 <- par_postJ1_metr_norm1[i-1,1]
  beta_1 <- par_postJ1_metr_norm1[i-1,2]
  phi1_1 <- par_postJ1_metr_norm1[i-1,3]
  
  zeta0_1 <- par_postJ1_metr_norm1[i-1,4]
  zeta1_1 <- par_postJ1_metr_norm1[i-1,5]
  zeta2_1 <- par_postJ1_metr_norm1[i-1,6]
  zeta3_1 <- par_postJ1_metr_norm1[i-1,7]
  zeta4_1 <- par_postJ1_metr_norm1[i-1,8]
  zeta5_1 <- par_postJ1_metr_norm1[i-1,9]
  zeta6_1 <- par_postJ1_metr_norm1[i-1,10]
  zeta7_1 <- par_postJ1_metr_norm1[i-1,11]
  zeta8_1 <- par_postJ1_metr_norm1[i-1,12]
  zeta9_1 <- par_postJ1_metr_norm1[i-1,13]
  zeta10_1 <- par_postJ1_metr_norm1[i-1,14]
  zeta11_1 <- par_postJ1_metr_norm1[i-1,15]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1
  
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  
  c10 <- DS1$G0C0 + DS1$G1C0
  c11 <- DS1$G0C1 + DS1$G1C1
  c12 <- DS1$G0C2 + DS1$G1C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  
  b2L <- DS2$G0BL + DS2$G1BL
  b2H <- DS2$G0BH + DS2$G1BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 12)
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  
  par_postJ1_metr_norm1[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0112)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.0112) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.0112))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ1_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.17)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_star) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_star) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.17) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.17))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ1_metr_norm1[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_metr_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ1_metr_norm2[i-1,1]
  beta_1 <- par_postJ1_metr_norm2[i-1,2]
  phi1_1 <- par_postJ1_metr_norm2[i-1,3]
  
  zeta0_1 <- par_postJ1_metr_norm2[i-1,4]
  zeta1_1 <- par_postJ1_metr_norm2[i-1,5]
  zeta2_1 <- par_postJ1_metr_norm2[i-1,6]
  zeta3_1 <- par_postJ1_metr_norm2[i-1,7]
  zeta4_1 <- par_postJ1_metr_norm2[i-1,8]
  zeta5_1 <- par_postJ1_metr_norm2[i-1,9]
  zeta6_1 <- par_postJ1_metr_norm2[i-1,10]
  zeta7_1 <- par_postJ1_metr_norm2[i-1,11]
  zeta8_1 <- par_postJ1_metr_norm2[i-1,12]
  zeta9_1 <- par_postJ1_metr_norm2[i-1,13]
  zeta10_1 <- par_postJ1_metr_norm2[i-1,14]
  zeta11_1 <- par_postJ1_metr_norm2[i-1,15]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta0_1 * bL + zeta6_1  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta1_1 * bL + zeta7_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x1 - 77.5))))                 * (zeta2_1 * bL + zeta8_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta3_1 * bL + zeta9_1  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta4_1 * bL + zeta10_1 * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x1 - 77.5)))) * (phi1_1 ^ d1) * (zeta5_1 * bL + zeta11_1 * bH) )
  
  DS1$PostG1C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1
  
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1          * exp(beta_1 * (x2 - 77.5))))                 * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * exp(beta_1 * (x2 - 77.5)))) * (phi1_1 ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2)  )
  
  DS2$PostG1BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  
  c10 <- DS1$G0C0 + DS1$G1C0
  c11 <- DS1$G0C1 + DS1$G1C1
  c12 <- DS1$G0C2 + DS1$G1C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  
  b2L <- DS2$G0BL + DS2$G1BL
  b2H <- DS2$G0BH + DS2$G1BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 12)
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  
  par_postJ1_metr_norm2[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0112)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.0112) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.0112))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ1_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.17)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_star) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_star) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.17) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.17))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ1_metr_norm2[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_metr_norm2[i, 2] <- beta
  
}

################################################# - J=2 - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.05, 0.09, 0.14, 0.05)
init_zeta1 <- c(12, 10, 2, 4, 20, 6, 2, 10, 4, 1, 2, 1, 1, 10, 1, 2, 4, 8)/100

#init_tau2 <- c(0.0005, 0.06, 0.5, 0.15)
#init_zeta2 <- c(0.10, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1)

# - Set prior distributions
phi_0_prior_shape <- 0.09
phi_1_prior_shape <- 0.5
phi_2_prior_shape <- 0.9

phi_scale <- 1

zeta_prior2 <- c(12, 10, 2, 4, 20, 6, 2, 10, 4, 1, 2, 1, 1, 10, 1, 2, 4, 8)
#zeta_prior2 <- c(9, 8, 7, 5, 6, 4, 1, 2, 3, 3, 2, 1, 4, 6, 5, 7, 8, 9)

# - Parameter storage

par_postJ2_gibbs_norm1 <- matrix(NA, 100000, 22)
#par_postJ2_gibbs_norm2 <- matrix(NA, 100000, 22)

#par_postJ2_metr_norm1 <- matrix(NA, 100000, 22)
#par_postJ2_metr_norm2 <- matrix(NA, 100000, 22)

par_postJ2_gibbs_norm1[1,c(1:4)] <- init_tau1
par_postJ2_gibbs_norm1[1,c(5:22)] <- init_zeta1

#par_postJ2_gibbs_norm2[1,c(1:4)] <- init_tau2
#par_postJ2_gibbs_norm2[1,c(5:22)] <- init_zeta2

#par_postJ2_metr_norm1[1,c(1:4)] <- init_tau1
#par_postJ2_metr_norm1[1,c(5:22)] <- init_zeta1

#par_postJ2_metr_norm2[1,c(1:4)] <- init_tau2
#par_postJ2_metr_norm2[1,c(5:22)] <- init_zeta2

colnames(par_postJ2_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")
#colnames(par_postJ2_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")

#colnames(par_postJ2_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")
#colnames(par_postJ2_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")


##### - MCMC with Gibbs steps

for(i in 43959:100000){
  
  phi0_1 <- par_postJ2_gibbs_norm1[i-1,1]
  beta_1 <- par_postJ2_gibbs_norm1[i-1,2]
  phi1_1 <- par_postJ2_gibbs_norm1[i-1,3]
  phi2_1 <- par_postJ2_gibbs_norm1[i-1,4]
  
  zeta0_1 <- par_postJ2_gibbs_norm1[i-1,5]
  zeta1_1 <- par_postJ2_gibbs_norm1[i-1,6]
  zeta2_1 <- par_postJ2_gibbs_norm1[i-1,7]
  zeta3_1 <- par_postJ2_gibbs_norm1[i-1,8]
  zeta4_1 <- par_postJ2_gibbs_norm1[i-1,9]
  zeta5_1 <- par_postJ2_gibbs_norm1[i-1,10]
  zeta6_1 <- par_postJ2_gibbs_norm1[i-1,11]
  zeta7_1 <- par_postJ2_gibbs_norm1[i-1,12]
  zeta8_1 <- par_postJ2_gibbs_norm1[i-1,13]
  zeta9_1 <- par_postJ2_gibbs_norm1[i-1,14]
  zeta10_1 <- par_postJ2_gibbs_norm1[i-1,15]
  zeta11_1 <- par_postJ2_gibbs_norm1[i-1,16]
  zeta12_1 <- par_postJ2_gibbs_norm1[i-1,17]
  zeta13_1 <- par_postJ2_gibbs_norm1[i-1,18]
  zeta14_1 <- par_postJ2_gibbs_norm1[i-1,19]
  zeta15_1 <- par_postJ2_gibbs_norm1[i-1,20]
  zeta16_1 <- par_postJ2_gibbs_norm1[i-1,21]
  zeta17_1 <- par_postJ2_gibbs_norm1[i-1,22]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta0_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta1_1 * bL + zeta10_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta2_1 * bL + zeta11_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta3_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta4_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta5_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta6_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta7_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  
  DS1$PostG2C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1
  
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 18)
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  
  par_postJ2_gibbs_norm1[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), phi_scale + H / log(nrow(DS_C)))
  
  par_postJ2_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5) ) * phi0 * (phi2_1 ^ g12)) + sum( (g21 + g22) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5) ) * phi0 * (phi2_1 ^ g22))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * (g11 + g12)) + sum(d2 * (g21 + g22))) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1)
  
  par_postJ2_gibbs_norm1[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( g12 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 ) + sum( g22 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 )
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + (sum(d1 * g12) + sum(d2 * g22)) / log(nrow(DS_C)), scale=1/(phi_scale + H2 / log(nrow(DS_C))), b=1)
  
  par_postJ2_gibbs_norm1[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_gibbs_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ2_gibbs_norm2[i-1,1]
  beta_1 <- par_postJ2_gibbs_norm2[i-1,2]
  phi1_1 <- par_postJ2_gibbs_norm2[i-1,3]
  phi2_1 <- par_postJ2_gibbs_norm2[i-1,4]
  
  zeta0_1 <- par_postJ2_gibbs_norm2[i-1,5]
  zeta1_1 <- par_postJ2_gibbs_norm2[i-1,6]
  zeta2_1 <- par_postJ2_gibbs_norm2[i-1,7]
  zeta3_1 <- par_postJ2_gibbs_norm2[i-1,8]
  zeta4_1 <- par_postJ2_gibbs_norm2[i-1,9]
  zeta5_1 <- par_postJ2_gibbs_norm2[i-1,10]
  zeta6_1 <- par_postJ2_gibbs_norm2[i-1,11]
  zeta7_1 <- par_postJ2_gibbs_norm2[i-1,12]
  zeta8_1 <- par_postJ2_gibbs_norm2[i-1,13]
  zeta9_1 <- par_postJ2_gibbs_norm2[i-1,14]
  zeta10_1 <- par_postJ2_gibbs_norm2[i-1,15]
  zeta11_1 <- par_postJ2_gibbs_norm2[i-1,16]
  zeta12_1 <- par_postJ2_gibbs_norm2[i-1,17]
  zeta13_1 <- par_postJ2_gibbs_norm2[i-1,18]
  zeta14_1 <- par_postJ2_gibbs_norm2[i-1,19]
  zeta15_1 <- par_postJ2_gibbs_norm2[i-1,20]
  zeta16_1 <- par_postJ2_gibbs_norm2[i-1,21]
  zeta17_1 <- par_postJ2_gibbs_norm2[i-1,22]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta0_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta1_1 * bL + zeta10_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta2_1 * bL + zeta11_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta3_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta4_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta5_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta6_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta7_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  
  DS1$PostG2C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1
  
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 18)
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  
  par_postJ2_gibbs_norm2[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), phi_scale + H / log(nrow(DS_C)))
  
  par_postJ2_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5) ) * phi0 * (phi2_1 ^ g12)) + sum( (g21 + g22) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5) ) * phi0 * (phi2_1 ^ g22))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * (g11 + g12)) + sum(d2 * (g21 + g22))) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1)
  
  par_postJ2_gibbs_norm2[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( g12 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 ) + sum( g22 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 )
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + (sum(d1 * g12) + sum(d2 * g22)) / log(nrow(DS_C)), scale=1/(phi_scale + H2 / log(nrow(DS_C))), b=1)
  
  par_postJ2_gibbs_norm2[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_gibbs_norm2[i, 2] <- beta
  
}


#for(i in 2:100000){
  
  phi0_1 <- par_postJ2_metr_norm1[i-1,1]
  beta_1 <- par_postJ2_metr_norm1[i-1,2]
  phi1_1 <- par_postJ2_metr_norm1[i-1,3]
  phi2_1 <- par_postJ2_metr_norm1[i-1,4]
  
  zeta0_1 <- par_postJ2_metr_norm1[i-1,5]
  zeta1_1 <- par_postJ2_metr_norm1[i-1,6]
  zeta2_1 <- par_postJ2_metr_norm1[i-1,7]
  zeta3_1 <- par_postJ2_metr_norm1[i-1,8]
  zeta4_1 <- par_postJ2_metr_norm1[i-1,9]
  zeta5_1 <- par_postJ2_metr_norm1[i-1,10]
  zeta6_1 <- par_postJ2_metr_norm1[i-1,11]
  zeta7_1 <- par_postJ2_metr_norm1[i-1,12]
  zeta8_1 <- par_postJ2_metr_norm1[i-1,13]
  zeta9_1 <- par_postJ2_metr_norm1[i-1,14]
  zeta10_1 <- par_postJ2_metr_norm1[i-1,15]
  zeta11_1 <- par_postJ2_metr_norm1[i-1,16]
  zeta12_1 <- par_postJ2_metr_norm1[i-1,17]
  zeta13_1 <- par_postJ2_metr_norm1[i-1,18]
  zeta14_1 <- par_postJ2_metr_norm1[i-1,19]
  zeta15_1 <- par_postJ2_metr_norm1[i-1,20]
  zeta16_1 <- par_postJ2_metr_norm1[i-1,21]
  zeta17_1 <- par_postJ2_metr_norm1[i-1,22]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta0_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta1_1 * bL + zeta10_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta2_1 * bL + zeta11_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta3_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta4_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta5_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta6_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta7_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  
  DS1$PostG2C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1
  
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 18)
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  
  par_postJ2_metr_norm1[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.004)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.004) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.004))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ2_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_star) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_1   ) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_star) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_1   ) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.049))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ2_metr_norm1[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_star ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_star) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_1    ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_star ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_star) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_1    ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.112))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi1_unif > r, phi2_1, phi2_star)
  
  par_postJ2_metr_norm1[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_metr_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ2_metr_norm2[i-1,1]
  beta_1 <- par_postJ2_metr_norm2[i-1,2]
  phi1_1 <- par_postJ2_metr_norm2[i-1,3]
  phi2_1 <- par_postJ2_metr_norm2[i-1,4]
  
  zeta0_1 <- par_postJ2_metr_norm2[i-1,5]
  zeta1_1 <- par_postJ2_metr_norm2[i-1,6]
  zeta2_1 <- par_postJ2_metr_norm2[i-1,7]
  zeta3_1 <- par_postJ2_metr_norm2[i-1,8]
  zeta4_1 <- par_postJ2_metr_norm2[i-1,9]
  zeta5_1 <- par_postJ2_metr_norm2[i-1,10]
  zeta6_1 <- par_postJ2_metr_norm2[i-1,11]
  zeta7_1 <- par_postJ2_metr_norm2[i-1,12]
  zeta8_1 <- par_postJ2_metr_norm2[i-1,13]
  zeta9_1 <- par_postJ2_metr_norm2[i-1,14]
  zeta10_1 <- par_postJ2_metr_norm2[i-1,15]
  zeta11_1 <- par_postJ2_metr_norm2[i-1,16]
  zeta12_1 <- par_postJ2_metr_norm2[i-1,17]
  zeta13_1 <- par_postJ2_metr_norm2[i-1,18]
  zeta14_1 <- par_postJ2_metr_norm2[i-1,19]
  zeta15_1 <- par_postJ2_metr_norm2[i-1,20]
  zeta16_1 <- par_postJ2_metr_norm2[i-1,21]
  zeta17_1 <- par_postJ2_metr_norm2[i-1,22]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta0_1 * bL + zeta9_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta1_1 * bL + zeta10_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta2_1 * bL + zeta11_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta3_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta4_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta5_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta6_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta7_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1 + zeta10_1 + zeta11_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) )
  
  
  DS1$PostG2C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1
  
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1)  * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 18)
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  
  par_postJ2_metr_norm2[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.004)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.004) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.004))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ2_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_star) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_1   ) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_star) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_1   ) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.049))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ2_metr_norm2[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_star ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_star) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_1    ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_star ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_star) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_1    ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.112))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi1_unif > r, phi2_1, phi2_star)
  
  par_postJ2_metr_norm2[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_metr_norm2[i, 2] <- beta
  
}


######################################################## - J=3  - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00004, 0.1, 0.02, 0.03, 0.02)
init_zeta1 <- c(12, 10, 2, 4, 10, 4, 1, 10, 4, 1, 10, 2, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 3, 7)/100

#init_tau2 <- c(0.04, 0.06, 0.8, 0.5, 0.6)
#init_zeta2 <- c(0.06, rep(0.04, 22), 0.06)

# - Set prior distributions
phi_0_prior_shape <- 0.09
phi_1_prior_shape <- 0.5
phi_2_prior_shape <- 0.7
phi_3_prior_shape <- 0.9

phi_scale <- 1

zeta_prior3 <- c(12, 10, 2, 4, 10, 4, 1, 10, 4, 1, 10, 2, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 3, 7)

# - Parameter storage

par_postJ3_gibbs_norm1 <- matrix(NA, 100000, 29)
#par_postJ3_gibbs_norm2 <- matrix(NA, 100000, 29)

#par_postJ3_metr_norm1 <- matrix(NA, 100000, 29)
#par_postJ3_metr_norm2 <- matrix(NA, 100000, 29)

par_postJ3_gibbs_norm1[1,c(1:5)] <- init_tau1
par_postJ3_gibbs_norm1[1,c(6:29)] <- init_zeta1

#par_postJ3_gibbs_norm2[1,c(1:5)] <- init_tau2
#par_postJ3_gibbs_norm2[1,c(6:29)] <- init_zeta2

#par_postJ3_metr_norm1[1,c(1:5)] <- init_tau1
#par_postJ3_metr_norm1[1,c(6:29)] <- init_zeta1

#par_postJ3_metr_norm2[1,c(1:5)] <- init_tau2
#par_postJ3_metr_norm2[1,c(6:29)] <- init_zeta2

colnames(par_postJ3_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")
#colnames(par_postJ3_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")

#colnames(par_postJ3_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")
#colnames(par_postJ3_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")


##### - MCMC with Gibbs steps

for(i in 3939:100000){
  
  phi0_1 <- par_postJ3_gibbs_norm1[i-1,1]
  beta_1 <- par_postJ3_gibbs_norm1[i-1,2]
  phi1_1 <- par_postJ3_gibbs_norm1[i-1,3]
  phi2_1 <- par_postJ3_gibbs_norm1[i-1,4]
  phi3_1 <- par_postJ3_gibbs_norm1[i-1,5]
  
  
  zeta0_1 <- par_postJ3_gibbs_norm1[i-1,6]
  zeta1_1 <- par_postJ3_gibbs_norm1[i-1,7]
  zeta2_1 <- par_postJ3_gibbs_norm1[i-1,8]
  zeta3_1 <- par_postJ3_gibbs_norm1[i-1,9]
  zeta4_1 <- par_postJ3_gibbs_norm1[i-1,10]
  zeta5_1 <- par_postJ3_gibbs_norm1[i-1,11]
  zeta6_1 <- par_postJ3_gibbs_norm1[i-1,12]
  zeta7_1 <- par_postJ3_gibbs_norm1[i-1,13]
  zeta8_1 <- par_postJ3_gibbs_norm1[i-1,14]
  zeta9_1 <- par_postJ3_gibbs_norm1[i-1,15]
  zeta10_1 <- par_postJ3_gibbs_norm1[i-1,16]
  zeta11_1 <- par_postJ3_gibbs_norm1[i-1,17]
  zeta12_1 <- par_postJ3_gibbs_norm1[i-1,18]
  zeta13_1 <- par_postJ3_gibbs_norm1[i-1,19]
  zeta14_1 <- par_postJ3_gibbs_norm1[i-1,20]
  zeta15_1 <- par_postJ3_gibbs_norm1[i-1,21]
  zeta16_1 <- par_postJ3_gibbs_norm1[i-1,22]
  zeta17_1 <- par_postJ3_gibbs_norm1[i-1,23]
  zeta18_1 <- par_postJ3_gibbs_norm1[i-1,24]
  zeta19_1 <- par_postJ3_gibbs_norm1[i-1,25]
  zeta20_1 <- par_postJ3_gibbs_norm1[i-1,26]
  zeta21_1 <- par_postJ3_gibbs_norm1[i-1,27]
  zeta22_1 <- par_postJ3_gibbs_norm1[i-1,28]
  zeta23_1 <- par_postJ3_gibbs_norm1[i-1,29]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta0_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta1_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta2_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta3_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta4_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta5_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta6_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta7_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta8_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta9_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta10_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1- DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL
  
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 24)
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  
  par_postJ3_gibbs_norm1[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ3_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum( (g21 + g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * (g11 + g12 + g13)) + sum(d2 * (g21 + g22 + g23))) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ3_gibbs_norm1[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g13)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g23))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + (sum(d1 * (g12 + g13)) + sum(d2 * (g22 + g23))) / log(nrow(DS_C)), scale=1/(phi_scale + H2 / log(nrow(DS_C))), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_2_prior_shape + sum(d * (g2 + g3)), rate = phi_scale + H2)#rgamma(1, phi_2_prior_shape + sum(d * (g2 + g3)), phi_scale + H2)
  
  par_postJ3_gibbs_norm1[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( g13 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2) + sum( g23 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2)
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + (sum(d1 * g13) + sum(d2 * g23)) / log(nrow(DS_C)), scale=1/(phi_scale + H3 / log(nrow(DS_C))), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ3_gibbs_norm1[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_gibbs_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ3_gibbs_norm2[i-1,1]
  beta_1 <- par_postJ3_gibbs_norm2[i-1,2]
  phi1_1 <- par_postJ3_gibbs_norm2[i-1,3]
  phi2_1 <- par_postJ3_gibbs_norm2[i-1,4]
  phi3_1 <- par_postJ3_gibbs_norm2[i-1,5]
  
  
  zeta0_1 <- par_postJ3_gibbs_norm2[i-1,6]
  zeta1_1 <- par_postJ3_gibbs_norm2[i-1,7]
  zeta2_1 <- par_postJ3_gibbs_norm2[i-1,8]
  zeta3_1 <- par_postJ3_gibbs_norm2[i-1,9]
  zeta4_1 <- par_postJ3_gibbs_norm2[i-1,10]
  zeta5_1 <- par_postJ3_gibbs_norm2[i-1,11]
  zeta6_1 <- par_postJ3_gibbs_norm2[i-1,12]
  zeta7_1 <- par_postJ3_gibbs_norm2[i-1,13]
  zeta8_1 <- par_postJ3_gibbs_norm2[i-1,14]
  zeta9_1 <- par_postJ3_gibbs_norm2[i-1,15]
  zeta10_1 <- par_postJ3_gibbs_norm2[i-1,16]
  zeta11_1 <- par_postJ3_gibbs_norm2[i-1,17]
  zeta12_1 <- par_postJ3_gibbs_norm2[i-1,18]
  zeta13_1 <- par_postJ3_gibbs_norm2[i-1,19]
  zeta14_1 <- par_postJ3_gibbs_norm2[i-1,20]
  zeta15_1 <- par_postJ3_gibbs_norm2[i-1,21]
  zeta16_1 <- par_postJ3_gibbs_norm2[i-1,22]
  zeta17_1 <- par_postJ3_gibbs_norm2[i-1,23]
  zeta18_1 <- par_postJ3_gibbs_norm2[i-1,24]
  zeta19_1 <- par_postJ3_gibbs_norm2[i-1,25]
  zeta20_1 <- par_postJ3_gibbs_norm2[i-1,26]
  zeta21_1 <- par_postJ3_gibbs_norm2[i-1,27]
  zeta22_1 <- par_postJ3_gibbs_norm2[i-1,28]
  zeta23_1 <- par_postJ3_gibbs_norm2[i-1,29]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta0_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta1_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta2_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta3_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta4_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta5_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta6_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta7_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta8_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta9_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta10_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1- DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL
  
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 24)
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  
  par_postJ3_gibbs_norm2[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ3_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum( (g21 + g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * (g11 + g12 + g13)) + sum(d2 * (g21 + g22 + g23))) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ3_gibbs_norm2[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g13)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g23))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + (sum(d1 * (g12 + g13)) + sum(d2 * (g22 + g23))) / log(nrow(DS_C)), scale=1/(phi_scale + H2 / log(nrow(DS_C))), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_2_prior_shape + sum(d * (g2 + g3)), rate = phi_scale + H2)#rgamma(1, phi_2_prior_shape + sum(d * (g2 + g3)), phi_scale + H2)
  
  par_postJ3_gibbs_norm2[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( g13 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2) + sum( g23 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2)
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + (sum(d1 * g13) + sum(d2 * g23)) / log(nrow(DS_C)), scale=1/(phi_scale + H3 / log(nrow(DS_C))), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ3_gibbs_norm2[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_gibbs_norm2[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ3_metr_norm1[i-1,1]
  beta_1 <- par_postJ3_metr_norm1[i-1,2]
  phi1_1 <- par_postJ3_metr_norm1[i-1,3]
  phi2_1 <- par_postJ3_metr_norm1[i-1,4]
  phi3_1 <- par_postJ3_metr_norm1[i-1,5]
  
  
  zeta0_1 <- par_postJ3_metr_norm1[i-1,6]
  zeta1_1 <- par_postJ3_metr_norm1[i-1,7]
  zeta2_1 <- par_postJ3_metr_norm1[i-1,8]
  zeta3_1 <- par_postJ3_metr_norm1[i-1,9]
  zeta4_1 <- par_postJ3_metr_norm1[i-1,10]
  zeta5_1 <- par_postJ3_metr_norm1[i-1,11]
  zeta6_1 <- par_postJ3_metr_norm1[i-1,12]
  zeta7_1 <- par_postJ3_metr_norm1[i-1,13]
  zeta8_1 <- par_postJ3_metr_norm1[i-1,14]
  zeta9_1 <- par_postJ3_metr_norm1[i-1,15]
  zeta10_1 <- par_postJ3_metr_norm1[i-1,16]
  zeta11_1 <- par_postJ3_metr_norm1[i-1,17]
  zeta12_1 <- par_postJ3_metr_norm1[i-1,18]
  zeta13_1 <- par_postJ3_metr_norm1[i-1,19]
  zeta14_1 <- par_postJ3_metr_norm1[i-1,20]
  zeta15_1 <- par_postJ3_metr_norm1[i-1,21]
  zeta16_1 <- par_postJ3_metr_norm1[i-1,22]
  zeta17_1 <- par_postJ3_metr_norm1[i-1,23]
  zeta18_1 <- par_postJ3_metr_norm1[i-1,24]
  zeta19_1 <- par_postJ3_metr_norm1[i-1,25]
  zeta20_1 <- par_postJ3_metr_norm1[i-1,26]
  zeta21_1 <- par_postJ3_metr_norm1[i-1,27]
  zeta22_1 <- par_postJ3_metr_norm1[i-1,28]
  zeta23_1 <- par_postJ3_metr_norm1[i-1,29]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta0_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta1_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta2_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta3_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta4_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta5_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta6_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta7_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta8_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta9_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta10_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1- DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL
  
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 24)
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  
  par_postJ3_metr_norm1[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0034)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi0_star, a=0, mean = phi0_1, sd = 0.0034) / dtruncnorm(phi0_1, a=0, mean = phi0_star, sd = 0.0034))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)  
  
  par_postJ3_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.058)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_star) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_1   ) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_star) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_1   ) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi1_star, a=0, mean = phi1_1, sd = 0.058) / dtruncnorm(phi1_1, a=0, mean = phi1_star, sd = 0.058))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)  
  
  par_postJ3_metr_norm1[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.115)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_star ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_star) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_1    ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_1   ) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_star ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_star) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_1    ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_1   ) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi2_star, a=0, mean = phi2_1, sd = 0.115) / dtruncnorm(phi2_1, a=0, mean = phi2_star, sd = 0.115))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)  
  
  par_postJ3_metr_norm1[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.125)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_star ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_star) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_1    ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_star ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_star) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_1    ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi3_star, a=0, mean = phi3_1, sd = 0.125) / dtruncnorm(phi3_1, a=0, mean = phi3_star, sd = 0.125))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)  
  
  par_postJ3_metr_norm1[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_metr_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ3_metr_norm2[i-1,1]
  beta_1 <- par_postJ3_metr_norm2[i-1,2]
  phi1_1 <- par_postJ3_metr_norm2[i-1,3]
  phi2_1 <- par_postJ3_metr_norm2[i-1,4]
  phi3_1 <- par_postJ3_metr_norm2[i-1,5]
  
  
  zeta0_1 <- par_postJ3_metr_norm2[i-1,6]
  zeta1_1 <- par_postJ3_metr_norm2[i-1,7]
  zeta2_1 <- par_postJ3_metr_norm2[i-1,8]
  zeta3_1 <- par_postJ3_metr_norm2[i-1,9]
  zeta4_1 <- par_postJ3_metr_norm2[i-1,10]
  zeta5_1 <- par_postJ3_metr_norm2[i-1,11]
  zeta6_1 <- par_postJ3_metr_norm2[i-1,12]
  zeta7_1 <- par_postJ3_metr_norm2[i-1,13]
  zeta8_1 <- par_postJ3_metr_norm2[i-1,14]
  zeta9_1 <- par_postJ3_metr_norm2[i-1,15]
  zeta10_1 <- par_postJ3_metr_norm2[i-1,16]
  zeta11_1 <- par_postJ3_metr_norm2[i-1,17]
  zeta12_1 <- par_postJ3_metr_norm2[i-1,18]
  zeta13_1 <- par_postJ3_metr_norm2[i-1,19]
  zeta14_1 <- par_postJ3_metr_norm2[i-1,20]
  zeta15_1 <- par_postJ3_metr_norm2[i-1,21]
  zeta16_1 <- par_postJ3_metr_norm2[i-1,22]
  zeta17_1 <- par_postJ3_metr_norm2[i-1,23]
  zeta18_1 <- par_postJ3_metr_norm2[i-1,24]
  zeta19_1 <- par_postJ3_metr_norm2[i-1,25]
  zeta20_1 <- par_postJ3_metr_norm2[i-1,26]
  zeta21_1 <- par_postJ3_metr_norm2[i-1,27]
  zeta22_1 <- par_postJ3_metr_norm2[i-1,28]
  zeta23_1 <- par_postJ3_metr_norm2[i-1,29]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta0_1 * bL + zeta12_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta1_1 * bL + zeta13_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * (zeta2_1 * bL + zeta14_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta3_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta4_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * (zeta5_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta6_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta7_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * (zeta8_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta9_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * (zeta10_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x1 - 77.5))))                                     * ((zeta0_1 + zeta1_1  + zeta2_1 ) * bL + (zeta12_1 + zeta13_1 + zeta14_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                    ^ d1) * ((zeta3_1 + zeta4_1  + zeta5_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d1) * ((zeta6_1 + zeta7_1  + zeta8_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d1) * ((zeta9_1 + zeta10_1 + zeta11_1) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) )
  
  DS1$PostG3C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1- DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                            * exp(beta_1 * (x2 - 77.5))))                                     * ((zeta0_1 + zeta12_1 ) * c0 + (zeta1_1  + zeta13_1) * c1 + (zeta2_1  + zeta14_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                   * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                    ^ d2) * ((zeta3_1 + zeta15_1) * c0 + (zeta4_1  + zeta16_1) * c1 + (zeta5_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1         ) ^ d2) * ((zeta6_1 + zeta18_1) * c0 + (zeta7_1  + zeta19_1) * c1 + (zeta8_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1) ^ d2) * ((zeta9_1 + zeta21_1) * c0 + (zeta10_1 + zeta22_1) * c1 + (zeta11_1 + zeta23_1) * c2))
  
  DS2$PostG3BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL
  
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 24)
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  
  par_postJ3_metr_norm2[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0034)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi0_star, a=0, mean = phi0_1, sd = 0.0034) / dtruncnorm(phi0_1, a=0, mean = phi0_star, sd = 0.0034))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)  
  
  par_postJ3_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.058)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_star) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_1   ) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_star) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_1   ) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi1_star, a=0, mean = phi1_1, sd = 0.058) / dtruncnorm(phi1_1, a=0, mean = phi1_star, sd = 0.058))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)  
  
  par_postJ3_metr_norm2[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.115)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_star ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_star) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_1    ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_1   ) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_star ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_star) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_1    ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_1   ) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi2_star, a=0, mean = phi2_1, sd = 0.115) / dtruncnorm(phi2_1, a=0, mean = phi2_star, sd = 0.115))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)  
  
  par_postJ3_metr_norm2[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.125)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_star ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_star) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_1    ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_star ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_star) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_1    ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(phi3_star, a=0, mean = phi3_1, sd = 0.125) / dtruncnorm(phi3_1, a=0, mean = phi3_star, sd = 0.125))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)  
  
  par_postJ3_metr_norm2[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_metr_norm2[i, 2] <- beta
  
}

######################################################## - J=4  - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00004, 0.1, 0.02, 0.03, 0.02, 0.9)
init_zeta1 <- c(10, 5, 1, 4, 10, 2, 4, 10, 4, 0.5, 10, 3, 0.5, 5, 2, 1, 2, 0.5, 1, 2, 0.5, 1, 4, 1.5, .5, 4, 2.5, 0.5, 4, 5)/100

#init_tau2 <- c(0.04, 0.06, 0.8, 0.5, 0.6, 0.04)
#init_zeta2 <- c(0.08, rep(0.04, 28), 0.08)

# - Set prior distributions
phi_0_prior_shape <- 0.09
phi_1_prior_shape <- 0.5
phi_2_prior_shape <- 0.75
phi_3_prior_shape <- 0.85
phi_4_prior_shape <- 0.9

phi_scale <- 1

zeta_prior4 <- c(10, 5, 1, 4, 10, 2, 4, 10, 4, 0.5, 10, 3, 0.5, 5, 2, 1, 2, 0.5, 1, 2, 0.5, 1, 4, 1.5, .5, 4, 2.5, 0.5, 4, 5)

# - Parameter storage

par_postJ4_gibbs_norm1 <- matrix(NA, 100000, 36)
#par_postJ4_gibbs_norm2 <- matrix(NA, 100000, 36)

#par_postJ4_metr_norm1 <- matrix(NA, 100000, 36)
#par_postJ4_metr_norm2 <- matrix(NA, 100000, 36)

par_postJ4_gibbs_norm1[1,c(1:6)] <- init_tau1
par_postJ4_gibbs_norm1[1,c(7:36)] <- init_zeta1

#par_postJ4_gibbs_norm2[1,c(1:6)] <- init_tau2
#par_postJ4_gibbs_norm2[1,c(7:36)] <- init_zeta2

#par_postJ4_metr_norm1[1,c(1:6)] <- init_tau1
#par_postJ4_metr_norm1[1,c(7:36)] <- init_zeta1

#par_postJ4_metr_norm2[1,c(1:6)] <- init_tau2
#par_postJ4_metr_norm2[1,c(7:36)] <- init_zeta2

colnames(par_postJ4_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")
#colnames(par_postJ4_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")

#colnames(par_postJ4_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")
#colnames(par_postJ4_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")


##### - MCMC with Gibbs steps

for(i in 2:100000){
  
  phi0_1 <- par_postJ4_gibbs_norm1[i-1,1]
  beta_1 <- par_postJ4_gibbs_norm1[i-1,2]
  phi1_1 <- par_postJ4_gibbs_norm1[i-1,3]
  phi2_1 <- par_postJ4_gibbs_norm1[i-1,4]
  phi3_1 <- par_postJ4_gibbs_norm1[i-1,5]
  phi4_1 <- par_postJ4_gibbs_norm1[i-1,6]
  
  zeta0_1 <- par_postJ4_gibbs_norm1[i-1,7]
  zeta1_1 <- par_postJ4_gibbs_norm1[i-1,8]
  zeta2_1 <- par_postJ4_gibbs_norm1[i-1,9]
  zeta3_1 <- par_postJ4_gibbs_norm1[i-1,10]
  zeta4_1 <- par_postJ4_gibbs_norm1[i-1,11]
  zeta5_1 <- par_postJ4_gibbs_norm1[i-1,12]
  zeta6_1 <- par_postJ4_gibbs_norm1[i-1,13]
  zeta7_1 <- par_postJ4_gibbs_norm1[i-1,14]
  zeta8_1 <- par_postJ4_gibbs_norm1[i-1,15]
  zeta9_1 <- par_postJ4_gibbs_norm1[i-1,16]
  zeta10_1 <- par_postJ4_gibbs_norm1[i-1,17]
  zeta11_1 <- par_postJ4_gibbs_norm1[i-1,18]
  zeta12_1 <- par_postJ4_gibbs_norm1[i-1,19]
  zeta13_1 <- par_postJ4_gibbs_norm1[i-1,20]
  zeta14_1 <- par_postJ4_gibbs_norm1[i-1,21]
  zeta15_1 <- par_postJ4_gibbs_norm1[i-1,22]
  zeta16_1 <- par_postJ4_gibbs_norm1[i-1,23]
  zeta17_1 <- par_postJ4_gibbs_norm1[i-1,24]
  zeta18_1 <- par_postJ4_gibbs_norm1[i-1,25]
  zeta19_1 <- par_postJ4_gibbs_norm1[i-1,26]
  zeta20_1 <- par_postJ4_gibbs_norm1[i-1,27]
  zeta21_1 <- par_postJ4_gibbs_norm1[i-1,28]
  zeta22_1 <- par_postJ4_gibbs_norm1[i-1,29]
  zeta23_1 <- par_postJ4_gibbs_norm1[i-1,30]
  zeta24_1 <- par_postJ4_gibbs_norm1[i-1,31]
  zeta25_1 <- par_postJ4_gibbs_norm1[i-1,32]
  zeta26_1 <- par_postJ4_gibbs_norm1[i-1,33]
  zeta27_1 <- par_postJ4_gibbs_norm1[i-1,34]
  zeta28_1 <- par_postJ4_gibbs_norm1[i-1,35]
  zeta29_1 <- par_postJ4_gibbs_norm1[i-1,36]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta0_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta1_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta2_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta3_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta4_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta5_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta6_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta7_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta8_1 * bL + zeta23_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta11_1 * bL + zeta26_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1 - DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1 - DS1$PostG3C2 - DS1$PostG4C0 - DS1$PostG4C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta21_1 * c0 + zeta22_1 * c1 + zeta23_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta24_1 * c0 + zeta25_1 * c1 + zeta26_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL - DS2$PostG3BH - DS2$PostG4BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2)), 1, 0)
  DS1$G4C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0)), 1, 0)
  DS1$G4C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0 + DS1$PostG4C1)), 1, 0)
  DS1$G4C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1 - DS1$G3C2 - DS1$G4C0 - DS1$G4C1
  
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  g14 <- DS1$G4C0 + DS1$G4C1 + DS1$G4C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0 + DS1$G4C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1 + DS1$G4C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2 + DS1$G4C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH)), 1, 0)
  DS2$G4BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH + DS2$PostG4BL)), 1, 0)
  DS2$G4BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL - DS2$G3BH - DS2$G4BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  g24 <- DS2$G4BL + DS2$G4BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL + DS2$G4BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH + DS2$G4BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 30)
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + (sum(bL * c10 * g14) + sum(b2L * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + (sum(bL * c11 * g14) + sum(b2L * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + (sum(bL * c12 * g14) + sum(b2L * c2 * g24)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + (sum(bH * c10 * g14) + sum(b2H * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + (sum(bH * c11 * g14) + sum(b2H * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + (sum(bH * c12 * g14) + sum(b2H * c2 * g24)) / log(nrow(DS_C)), 1)
  
  par_postJ4_gibbs_norm1[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ4_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g21 + g22 + g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * (g11 + g12 + g13 + g14)) + sum(d2 * (g21 + g22 + g23 + g24))) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ4_gibbs_norm1[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + (sum(d1 * (g12 + g13 + g14)) + sum(d2 * (g22 + g23 + g24))) / log(nrow(DS_C)), scale=1/(phi_scale + H2 / log(nrow(DS_C))), b=1)
  
  par_postJ4_gibbs_norm1[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( (g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g14)) + sum( (g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g24))
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + (sum(d1 * (g13 + g14)) + sum(d2 * (g23 + g24))) / log(nrow(DS_C)), scale=1/(phi_scale + H3 / log(nrow(DS_C))), b=1)
  
  par_postJ4_gibbs_norm1[i, 5] <- phi3
  
  # - phi_4
  
  H4 <- sum( g14 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * phi3) + sum( g24 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * phi3)
  
  phi4 <- rtgamma(1, shape=phi_4_prior_shape + (sum(d1 * g14) + sum(d2 * g24)) / log(nrow(DS_C)), scale=1/(phi_scale + H4 / log(nrow(DS_C))), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ4_gibbs_norm1[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.025)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.025) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.025))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_gibbs_norm1[i, 2] <- beta
  
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ4_gibbs_norm2[i-1,1]
  beta_1 <- par_postJ4_gibbs_norm2[i-1,2]
  phi1_1 <- par_postJ4_gibbs_norm2[i-1,3]
  phi2_1 <- par_postJ4_gibbs_norm2[i-1,4]
  phi3_1 <- par_postJ4_gibbs_norm2[i-1,5]
  phi4_1 <- par_postJ4_gibbs_norm2[i-1,6]
  
  zeta0_1 <- par_postJ4_gibbs_norm2[i-1,7]
  zeta1_1 <- par_postJ4_gibbs_norm2[i-1,8]
  zeta2_1 <- par_postJ4_gibbs_norm2[i-1,9]
  zeta3_1 <- par_postJ4_gibbs_norm2[i-1,10]
  zeta4_1 <- par_postJ4_gibbs_norm2[i-1,11]
  zeta5_1 <- par_postJ4_gibbs_norm2[i-1,12]
  zeta6_1 <- par_postJ4_gibbs_norm2[i-1,13]
  zeta7_1 <- par_postJ4_gibbs_norm2[i-1,14]
  zeta8_1 <- par_postJ4_gibbs_norm2[i-1,15]
  zeta9_1 <- par_postJ4_gibbs_norm2[i-1,16]
  zeta10_1 <- par_postJ4_gibbs_norm2[i-1,17]
  zeta11_1 <- par_postJ4_gibbs_norm2[i-1,18]
  zeta12_1 <- par_postJ4_gibbs_norm2[i-1,19]
  zeta13_1 <- par_postJ4_gibbs_norm2[i-1,20]
  zeta14_1 <- par_postJ4_gibbs_norm2[i-1,21]
  zeta15_1 <- par_postJ4_gibbs_norm2[i-1,22]
  zeta16_1 <- par_postJ4_gibbs_norm2[i-1,23]
  zeta17_1 <- par_postJ4_gibbs_norm2[i-1,24]
  zeta18_1 <- par_postJ4_gibbs_norm2[i-1,25]
  zeta19_1 <- par_postJ4_gibbs_norm2[i-1,26]
  zeta20_1 <- par_postJ4_gibbs_norm2[i-1,27]
  zeta21_1 <- par_postJ4_gibbs_norm2[i-1,28]
  zeta22_1 <- par_postJ4_gibbs_norm2[i-1,29]
  zeta23_1 <- par_postJ4_gibbs_norm2[i-1,30]
  zeta24_1 <- par_postJ4_gibbs_norm2[i-1,31]
  zeta25_1 <- par_postJ4_gibbs_norm2[i-1,32]
  zeta26_1 <- par_postJ4_gibbs_norm2[i-1,33]
  zeta27_1 <- par_postJ4_gibbs_norm2[i-1,34]
  zeta28_1 <- par_postJ4_gibbs_norm2[i-1,35]
  zeta29_1 <- par_postJ4_gibbs_norm2[i-1,36]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta0_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta1_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta2_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta3_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta4_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta5_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta6_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta7_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta8_1 * bL + zeta23_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta11_1 * bL + zeta26_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1 - DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1 - DS1$PostG3C2 - DS1$PostG4C0 - DS1$PostG4C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta21_1 * c0 + zeta22_1 * c1 + zeta23_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta24_1 * c0 + zeta25_1 * c1 + zeta26_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL - DS2$PostG3BH - DS2$PostG4BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2)), 1, 0)
  DS1$G4C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0)), 1, 0)
  DS1$G4C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0 + DS1$PostG4C1)), 1, 0)
  DS1$G4C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1 - DS1$G3C2 - DS1$G4C0 - DS1$G4C1
  
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  g14 <- DS1$G4C0 + DS1$G4C1 + DS1$G4C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0 + DS1$G4C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1 + DS1$G4C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2 + DS1$G4C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH)), 1, 0)
  DS2$G4BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH + DS2$PostG4BL)), 1, 0)
  DS2$G4BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL - DS2$G3BH - DS2$G4BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  g24 <- DS2$G4BL + DS2$G4BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL + DS2$G4BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH + DS2$G4BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 30)
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + (sum(bL * c10 * g14) + sum(b2L * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + (sum(bL * c11 * g14) + sum(b2L * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + (sum(bL * c12 * g14) + sum(b2L * c2 * g24)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + (sum(bH * c10 * g14) + sum(b2H * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + (sum(bH * c11 * g14) + sum(b2H * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + (sum(bH * c12 * g14) + sum(b2H * c2 * g24)) / log(nrow(DS_C)), 1)
  
  par_postJ4_gibbs_norm2[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi0 <- rgamma(1, phi_0_prior_shape + (sum(d1) + sum(d2)) / log(nrow(DS_C)), scale=1/(phi_scale + H / log(nrow(DS_C))))
  
  par_postJ4_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g21 + g22 + g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + (sum(d1 * (g11 + g12 + g13 + g14)) + sum(d2 * (g21 + g22 + g23 + g24))) / log(nrow(DS_C)), scale=1/(phi_scale + H1 / log(nrow(DS_C))), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ4_gibbs_norm2[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + (sum(d1 * (g12 + g13 + g14)) + sum(d2 * (g22 + g23 + g24))) / log(nrow(DS_C)), scale=1/(phi_scale + H2 / log(nrow(DS_C))), b=1)
  
  par_postJ4_gibbs_norm2[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( (g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g14)) + sum( (g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g24))
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + (sum(d1 * (g13 + g14)) + sum(d2 * (g23 + g24))) / log(nrow(DS_C)), scale=1/(phi_scale + H3 / log(nrow(DS_C))), b=1)
  
  par_postJ4_gibbs_norm2[i, 5] <- phi3
  
  # - phi_4
  
  H4 <- sum( g14 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * phi3) + sum( g24 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * phi3)
  
  phi4 <- rtgamma(1, shape=phi_4_prior_shape + (sum(d1 * g14) + sum(d2 * g24)) / log(nrow(DS_C)), scale=1/(phi_scale + H4 / log(nrow(DS_C))), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ4_gibbs_norm2[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_gibbs_norm2[i, 2] <- beta
  
  
}


#for(i in 2:100000){
  
  phi0_1 <- par_postJ4_metr_norm1[i-1,1]
  beta_1 <- par_postJ4_metr_norm1[i-1,2]
  phi1_1 <- par_postJ4_metr_norm1[i-1,3]
  phi2_1 <- par_postJ4_metr_norm1[i-1,4]
  phi3_1 <- par_postJ4_metr_norm1[i-1,5]
  phi4_1 <- par_postJ4_metr_norm1[i-1,6]
  
  zeta0_1 <- par_postJ4_metr_norm1[i-1,7]
  zeta1_1 <- par_postJ4_metr_norm1[i-1,8]
  zeta2_1 <- par_postJ4_metr_norm1[i-1,9]
  zeta3_1 <- par_postJ4_metr_norm1[i-1,10]
  zeta4_1 <- par_postJ4_metr_norm1[i-1,11]
  zeta5_1 <- par_postJ4_metr_norm1[i-1,12]
  zeta6_1 <- par_postJ4_metr_norm1[i-1,13]
  zeta7_1 <- par_postJ4_metr_norm1[i-1,14]
  zeta8_1 <- par_postJ4_metr_norm1[i-1,15]
  zeta9_1 <- par_postJ4_metr_norm1[i-1,16]
  zeta10_1 <- par_postJ4_metr_norm1[i-1,17]
  zeta11_1 <- par_postJ4_metr_norm1[i-1,18]
  zeta12_1 <- par_postJ4_metr_norm1[i-1,19]
  zeta13_1 <- par_postJ4_metr_norm1[i-1,20]
  zeta14_1 <- par_postJ4_metr_norm1[i-1,21]
  zeta15_1 <- par_postJ4_metr_norm1[i-1,22]
  zeta16_1 <- par_postJ4_metr_norm1[i-1,23]
  zeta17_1 <- par_postJ4_metr_norm1[i-1,24]
  zeta18_1 <- par_postJ4_metr_norm1[i-1,25]
  zeta19_1 <- par_postJ4_metr_norm1[i-1,26]
  zeta20_1 <- par_postJ4_metr_norm1[i-1,27]
  zeta21_1 <- par_postJ4_metr_norm1[i-1,28]
  zeta22_1 <- par_postJ4_metr_norm1[i-1,29]
  zeta23_1 <- par_postJ4_metr_norm1[i-1,30]
  zeta24_1 <- par_postJ4_metr_norm1[i-1,31]
  zeta25_1 <- par_postJ4_metr_norm1[i-1,32]
  zeta26_1 <- par_postJ4_metr_norm1[i-1,33]
  zeta27_1 <- par_postJ4_metr_norm1[i-1,34]
  zeta28_1 <- par_postJ4_metr_norm1[i-1,35]
  zeta29_1 <- par_postJ4_metr_norm1[i-1,36]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta0_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta1_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta2_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta3_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta4_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta5_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta6_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta7_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta8_1 * bL + zeta23_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta11_1 * bL + zeta26_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1 - DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1 - DS1$PostG3C2 - DS1$PostG4C0 - DS1$PostG4C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta21_1 * c0 + zeta22_1 * c1 + zeta23_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta24_1 * c0 + zeta25_1 * c1 + zeta26_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL - DS2$PostG3BH - DS2$PostG4BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2)), 1, 0)
  DS1$G4C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0)), 1, 0)
  DS1$G4C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0 + DS1$PostG4C1)), 1, 0)
  DS1$G4C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1 - DS1$G3C2 - DS1$G4C0 - DS1$G4C1
  
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  g14 <- DS1$G4C0 + DS1$G4C1 + DS1$G4C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0 + DS1$G4C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1 + DS1$G4C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2 + DS1$G4C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH)), 1, 0)
  DS2$G4BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH + DS2$PostG4BL)), 1, 0)
  DS2$G4BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL - DS2$G3BH - DS2$G4BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  g24 <- DS2$G4BL + DS2$G4BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL + DS2$G4BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH + DS2$G4BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 30)
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + (sum(bL * c10 * g14) + sum(b2L * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + (sum(bL * c11 * g14) + sum(b2L * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + (sum(bL * c12 * g14) + sum(b2L * c2 * g24)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + (sum(bH * c10 * g14) + sum(b2H * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + (sum(bH * c11 * g14) + sum(b2H * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + (sum(bH * c12 * g14) + sum(b2H * c2 * g24)) / log(nrow(DS_C)), 1)
  
  par_postJ4_metr_norm1[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.00325)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.00325) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.00325))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ4_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.09)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_star) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_1   ) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_star) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_1   ) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.09) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.09))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ4_metr_norm1[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.11)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_star ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_star) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_1    ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_1   ) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_star ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_star) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_1    ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_1   ) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.11) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.11))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)
  
  par_postJ4_metr_norm1[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.0975)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_star ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_star) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_1    ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_1   ) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_star ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_star) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_1    ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_1   ) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi3_star, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.0975) / dtruncnorm(phi3_1, a=(10^(-20)), b=1, mean = phi3_star, sd = 0.0975))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)
  
  par_postJ4_metr_norm1[i, 5] <- phi3
  
  # - phi_4
  
  phi4_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.18)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_star ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_star) + beta_1 * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_1    ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_star ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_star) + beta_1 * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_1    ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi4_star, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.18) / dtruncnorm(phi4_1, a=(10^(-20)), b=1, mean = phi4_star, sd = 0.18))
  
  phi4_unif <- runif(1, 0, 1)
  
  phi4 <- ifelse( phi4_unif > r, phi4_1, phi4_star)
  
  par_postJ4_metr_norm1[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
            sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
            sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
            sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_metr_norm1[i, 2] <- beta
  
}

#for(i in 2:100000){
  
  phi0_1 <- par_postJ4_metr_norm2[i-1,1]
  beta_1 <- par_postJ4_metr_norm2[i-1,2]
  phi1_1 <- par_postJ4_metr_norm2[i-1,3]
  phi2_1 <- par_postJ4_metr_norm2[i-1,4]
  phi3_1 <- par_postJ4_metr_norm2[i-1,5]
  phi4_1 <- par_postJ4_metr_norm2[i-1,6]
  
  zeta0_1 <- par_postJ4_metr_norm2[i-1,7]
  zeta1_1 <- par_postJ4_metr_norm2[i-1,8]
  zeta2_1 <- par_postJ4_metr_norm2[i-1,9]
  zeta3_1 <- par_postJ4_metr_norm2[i-1,10]
  zeta4_1 <- par_postJ4_metr_norm2[i-1,11]
  zeta5_1 <- par_postJ4_metr_norm2[i-1,12]
  zeta6_1 <- par_postJ4_metr_norm2[i-1,13]
  zeta7_1 <- par_postJ4_metr_norm2[i-1,14]
  zeta8_1 <- par_postJ4_metr_norm2[i-1,15]
  zeta9_1 <- par_postJ4_metr_norm2[i-1,16]
  zeta10_1 <- par_postJ4_metr_norm2[i-1,17]
  zeta11_1 <- par_postJ4_metr_norm2[i-1,18]
  zeta12_1 <- par_postJ4_metr_norm2[i-1,19]
  zeta13_1 <- par_postJ4_metr_norm2[i-1,20]
  zeta14_1 <- par_postJ4_metr_norm2[i-1,21]
  zeta15_1 <- par_postJ4_metr_norm2[i-1,22]
  zeta16_1 <- par_postJ4_metr_norm2[i-1,23]
  zeta17_1 <- par_postJ4_metr_norm2[i-1,24]
  zeta18_1 <- par_postJ4_metr_norm2[i-1,25]
  zeta19_1 <- par_postJ4_metr_norm2[i-1,26]
  zeta20_1 <- par_postJ4_metr_norm2[i-1,27]
  zeta21_1 <- par_postJ4_metr_norm2[i-1,28]
  zeta22_1 <- par_postJ4_metr_norm2[i-1,29]
  zeta23_1 <- par_postJ4_metr_norm2[i-1,30]
  zeta24_1 <- par_postJ4_metr_norm2[i-1,31]
  zeta25_1 <- par_postJ4_metr_norm2[i-1,32]
  zeta26_1 <- par_postJ4_metr_norm2[i-1,33]
  zeta27_1 <- par_postJ4_metr_norm2[i-1,34]
  zeta28_1 <- par_postJ4_metr_norm2[i-1,35]
  zeta29_1 <- par_postJ4_metr_norm2[i-1,36]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta0_1 * bL + zeta15_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta1_1 * bL + zeta16_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * (zeta2_1 * bL + zeta17_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta3_1 * bL + zeta18_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta4_1 * bL + zeta19_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * (zeta5_1 * bL + zeta20_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta6_1 * bL + zeta21_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta7_1 * bL + zeta22_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG2C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * (zeta8_1 * bL + zeta23_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG3C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * (zeta11_1 * bL + zeta26_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta9_1 * bL + zeta24_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * (zeta10_1 * bL + zeta25_1  * bH) / 
    ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x1 - 77.5))))                                              * ((zeta0_1  + zeta1_1  + zeta2_1 ) * bL + (zeta15_1 + zeta16_1 + zeta17_1)  * bH) + 
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1                             ^ d1) * ((zeta3_1  + zeta4_1  + zeta5_1 ) * bL + (zeta18_1 + zeta19_1 + zeta20_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d1) * ((zeta6_1  + zeta7_1  + zeta8_1 ) * bL + (zeta21_1 + zeta22_1 + zeta23_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d1) * ((zeta9_1  + zeta10_1 + zeta11_1) * bL + (zeta24_1 + zeta25_1 + zeta26_1)  * bH) +
        exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d1) * ((zeta12_1 + zeta13_1 + zeta14_1) * bL + (zeta27_1 + zeta28_1 + zeta29_1)  * bH) )
  
  DS1$PostG4C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1 - DS1$PostG2C2 - DS1$PostG3C0 - DS1$PostG3C1 - DS1$PostG3C2 - DS1$PostG4C0 - DS1$PostG4C1
  
  # - dataset 2
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * (zeta15_1 * c0 + zeta16_1 * c1 + zeta17_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * (zeta18_1 * c0 + zeta19_1 * c1 + zeta20_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG2BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * (zeta21_1 * c0 + zeta22_1 * c1 + zeta23_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG3BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * (zeta24_1 * c0 + zeta25_1 * c1 + zeta26_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
    ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                                     * exp(beta_1 * (x2 - 77.5))))                                              * ((zeta0_1  + zeta15_1) * c0 + (zeta1_1  + zeta16_1) * c1 + (zeta2_1  + zeta17_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1                            * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1                             ^ d2) * ((zeta3_1  + zeta18_1) * c0 + (zeta4_1  + zeta19_1) * c1 + (zeta5_1  + zeta20_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1                   * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1                  ) ^ d2) * ((zeta6_1  + zeta21_1) * c0 + (zeta7_1  + zeta22_1) * c1 + (zeta8_1  + zeta23_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1          * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1         ) ^ d2) * ((zeta9_1  + zeta24_1) * c0 + (zeta10_1 + zeta25_1) * c1 + (zeta11_1 + zeta26_1) * c2) +
        exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * phi3_1 * phi4_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1 * phi3_1 * phi4_1) ^ d2) * ((zeta12_1 + zeta27_1) * c0 + (zeta13_1 + zeta28_1) * c1 + (zeta14_1 + zeta29_1) * c2))
  
  DS2$PostG4BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL - DS2$PostG2BH - DS2$PostG3BL - DS2$PostG3BH - DS2$PostG4BL
  
  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  DS1$G0C0 <- ifelse(g_unif > DS1$PostG0C0, 0, 1)
  DS1$G0C1 <- ifelse((g_unif > DS1$PostG0C0 & g_unif < (DS1$PostG0C0 + DS1$PostG0C1)), 1, 0)
  DS1$G0C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)), 1, 0)
  DS1$G1C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0)), 1, 0)
  DS1$G1C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1)), 1, 0)
  DS1$G1C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2)), 1, 0)
  DS1$G2C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0)), 1, 0)
  DS1$G2C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1)), 1, 0)
  DS1$G2C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2)), 1, 0)
  DS1$G3C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0)), 1, 0)
  DS1$G3C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1)), 1, 0)
  DS1$G3C2 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2)), 1, 0)
  DS1$G4C0 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0)), 1, 0)
  DS1$G4C1 <- ifelse((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0) & g_unif < (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2 + DS1$PostG2C0 + DS1$PostG2C1 + DS1$PostG2C2 + DS1$PostG3C0 + DS1$PostG3C1 + DS1$PostG3C2 + DS1$PostG4C0 + DS1$PostG4C1)), 1, 0)
  DS1$G4C2 <- 1 - DS1$G0C0 - DS1$G0C1 - DS1$G0C2 - DS1$G1C0 - DS1$G1C1 - DS1$G1C2 - DS1$G2C0 - DS1$G2C1 - DS1$G2C2 - DS1$G3C0 - DS1$G3C1 - DS1$G3C2 - DS1$G4C0 - DS1$G4C1
  
  
  g10 <- DS1$G0C0 + DS1$G0C1 + DS1$G0C2
  g11 <- DS1$G1C0 + DS1$G1C1 + DS1$G1C2
  g12 <- DS1$G2C0 + DS1$G2C1 + DS1$G2C2
  g13 <- DS1$G3C0 + DS1$G3C1 + DS1$G3C2
  g14 <- DS1$G4C0 + DS1$G4C1 + DS1$G4C2
  
  c10 <- DS1$G0C0 + DS1$G1C0 + DS1$G2C0 + DS1$G3C0 + DS1$G4C0
  c11 <- DS1$G0C1 + DS1$G1C1 + DS1$G2C1 + DS1$G3C1 + DS1$G4C1
  c12 <- DS1$G0C2 + DS1$G1C2 + DS1$G2C2 + DS1$G3C2 + DS1$G4C2
  
  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  DS2$G0BL <- ifelse(g_unif > DS2$PostG0BL, 0, 1)
  DS2$G0BH <- ifelse((g_unif > DS2$PostG0BL & g_unif < (DS2$PostG0BL + DS2$PostG0BH)), 1, 0)
  DS2$G1BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL)), 1, 0)
  DS2$G1BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH)), 1, 0)
  DS2$G2BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL)), 1, 0)
  DS2$G2BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH)), 1, 0)
  DS2$G3BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL)), 1, 0)
  DS2$G3BH <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH)), 1, 0)
  DS2$G4BL <- ifelse((g_unif > (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH) & g_unif < (DS2$PostG0BL + DS2$PostG0BH + DS2$PostG1BL + DS2$PostG1BH + DS2$PostG2BL + DS2$PostG2BH + DS2$PostG3BL + DS2$PostG3BH + DS2$PostG4BL)), 1, 0)
  DS2$G4BH <- 1 - DS2$G0BL - DS2$G0BH - DS2$G1BL - DS2$G1BH - DS2$G2BL - DS2$G2BH - DS2$G3BL - DS2$G3BH - DS2$G4BL
  
  g20 <- DS2$G0BL + DS2$G0BH
  g21 <- DS2$G1BL + DS2$G1BH
  g22 <- DS2$G2BL + DS2$G2BH
  g23 <- DS2$G3BL + DS2$G3BH
  g24 <- DS2$G4BL + DS2$G4BH
  
  b2L <- DS2$G0BL + DS2$G1BL + DS2$G2BL + DS2$G3BL + DS2$G4BL
  b2H <- DS2$G0BH + DS2$G1BH + DS2$G2BH + DS2$G3BH + DS2$G4BH
  
  
  # - Generate value from Dirichlet
  
  gamma_zeta <- rep(0, 30)
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + (sum(bL * c10 * g10) + sum(b2L * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + (sum(bL * c11 * g10) + sum(b2L * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + (sum(bL * c12 * g10) + sum(b2L * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + (sum(bL * c10 * g11) + sum(b2L * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + (sum(bL * c11 * g11) + sum(b2L * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + (sum(bL * c12 * g11) + sum(b2L * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + (sum(bL * c10 * g12) + sum(b2L * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + (sum(bL * c11 * g12) + sum(b2L * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + (sum(bL * c12 * g12) + sum(b2L * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + (sum(bL * c10 * g13) + sum(b2L * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + (sum(bL * c11 * g13) + sum(b2L * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + (sum(bL * c12 * g13) + sum(b2L * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + (sum(bL * c10 * g14) + sum(b2L * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + (sum(bL * c11 * g14) + sum(b2L * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + (sum(bL * c12 * g14) + sum(b2L * c2 * g24)) / log(nrow(DS_C)), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + (sum(bH * c10 * g10) + sum(b2H * c0 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + (sum(bH * c11 * g10) + sum(b2H * c1 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + (sum(bH * c12 * g10) + sum(b2H * c2 * g20)) / log(nrow(DS_C)), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + (sum(bH * c10 * g11) + sum(b2H * c0 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + (sum(bH * c11 * g11) + sum(b2H * c1 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + (sum(bH * c12 * g11) + sum(b2H * c2 * g21)) / log(nrow(DS_C)), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + (sum(bH * c10 * g12) + sum(b2H * c0 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + (sum(bH * c11 * g12) + sum(b2H * c1 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + (sum(bH * c12 * g12) + sum(b2H * c2 * g22)) / log(nrow(DS_C)), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + (sum(bH * c10 * g13) + sum(b2H * c0 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + (sum(bH * c11 * g13) + sum(b2H * c1 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + (sum(bH * c12 * g13) + sum(b2H * c2 * g23)) / log(nrow(DS_C)), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + (sum(bH * c10 * g14) + sum(b2H * c0 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + (sum(bH * c11 * g14) + sum(b2H * c1 * g24)) / log(nrow(DS_C)), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + (sum(bH * c12 * g14) + sum(b2H * c2 * g24)) / log(nrow(DS_C)), 1)
  
  par_postJ4_metr_norm2[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.00325)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.00325) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.00325))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ4_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.09)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_star) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_1   ) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_star) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_1   ) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.09) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.09))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ4_metr_norm2[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.11)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_star ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_star) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_1    ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_1   ) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_star ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_star) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_1    ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_1   ) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.11) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.11))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)
  
  par_postJ4_metr_norm2[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.0975)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_star ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_star) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_1    ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_1   ) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_star ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_star) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_1    ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_1   ) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi3_star, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.0975) / dtruncnorm(phi3_1, a=(10^(-20)), b=1, mean = phi3_star, sd = 0.0975))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)
  
  par_postJ4_metr_norm2[i, 5] <- phi3
  
  # - phi_4
  
  phi4_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.18)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_star ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_star) + beta_1 * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_1    ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_star ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_star) + beta_1 * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_1    ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_1   ) + beta_1 * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(phi4_star, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.18) / dtruncnorm(phi4_1, a=(10^(-20)), b=1, mean = phi4_star, sd = 0.18))
  
  phi4_unif <- runif(1, 0, 1)
  
  phi4 <- ifelse( phi4_unif > r, phi4_1, phi4_star)
  
  par_postJ4_metr_norm2[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp((sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
              sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
              sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
              sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / log(nrow(DS_C)) ) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_metr_norm2[i, 2] <- beta
  
}

## - Thinned sequence

# - J = 0

par_postJ0_gn1_th <- par_postJ0_gibbs_norm1[seq(10010, 100000, 10),]
#par_postJ0_gn2_th <- par_postJ0_gibbs_norm2[seq(10010, 100000, 10),]
#par_postJ0_mn1_th <- par_postJ0_metr_norm1[seq(10010, 100000, 10),]
#par_postJ0_mn2_th <- par_postJ0_metr_norm2[seq(10010, 100000, 10),]

# - J = 1

par_postJ1_gn1_th <- par_postJ1_gibbs_norm1[seq(10010, 100000, 10),]
#par_postJ1_gn2_th <- par_postJ1_gibbs_norm2[seq(10010, 100000, 10),]
#par_postJ1_mn1_th <- par_postJ1_metr_norm1[seq(10010, 100000, 10),]
#par_postJ1_mn2_th <- par_postJ1_metr_norm2[seq(10010, 100000, 10),]

# - J = 2

par_postJ2_gn1_th <- par_postJ2_gibbs_norm1[seq(10010, 100000, 10),]
#par_postJ2_gn2_th <- par_postJ2_gibbs_norm2[seq(10010, 100000, 10),]
#par_postJ2_mn1_th <- par_postJ2_metr_norm1[seq(10010, 100000, 10),]
#par_postJ2_mn2_th <- par_postJ2_metr_norm2[seq(10010, 100000, 10),]

# - J = 3

par_postJ3_gn1_th <- par_postJ3_gibbs_norm1[seq(10010, 100000, 10),]
#par_postJ3_gn2_th <- par_postJ3_gibbs_norm2[seq(10010, 100000, 10),]
#par_postJ3_mn1_th <- par_postJ3_metr_norm1[seq(10010, 100000, 10),]
#par_postJ3_mn2_th <- par_postJ3_metr_norm2[seq(10010, 100000, 10),]

# - J = 4

par_postJ4_gn1_th <- par_postJ4_gibbs_norm1[seq(10010, 100000, 10),]
#par_postJ4_gn2_th <- par_postJ4_gibbs_norm2[seq(10010, 100000, 10),]
#par_postJ4_mn1_th <- par_postJ4_metr_norm1[seq(10010, 100000, 10),]
#par_postJ4_mn2_th <- par_postJ4_metr_norm2[seq(10010, 100000, 10),]

######################### - WBIC calculation - #####################################

# - J=0
S_sample_J0_gn1 <- par_postJ0_gn1_th[sample(nrow(par_postJ0_gn1_th), 1000),]
#S_sample_J0_gn2 <- par_postJ0_gn2_th[sample(nrow(par_postJ0_gn2_th), 1000),]
#S_sample_J0_mn1 <- par_postJ0_mn1_th[sample(nrow(par_postJ0_mn1_th), 1000),]
#S_sample_J0_mn2 <- par_postJ0_mn2_th[sample(nrow(par_postJ0_mn2_th), 1000),]

# - gn1
avg_log_J0_gn1 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_gn1[,1]
b <- S_sample_J0_gn1[,2]

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_gn1[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

WBIC_J0_gn1 <- -sum(avg_log_J0_gn1)

# - gn2
avg_log_J0_gn2 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_gn2[,1]
b <- S_sample_J0_gn2[,2]

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_gn2[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

WBIC_J0_gn2 <- -sum(avg_log_J0_gn2)

# - mn1
avg_log_J0_mn1 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_mn1[,1]
b <- S_sample_J0_mn1[,2]

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_mn1[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

WBIC_J0_mn1 <- -sum(avg_log_J0_mn1)

# - mn2
avg_log_J0_mn2 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_mn2[,1]
b <- S_sample_J0_mn2[,2]

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_mn2[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

WBIC_J0_mn2 <- -sum(avg_log_J0_mn2)


# - J=1
S_sample_J1_gn1 <- par_postJ1_gn1_th#[sample(nrow(par_postJ1_gn1_th), 1000),]
#S_sample_J1_gn2 <- par_postJ1_gn2_th[sample(nrow(par_postJ1_gn2_th), 1000),]
#S_sample_J1_mn1 <- par_postJ1_mn1_th[sample(nrow(par_postJ1_mn1_th), 1000),]
#S_sample_J1_mn2 <- par_postJ1_mn2_th[sample(nrow(par_postJ1_mn2_th), 1000),]

# - gn1

avg_log_J1_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J1_P2_gn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J1_gn1[,1]
phi1 <- S_sample_J1_gn1[,3]
b <- S_sample_J1_gn1[,2]

z0 <- S_sample_J1_gn1[,4]
z1 <- S_sample_J1_gn1[,5]
z2 <- S_sample_J1_gn1[,6]
z3 <- S_sample_J1_gn1[,7]
z4 <- S_sample_J1_gn1[,8]
z5 <- S_sample_J1_gn1[,9]
z6 <- S_sample_J1_gn1[,10]
z7 <- S_sample_J1_gn1[,11]
z8 <- S_sample_J1_gn1[,12]
z9 <- S_sample_J1_gn1[,13]
z10 <- S_sample_J1_gn1[,14]
z11 <- S_sample_J1_gn1[,15]

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

WBIC_J1_gn1 <- sum(avg_log_J1_P1_gn1) + sum(avg_log_J1_P2_gn1)

# - gn2

avg_log_J1_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J1_P2_gn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J1_gn2[,1]
phi1 <- S_sample_J1_gn2[,3]
b <- S_sample_J1_gn2[,2]

z0 <- S_sample_J1_gn2[,4]
z1 <- S_sample_J1_gn2[,5]
z2 <- S_sample_J1_gn2[,6]
z3 <- S_sample_J1_gn2[,7]
z4 <- S_sample_J1_gn2[,8]
z5 <- S_sample_J1_gn2[,9]
z6 <- S_sample_J1_gn2[,10]
z7 <- S_sample_J1_gn2[,11]
z8 <- S_sample_J1_gn2[,12]
z9 <- S_sample_J1_gn2[,13]
z10 <- S_sample_J1_gn2[,14]
z11 <- S_sample_J1_gn2[,15]

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

WBIC_J1_gn2 <- sum(avg_log_J1_P1_gn2) + sum(avg_log_J1_P2_gn2)

# - mn1

avg_log_J1_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J1_P2_mn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J1_mn1[,1]
phi1 <- S_sample_J1_mn1[,3]
b <- S_sample_J1_mn1[,2]

z0 <- S_sample_J1_mn1[,4]
z1 <- S_sample_J1_mn1[,5]
z2 <- S_sample_J1_mn1[,6]
z3 <- S_sample_J1_mn1[,7]
z4 <- S_sample_J1_mn1[,8]
z5 <- S_sample_J1_mn1[,9]
z6 <- S_sample_J1_mn1[,10]
z7 <- S_sample_J1_mn1[,11]
z8 <- S_sample_J1_mn1[,12]
z9 <- S_sample_J1_mn1[,13]
z10 <- S_sample_J1_mn1[,14]
z11 <- S_sample_J1_mn1[,15]

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

WBIC_J1_mn1 <- sum(avg_log_J1_P1_mn1) + sum(avg_log_J1_P2_mn1)

# - mn2

avg_log_J1_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J1_P2_mn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J1_mn2[,1]
phi1 <- S_sample_J1_mn2[,3]
b <- S_sample_J1_mn2[,2]

z0 <- S_sample_J1_mn2[,4]
z1 <- S_sample_J1_mn2[,5]
z2 <- S_sample_J1_mn2[,6]
z3 <- S_sample_J1_mn2[,7]
z4 <- S_sample_J1_mn2[,8]
z5 <- S_sample_J1_mn2[,9]
z6 <- S_sample_J1_mn2[,10]
z7 <- S_sample_J1_mn2[,11]
z8 <- S_sample_J1_mn2[,12]
z9 <- S_sample_J1_mn2[,13]
z10 <- S_sample_J1_mn2[,14]
z11 <- S_sample_J1_mn2[,15]

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

WBIC_J1_mn2 <- sum(avg_log_J1_P1_mn2) + sum(avg_log_J1_P2_mn2)

# - J=2
S_sample_J2_gn1 <- par_postJ2_gn1_th#[sample(nrow(par_postJ2_gn1_th), 1000),]
#S_sample_J2_gn2 <- par_postJ2_gn2_th[sample(nrow(par_postJ2_gn2_th), 1000),]
#S_sample_J2_mn1 <- par_postJ2_mn1_th[sample(nrow(par_postJ2_mn1_th), 1000),]
#S_sample_J2_mn2 <- par_postJ2_mn2_th[sample(nrow(par_postJ2_mn2_th), 1000),]

# - gn1

avg_log_J2_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J2_P2_gn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J2_gn1[,1]
phi1 <- S_sample_J2_gn1[,3]
phi2 <- S_sample_J2_gn1[,4]
b <- S_sample_J2_gn1[,2]

z0 <- S_sample_J2_gn1[,5]
z1 <- S_sample_J2_gn1[,6]
z2 <- S_sample_J2_gn1[,7]
z3 <- S_sample_J2_gn1[,8]
z4 <- S_sample_J2_gn1[,9]
z5 <- S_sample_J2_gn1[,10]
z6 <- S_sample_J2_gn1[,11]
z7 <- S_sample_J2_gn1[,12]
z8 <- S_sample_J2_gn1[,13]
z9 <- S_sample_J2_gn1[,14]
z10 <- S_sample_J2_gn1[,15]
z11 <- S_sample_J2_gn1[,16]
z12 <- S_sample_J2_gn1[,17]
z13 <- S_sample_J2_gn1[,18]
z14 <- S_sample_J2_gn1[,19]
z15 <- S_sample_J2_gn1[,20]
z16 <- S_sample_J2_gn1[,21]
z17 <- S_sample_J2_gn1[,22]

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

WBIC_J2_gn1 <- sum(avg_log_J2_P1_gn1) + sum(avg_log_J2_P2_gn1)

# - gn2

avg_log_J2_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J2_P2_gn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J2_gn2[,1]
phi1 <- S_sample_J2_gn2[,3]
phi2 <- S_sample_J2_gn2[,4]
b <- S_sample_J2_gn2[,2]

z0 <- S_sample_J2_gn2[,5]
z1 <- S_sample_J2_gn2[,6]
z2 <- S_sample_J2_gn2[,7]
z3 <- S_sample_J2_gn2[,8]
z4 <- S_sample_J2_gn2[,9]
z5 <- S_sample_J2_gn2[,10]
z6 <- S_sample_J2_gn2[,11]
z7 <- S_sample_J2_gn2[,12]
z8 <- S_sample_J2_gn2[,13]
z9 <- S_sample_J2_gn2[,14]
z10 <- S_sample_J2_gn2[,15]
z11 <- S_sample_J2_gn2[,16]
z12 <- S_sample_J2_gn2[,17]
z13 <- S_sample_J2_gn2[,18]
z14 <- S_sample_J2_gn2[,19]
z15 <- S_sample_J2_gn2[,20]
z16 <- S_sample_J2_gn2[,21]
z17 <- S_sample_J2_gn2[,22]

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

WBIC_J2_gn2 <- sum(avg_log_J2_P1_gn2) + sum(avg_log_J2_P2_gn2)

# - mn1

avg_log_J2_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J2_P2_mn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J2_mn1[,1]
phi1 <- S_sample_J2_mn1[,3]
phi2 <- S_sample_J2_mn1[,4]
b <- S_sample_J2_mn1[,2]

z0 <- S_sample_J2_mn1[,5]
z1 <- S_sample_J2_mn1[,6]
z2 <- S_sample_J2_mn1[,7]
z3 <- S_sample_J2_mn1[,8]
z4 <- S_sample_J2_mn1[,9]
z5 <- S_sample_J2_mn1[,10]
z6 <- S_sample_J2_mn1[,11]
z7 <- S_sample_J2_mn1[,12]
z8 <- S_sample_J2_mn1[,13]
z9 <- S_sample_J2_mn1[,14]
z10 <- S_sample_J2_mn1[,15]
z11 <- S_sample_J2_mn1[,16]
z12 <- S_sample_J2_mn1[,17]
z13 <- S_sample_J2_mn1[,18]
z14 <- S_sample_J2_mn1[,19]
z15 <- S_sample_J2_mn1[,20]
z16 <- S_sample_J2_mn1[,21]
z17 <- S_sample_J2_mn1[,22]

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

WBIC_J2_mn1 <- sum(avg_log_J2_P1_mn1) + sum(avg_log_J2_P2_mn1)

# - mn2

avg_log_J2_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J2_P2_mn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J2_mn2[,1]
phi1 <- S_sample_J2_mn2[,3]
phi2 <- S_sample_J2_mn2[,4]
b <- S_sample_J2_mn2[,2]

z0 <- S_sample_J2_mn2[,5]
z1 <- S_sample_J2_mn2[,6]
z2 <- S_sample_J2_mn2[,7]
z3 <- S_sample_J2_mn2[,8]
z4 <- S_sample_J2_mn2[,9]
z5 <- S_sample_J2_mn2[,10]
z6 <- S_sample_J2_mn2[,11]
z7 <- S_sample_J2_mn2[,12]
z8 <- S_sample_J2_mn2[,13]
z9 <- S_sample_J2_mn2[,14]
z10 <- S_sample_J2_mn2[,15]
z11 <- S_sample_J2_mn2[,16]
z12 <- S_sample_J2_mn2[,17]
z13 <- S_sample_J2_mn2[,18]
z14 <- S_sample_J2_mn2[,19]
z15 <- S_sample_J2_mn2[,20]
z16 <- S_sample_J2_mn2[,21]
z17 <- S_sample_J2_mn2[,22]

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

WBIC_J2_mn2 <- sum(avg_log_J2_P1_mn2) + sum(avg_log_J2_P2_mn2)

# - J=3
S_sample_J3_gn1 <- par_postJ3_gn1_th#[sample(nrow(par_postJ3_gn1_th), 1000),]
#S_sample_J3_gn2 <- par_postJ3_gn2_th[sample(nrow(par_postJ3_gn2_th), 1000),]
#S_sample_J3_mn1 <- par_postJ3_mn1_th[sample(nrow(par_postJ3_mn1_th), 1000),]
#S_sample_J3_mn2 <- par_postJ3_mn2_th[sample(nrow(par_postJ3_mn2_th), 1000),]

# - gn1

avg_log_J3_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J3_P2_gn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J3_gn1[,1]
phi1 <- S_sample_J3_gn1[,3]
phi2 <- S_sample_J3_gn1[,4]
phi3 <- S_sample_J3_gn1[,5]
b <- S_sample_J3_gn1[,2]

z0 <- S_sample_J3_gn1[,6]
z1 <- S_sample_J3_gn1[,7]
z2 <- S_sample_J3_gn1[,8]
z3 <- S_sample_J3_gn1[,9]
z4 <- S_sample_J3_gn1[,10]
z5 <- S_sample_J3_gn1[,11]
z6 <- S_sample_J3_gn1[,12]
z7 <- S_sample_J3_gn1[,13]
z8 <- S_sample_J3_gn1[,14]
z9 <- S_sample_J3_gn1[,15]
z10 <- S_sample_J3_gn1[,16]
z11 <- S_sample_J3_gn1[,17]
z12 <- S_sample_J3_gn1[,18]
z13 <- S_sample_J3_gn1[,19]
z14 <- S_sample_J3_gn1[,20]
z15 <- S_sample_J3_gn1[,21]
z16 <- S_sample_J3_gn1[,22]
z17 <- S_sample_J3_gn1[,23]
z18 <- S_sample_J3_gn1[,24]
z19 <- S_sample_J3_gn1[,25]
z20 <- S_sample_J3_gn1[,26]
z21 <- S_sample_J3_gn1[,27]
z22 <- S_sample_J3_gn1[,28]
z23 <- S_sample_J3_gn1[,29]

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

WBIC_J3_gn1 <- sum(avg_log_J3_P1_gn1) + sum(avg_log_J3_P2_gn1)

# - gn2

avg_log_J3_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J3_P2_gn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J3_gn2[,1]
phi1 <- S_sample_J3_gn2[,3]
phi2 <- S_sample_J3_gn2[,4]
phi3 <- S_sample_J3_gn2[,5]
b <- S_sample_J3_gn2[,2]

z0 <- S_sample_J3_gn2[,6]
z1 <- S_sample_J3_gn2[,7]
z2 <- S_sample_J3_gn2[,8]
z3 <- S_sample_J3_gn2[,9]
z4 <- S_sample_J3_gn2[,10]
z5 <- S_sample_J3_gn2[,11]
z6 <- S_sample_J3_gn2[,12]
z7 <- S_sample_J3_gn2[,13]
z8 <- S_sample_J3_gn2[,14]
z9 <- S_sample_J3_gn2[,15]
z10 <- S_sample_J3_gn2[,16]
z11 <- S_sample_J3_gn2[,17]
z12 <- S_sample_J3_gn2[,18]
z13 <- S_sample_J3_gn2[,19]
z14 <- S_sample_J3_gn2[,20]
z15 <- S_sample_J3_gn2[,21]
z16 <- S_sample_J3_gn2[,22]
z17 <- S_sample_J3_gn2[,23]
z18 <- S_sample_J3_gn2[,24]
z19 <- S_sample_J3_gn2[,25]
z20 <- S_sample_J3_gn2[,26]
z21 <- S_sample_J3_gn2[,27]
z22 <- S_sample_J3_gn2[,28]
z23 <- S_sample_J3_gn2[,29]

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

WBIC_J3_gn2 <- sum(avg_log_J3_P1_gn2) + sum(avg_log_J3_P2_gn2)

# - mn1

avg_log_J3_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J3_P2_mn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J3_mn1[,1]
phi1 <- S_sample_J3_mn1[,3]
phi2 <- S_sample_J3_mn1[,4]
phi3 <- S_sample_J3_mn1[,5]
b <- S_sample_J3_mn1[,2]

z0 <- S_sample_J3_mn1[,6]
z1 <- S_sample_J3_mn1[,7]
z2 <- S_sample_J3_mn1[,8]
z3 <- S_sample_J3_mn1[,9]
z4 <- S_sample_J3_mn1[,10]
z5 <- S_sample_J3_mn1[,11]
z6 <- S_sample_J3_mn1[,12]
z7 <- S_sample_J3_mn1[,13]
z8 <- S_sample_J3_mn1[,14]
z9 <- S_sample_J3_mn1[,15]
z10 <- S_sample_J3_mn1[,16]
z11 <- S_sample_J3_mn1[,17]
z12 <- S_sample_J3_mn1[,18]
z13 <- S_sample_J3_mn1[,19]
z14 <- S_sample_J3_mn1[,20]
z15 <- S_sample_J3_mn1[,21]
z16 <- S_sample_J3_mn1[,22]
z17 <- S_sample_J3_mn1[,23]
z18 <- S_sample_J3_mn1[,24]
z19 <- S_sample_J3_mn1[,25]
z20 <- S_sample_J3_mn1[,26]
z21 <- S_sample_J3_mn1[,27]
z22 <- S_sample_J3_mn1[,28]
z23 <- S_sample_J3_mn1[,29]

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

WBIC_J3_mn1 <- sum(avg_log_J3_P1_mn1) + sum(avg_log_J3_P2_mn1)

# - mn2

avg_log_J3_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J3_P2_mn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J3_mn2[,1]
phi1 <- S_sample_J3_mn2[,3]
phi2 <- S_sample_J3_mn2[,4]
phi3 <- S_sample_J3_mn2[,5]
b <- S_sample_J3_mn2[,2]

z0 <- S_sample_J3_mn2[,6]
z1 <- S_sample_J3_mn2[,7]
z2 <- S_sample_J3_mn2[,8]
z3 <- S_sample_J3_mn2[,9]
z4 <- S_sample_J3_mn2[,10]
z5 <- S_sample_J3_mn2[,11]
z6 <- S_sample_J3_mn2[,12]
z7 <- S_sample_J3_mn2[,13]
z8 <- S_sample_J3_mn2[,14]
z9 <- S_sample_J3_mn2[,15]
z10 <- S_sample_J3_mn2[,16]
z11 <- S_sample_J3_mn2[,17]
z12 <- S_sample_J3_mn2[,18]
z13 <- S_sample_J3_mn2[,19]
z14 <- S_sample_J3_mn2[,20]
z15 <- S_sample_J3_mn2[,21]
z16 <- S_sample_J3_mn2[,22]
z17 <- S_sample_J3_mn2[,23]
z18 <- S_sample_J3_mn2[,24]
z19 <- S_sample_J3_mn2[,25]
z20 <- S_sample_J3_mn2[,26]
z21 <- S_sample_J3_mn2[,27]
z22 <- S_sample_J3_mn2[,28]
z23 <- S_sample_J3_mn2[,29]

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

WBIC_J3_mn2 <- sum(avg_log_J3_P1_mn2) + sum(avg_log_J3_P2_mn2)

# - J=4
S_sample_J4_gn1 <- par_postJ4_gn1_th#[sample(nrow(par_postJ4_gn1_th), 1000),]
#S_sample_J4_gn2 <- par_postJ4_gn2_th[sample(nrow(par_postJ4_gn2_th), 1000),]
#S_sample_J4_mn1 <- par_postJ4_mn1_th[sample(nrow(par_postJ4_mn1_th), 1000),]
#S_sample_J4_mn2 <- par_postJ4_mn2_th[sample(nrow(par_postJ4_mn2_th), 1000),]

# - gn1

avg_log_J4_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J4_P2_gn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J4_gn1[,1]
phi1 <- S_sample_J4_gn1[,3]
phi2 <- S_sample_J4_gn1[,4]
phi3 <- S_sample_J4_gn1[,5]
phi4 <- S_sample_J4_gn1[,6]
b <- S_sample_J4_gn1[,2]

z0 <- S_sample_J4_gn1[,7]
z1 <- S_sample_J4_gn1[,8]
z2 <- S_sample_J4_gn1[,9]
z3 <- S_sample_J4_gn1[,10]
z4 <- S_sample_J4_gn1[,11]
z5 <- S_sample_J4_gn1[,12]
z6 <- S_sample_J4_gn1[,13]
z7 <- S_sample_J4_gn1[,14]
z8 <- S_sample_J4_gn1[,15]
z9 <- S_sample_J4_gn1[,16]

z10 <- S_sample_J4_gn1[,17]
z11 <- S_sample_J4_gn1[,18]
z12 <- S_sample_J4_gn1[,19]
z13 <- S_sample_J4_gn1[,20]
z14 <- S_sample_J4_gn1[,21]
z15 <- S_sample_J4_gn1[,22]
z16 <- S_sample_J4_gn1[,23]
z17 <- S_sample_J4_gn1[,24]
z18 <- S_sample_J4_gn1[,25]
z19 <- S_sample_J4_gn1[,26]
z20 <- S_sample_J4_gn1[,27]
z21 <- S_sample_J4_gn1[,28]
z22 <- S_sample_J4_gn1[,29]
z23 <- S_sample_J4_gn1[,30]
z24 <- S_sample_J4_gn1[,31]
z25 <- S_sample_J4_gn1[,32]
z26 <- S_sample_J4_gn1[,33]
z27 <- S_sample_J4_gn1[,34]
z28 <- S_sample_J4_gn1[,35]
z29 <- S_sample_J4_gn1[,36]

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}


for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

WBIC_J4_gn1 <- sum(avg_log_J4_P1_gn1) + sum(avg_log_J4_P2_gn1)

# - gn2

avg_log_J4_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J4_P2_gn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J4_gn2[,1]
phi1 <- S_sample_J4_gn2[,3]
phi2 <- S_sample_J4_gn2[,4]
phi3 <- S_sample_J4_gn2[,5]
phi4 <- S_sample_J4_gn2[,6]
b <- S_sample_J4_gn2[,2]

z0 <- S_sample_J4_gn2[,7]
z1 <- S_sample_J4_gn2[,8]
z2 <- S_sample_J4_gn2[,9]
z3 <- S_sample_J4_gn2[,10]
z4 <- S_sample_J4_gn2[,11]
z5 <- S_sample_J4_gn2[,12]
z6 <- S_sample_J4_gn2[,13]
z7 <- S_sample_J4_gn2[,14]
z8 <- S_sample_J4_gn2[,15]
z9 <- S_sample_J4_gn2[,16]

z10 <- S_sample_J4_gn2[,17]
z11 <- S_sample_J4_gn2[,18]
z12 <- S_sample_J4_gn2[,19]
z13 <- S_sample_J4_gn2[,20]
z14 <- S_sample_J4_gn2[,21]
z15 <- S_sample_J4_gn2[,22]
z16 <- S_sample_J4_gn2[,23]
z17 <- S_sample_J4_gn2[,24]
z18 <- S_sample_J4_gn2[,25]
z19 <- S_sample_J4_gn2[,26]
z20 <- S_sample_J4_gn2[,27]
z21 <- S_sample_J4_gn2[,28]
z22 <- S_sample_J4_gn2[,29]
z23 <- S_sample_J4_gn2[,30]
z24 <- S_sample_J4_gn2[,31]
z25 <- S_sample_J4_gn2[,32]
z26 <- S_sample_J4_gn2[,33]
z27 <- S_sample_J4_gn2[,34]
z28 <- S_sample_J4_gn2[,35]
z29 <- S_sample_J4_gn2[,36]

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}


for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

WBIC_J4_gn2 <- sum(avg_log_J4_P1_gn2) + sum(avg_log_J4_P2_gn2)

# - mn1

avg_log_J4_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J4_P2_mn1 <- rep(0, nrow(DS2))

phi0 <- S_sample_J4_mn1[,1]
phi1 <- S_sample_J4_mn1[,3]
phi2 <- S_sample_J4_mn1[,4]
phi3 <- S_sample_J4_mn1[,5]
phi4 <- S_sample_J4_mn1[,6]
b <- S_sample_J4_mn1[,2]

z0 <- S_sample_J4_mn1[,7]
z1 <- S_sample_J4_mn1[,8]
z2 <- S_sample_J4_mn1[,9]
z3 <- S_sample_J4_mn1[,10]
z4 <- S_sample_J4_mn1[,11]
z5 <- S_sample_J4_mn1[,12]
z6 <- S_sample_J4_mn1[,13]
z7 <- S_sample_J4_mn1[,14]
z8 <- S_sample_J4_mn1[,15]
z9 <- S_sample_J4_mn1[,16]

z10 <- S_sample_J4_mn1[,17]
z11 <- S_sample_J4_mn1[,18]
z12 <- S_sample_J4_mn1[,19]
z13 <- S_sample_J4_mn1[,20]
z14 <- S_sample_J4_mn1[,21]
z15 <- S_sample_J4_mn1[,22]
z16 <- S_sample_J4_mn1[,23]
z17 <- S_sample_J4_mn1[,24]
z18 <- S_sample_J4_mn1[,25]
z19 <- S_sample_J4_mn1[,26]
z20 <- S_sample_J4_mn1[,27]
z21 <- S_sample_J4_mn1[,28]
z22 <- S_sample_J4_mn1[,29]
z23 <- S_sample_J4_mn1[,30]
z24 <- S_sample_J4_mn1[,31]
z25 <- S_sample_J4_mn1[,32]
z26 <- S_sample_J4_mn1[,33]
z27 <- S_sample_J4_mn1[,34]
z28 <- S_sample_J4_mn1[,35]
z29 <- S_sample_J4_mn1[,36]

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}


for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

WBIC_J4_mn1 <- sum(avg_log_J4_P1_mn1) + sum(avg_log_J4_P2_mn1)

## - mn2

avg_log_J4_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J4_P2_mn2 <- rep(0, nrow(DS2))

phi0 <- S_sample_J4_mn2[,1]
phi1 <- S_sample_J4_mn2[,3]
phi2 <- S_sample_J4_mn2[,4]
phi3 <- S_sample_J4_mn2[,5]
phi4 <- S_sample_J4_mn2[,6]
b <- S_sample_J4_mn2[,2]

z0 <- S_sample_J4_mn2[,7]
z1 <- S_sample_J4_mn2[,8]
z2 <- S_sample_J4_mn2[,9]
z3 <- S_sample_J4_mn2[,10]
z4 <- S_sample_J4_mn2[,11]
z5 <- S_sample_J4_mn2[,12]
z6 <- S_sample_J4_mn2[,13]
z7 <- S_sample_J4_mn2[,14]
z8 <- S_sample_J4_mn2[,15]
z9 <- S_sample_J4_mn2[,16]

z10 <- S_sample_J4_mn2[,17]
z11 <- S_sample_J4_mn2[,18]
z12 <- S_sample_J4_mn2[,19]
z13 <- S_sample_J4_mn2[,20]
z14 <- S_sample_J4_mn2[,21]
z15 <- S_sample_J4_mn2[,22]
z16 <- S_sample_J4_mn2[,23]
z17 <- S_sample_J4_mn2[,24]
z18 <- S_sample_J4_mn2[,25]
z19 <- S_sample_J4_mn2[,26]
z20 <- S_sample_J4_mn2[,27]
z21 <- S_sample_J4_mn2[,28]
z22 <- S_sample_J4_mn2[,29]
z23 <- S_sample_J4_mn2[,30]
z24 <- S_sample_J4_mn2[,31]
z25 <- S_sample_J4_mn2[,32]
z26 <- S_sample_J4_mn2[,33]
z27 <- S_sample_J4_mn2[,34]
z28 <- S_sample_J4_mn2[,35]
z29 <- S_sample_J4_mn2[,36]

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}


for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

WBIC_J4_mn2 <- sum(avg_log_J4_P1_mn2) + sum(avg_log_J4_P2_mn2)




