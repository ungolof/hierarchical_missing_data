#################### - Acceptance/Rejection - ###############################

# - J=0
ar_beta_J0_C11 <- 0
beta_J0_C11 <- par_postJ0_gibbs_norm1[,2]
ar_beta_J0_C12 <- 0
beta_J0_C12 <- par_postJ0_gibbs_norm2[,2]

for(i in 90001:100000){
  ar_beta_J0_C11   <- ifelse(  beta_J0_C11[i]==  beta_J0_C11[i-1], ar_beta_J0_C11, ar_beta_J0_C11 + 1) 
  
  #ar_beta_J0_C12   <- ifelse(  beta_J0_C12[i]==  beta_J0_C12[i-1], ar_beta_J0_C12, ar_beta_J0_C12 + 1) 
  
}


ar_alpha_J0_C21 <- 0
alpha_J0_C21 <- par_postJ0_metr_norm1[,1]
ar_alpha_J0_C22 <- 0
alpha_J0_C22 <- par_postJ0_metr_norm2[,1]
ar_beta_J0_C21 <- 0
beta_J0_C21 <- par_postJ0_metr_norm1[,2]
ar_beta_J0_C22 <- 0
beta_J0_C22 <- par_postJ0_metr_norm2[,2]

for(i in 90001:100000){
  ar_alpha_J0_C21  <- ifelse( alpha_J0_C21[i]== alpha_J0_C21[i-1], ar_alpha_J0_C21, ar_alpha_J0_C21 + 1) 
  ar_beta_J0_C21   <- ifelse(  beta_J0_C21[i]==  beta_J0_C21[i-1], ar_beta_J0_C21, ar_beta_J0_C21 + 1) 
  
  ar_alpha_J0_C22  <- ifelse( alpha_J0_C22[i]== alpha_J0_C22[i-1], ar_alpha_J0_C22, ar_alpha_J0_C22 + 1) 
  ar_beta_J0_C22   <- ifelse(  beta_J0_C22[i]==  beta_J0_C22[i-1], ar_beta_J0_C22, ar_beta_J0_C22 + 1) 
  
}
ar_alpha_J0_C21
ar_beta_J0_C21

# - J=1
ar_beta_J1_C11 <- 0
beta_J1_C11 <- par_postJ1_gibbs_norm1[,2]

ar_beta_J1_C12 <- 0
beta_J1_C12 <- par_postJ1_gibbs_norm2[,2]

for(i in 90001:100000){
  ar_beta_J1_C11   <- ifelse(  beta_J1_C11[i]==  beta_J1_C11[i-1], ar_beta_J1_C11, ar_beta_J1_C11 + 1) 

  #ar_beta_J1_C12   <- ifelse(  beta_J1_C12[i]==  beta_J1_C12[i-1], ar_beta_J1_C12, ar_beta_J1_C12 + 1) 
}

ar_alpha_J1_C21 <- 0
alpha_J1_C21 <- par_postJ1_metr_norm1[,1]
ar_beta_J1_C21 <- 0
beta_J1_C21 <- par_postJ1_metr_norm1[,2]
ar_gamma1_J1_C21 <- 0
gamma1_J1_C21 <- par_postJ1_metr_norm1[,3]

ar_alpha_J1_C22 <- 0
alpha_J1_C22 <- par_postJ1_metr_norm2[,1]
ar_beta_J1_C22 <- 0
beta_J1_C22 <- par_postJ1_metr_norm2[,2]
ar_gamma1_J1_C22 <- 0
gamma1_J1_C22 <- par_postJ1_metr_norm2[,3]

for(i in 9001:10000){
  ar_alpha_J1_C21  <- ifelse( alpha_J1_C21[i]== alpha_J1_C21[i-1], ar_alpha_J1_C21, ar_alpha_J1_C21 + 1) 
  ar_beta_J1_C21   <- ifelse(  beta_J1_C21[i]==  beta_J1_C21[i-1], ar_beta_J1_C21, ar_beta_J1_C21 + 1) 
  ar_gamma1_J1_C21 <- ifelse(gamma1_J1_C21[i]==gamma1_J1_C21[i-1], ar_gamma1_J1_C21, ar_gamma1_J1_C21 + 1) 
  
  ar_alpha_J1_C22  <- ifelse( alpha_J1_C22[i]== alpha_J1_C22[i-1], ar_alpha_J1_C22, ar_alpha_J1_C22 + 1) 
  ar_beta_J1_C22   <- ifelse(  beta_J1_C22[i]==  beta_J1_C22[i-1], ar_beta_J1_C22, ar_beta_J1_C22 + 1) 
  ar_gamma1_J1_C22 <- ifelse(gamma1_J1_C22[i]==gamma1_J1_C22[i-1], ar_gamma1_J1_C22, ar_gamma1_J1_C22 + 1) 
}

ar_alpha_J1_C21
ar_beta_J1_C21
ar_gamma1_J1_C21

# - J=2
ar_beta_J2_C11 <- 0
beta_J2_C11 <- par_postJ2_gibbs_norm1[,2]

ar_beta_J2_C12 <- 0
beta_J2_C12 <- par_postJ2_gibbs_norm2[,2]

for(i in 2:100000){
  ar_beta_J2_C11   <- ifelse(  beta_J2_C11[i]==  beta_J2_C11[i-1], ar_beta_J2_C11, ar_beta_J2_C11 + 1) 

  ar_beta_J2_C12   <- ifelse(  beta_J2_C12[i]==  beta_J2_C12[i-1], ar_beta_J2_C12, ar_beta_J2_C12 + 1) 
}


ar_alpha_J2_C21 <- 0
alpha_J2_C21 <- par_postJ2_metr_norm1[,1]
ar_beta_J2_C21 <- 0
beta_J2_C21 <- par_postJ2_metr_norm1[,2]
ar_gamma1_J2_C21 <- 0
gamma1_J2_C21 <- par_postJ2_metr_norm1[,3]
ar_gamma2_J2_C21 <- 0
gamma2_J2_C21 <- par_postJ2_metr_norm1[,4]


ar_alpha_J2_C22 <- 0
alpha_J2_C22 <- par_postJ2_metr_norm2[,1]
ar_beta_J2_C22 <- 0
beta_J2_C22 <- par_postJ2_metr_norm2[,2]
ar_gamma1_J2_C22 <- 0
gamma1_J2_C22 <- par_postJ2_metr_norm2[,3]
ar_gamma2_J2_C22 <- 0
gamma2_J2_C22 <- par_postJ2_metr_norm2[,4]

for(i in 2:20000){
  ar_alpha_J2_C21  <- ifelse( alpha_J2_C21[i]== alpha_J2_C21[i-1], ar_alpha_J2_C21, ar_alpha_J2_C21 + 1) 
  ar_beta_J2_C21   <- ifelse(  beta_J2_C21[i]==  beta_J2_C21[i-1], ar_beta_J2_C21, ar_beta_J2_C21 + 1) 

  ar_gamma2_J2_C21 <- ifelse(gamma2_J2_C21[i]==gamma2_J2_C21[i-1], ar_gamma2_J2_C21, ar_gamma2_J2_C21 + 1)   
  ar_gamma1_J2_C21 <- ifelse(gamma1_J2_C21[i]==gamma1_J2_C21[i-1], ar_gamma1_J2_C21, ar_gamma1_J2_C21 + 1) 

  ar_alpha_J2_C22  <- ifelse( alpha_J2_C22[i]== alpha_J2_C22[i-1], ar_alpha_J2_C22, ar_alpha_J2_C22 + 1) 
  ar_beta_J2_C22   <- ifelse(  beta_J2_C22[i]==  beta_J2_C22[i-1], ar_beta_J2_C22, ar_beta_J2_C22 + 1) 
  ar_gamma2_J2_C22 <- ifelse(gamma2_J2_C22[i]==gamma2_J2_C22[i-1], ar_gamma2_J2_C22, ar_gamma2_J2_C22 + 1)   
  ar_gamma1_J2_C22 <- ifelse(gamma1_J2_C22[i]==gamma1_J2_C22[i-1], ar_gamma1_J2_C22, ar_gamma1_J2_C22 + 1) 

  }

# - J=3
ar_beta_J3_C11 <- 0
beta_J3_C11 <- par_postJ3_gibbs_norm1[,2]

ar_beta_J3_C12 <- 0
beta_J3_C12 <- par_postJ3_gibbs_norm2[,2]


for(i in 2:100000){
  ar_beta_J3_C11   <- ifelse(  beta_J3_C11[i]==  beta_J3_C11[i-1], ar_beta_J3_C11, ar_beta_J3_C11 + 1) 
 
  ar_beta_J3_C12   <- ifelse(  beta_J3_C12[i]==  beta_J3_C12[i-1], ar_beta_J3_C12, ar_beta_J3_C12 + 1) 
 
}

ar_alpha_J3_C21 <- 0
alpha_J3_C21 <- par_postJ3_metr_norm1[,1]
ar_beta_J3_C21 <- 0
beta_J3_C21 <- par_postJ3_metr_norm1[,2]
ar_gamma1_J3_C21 <- 0
gamma1_J3_C21 <- par_postJ3_metr_norm1[,3]
ar_gamma2_J3_C21 <- 0
gamma2_J3_C21 <- par_postJ3_metr_norm1[,4]
ar_gamma3_J3_C21 <- 0
gamma3_J3_C21 <- par_postJ3_metr_norm1[,5]

ar_alpha_J3_C22 <- 0
alpha_J3_C22 <- par_postJ3_metr_norm2[,1]
ar_beta_J3_C22 <- 0
beta_J3_C22 <- par_postJ3_metr_norm2[,2]
ar_gamma1_J3_C22 <- 0
gamma1_J3_C22 <- par_postJ3_metr_norm2[,3]
ar_gamma2_J3_C22 <- 0
gamma2_J3_C22 <- par_postJ3_metr_norm2[,4]
ar_gamma3_J3_C22 <- 0
gamma3_J3_C22 <- par_postJ3_metr_norm2[,5]

for(i in 2:100000){
  ar_alpha_J3_C21  <- ifelse( alpha_J3_C21[i]== alpha_J3_C21[i-1],  ar_alpha_J3_C21,  ar_alpha_J3_C21 + 1) 
  ar_beta_J3_C21   <- ifelse(  beta_J3_C21[i]==  beta_J3_C21[i-1],   ar_beta_J3_C21,   ar_beta_J3_C21 + 1) 
  ar_gamma1_J3_C21 <- ifelse(gamma1_J3_C21[i]==gamma1_J3_C21[i-1], ar_gamma1_J3_C21, ar_gamma1_J3_C21 + 1) 
  ar_gamma2_J3_C21 <- ifelse(gamma2_J3_C21[i]==gamma2_J3_C21[i-1], ar_gamma2_J3_C21, ar_gamma2_J3_C21 + 1)   
  ar_gamma3_J3_C21 <- ifelse(gamma3_J3_C21[i]==gamma3_J3_C21[i-1], ar_gamma3_J3_C21, ar_gamma3_J3_C21 + 1)   

  ar_alpha_J3_C22  <- ifelse( alpha_J3_C22[i]== alpha_J3_C22[i-1], ar_alpha_J3_C22, ar_alpha_J3_C22 + 1) 
  ar_beta_J3_C22   <- ifelse(  beta_J3_C22[i]==  beta_J3_C22[i-1], ar_beta_J3_C22, ar_beta_J3_C22 + 1) 
  ar_gamma1_J3_C22 <- ifelse(gamma1_J3_C22[i]==gamma1_J3_C22[i-1], ar_gamma1_J3_C22, ar_gamma1_J3_C22 + 1) 
  ar_gamma2_J3_C22 <- ifelse(gamma2_J3_C22[i]==gamma2_J3_C22[i-1], ar_gamma2_J3_C22, ar_gamma2_J3_C22 + 1)   
  ar_gamma3_J3_C22 <- ifelse(gamma3_J3_C22[i]==gamma3_J3_C22[i-1], ar_gamma3_J3_C22, ar_gamma3_J3_C22 + 1)   

}

# - J=4
ar_beta_J4_C11 <- 0
beta_J4_C11 <- par_postJ4_gibbs_norm1[,2]

ar_beta_J4_C12 <- 0
beta_J4_C12 <- par_postJ4_gibbs_norm2[,2]


for(i in 2:100000){
  ar_beta_J4_C11   <- ifelse(  beta_J4_C11[i]==  beta_J4_C11[i-1], ar_beta_J4_C11, ar_beta_J4_C11 + 1) 
  
  ar_beta_J4_C12   <- ifelse(  beta_J4_C12[i]==  beta_J4_C12[i-1], ar_beta_J4_C12, ar_beta_J4_C12 + 1) 
  
}


ar_alpha_J4_C21 <- 0
alpha_J4_C21 <- par_postJ4_metr_norm1[,1]
ar_beta_J4_C21 <- 0
beta_J4_C21 <- par_postJ4_metr_norm1[,2]
ar_gamma1_J4_C21 <- 0
gamma1_J4_C21 <- par_postJ4_metr_norm1[,3]
ar_gamma2_J4_C21 <- 0
gamma2_J4_C21 <- par_postJ4_metr_norm1[,4]
ar_gamma3_J4_C21 <- 0
gamma3_J4_C21 <- par_postJ4_metr_norm1[,5]
ar_gamma4_J4_C21 <- 0
gamma4_J4_C21 <- par_postJ4_metr_norm1[,6]

ar_alpha_J4_C22 <- 0
alpha_J4_C22 <- par_postJ4_metr_norm2[,1]
ar_beta_J4_C22 <- 0
beta_J4_C22 <- par_postJ4_metr_norm2[,2]
ar_gamma1_J4_C22 <- 0
gamma1_J4_C22 <- par_postJ4_metr_norm2[,3]
ar_gamma2_J4_C22 <- 0
gamma2_J4_C22 <- par_postJ4_metr_norm2[,4]
ar_gamma3_J4_C22 <- 0
gamma3_J4_C22 <- par_postJ4_metr_norm2[,5]
ar_gamma4_J4_C22 <- 0
gamma4_J4_C22 <- par_postJ4_metr_norm2[,6]

for(i in 2:100000){
  ar_alpha_J4_C21  <- ifelse( alpha_J4_C21[i]== alpha_J4_C21[i-1],  ar_alpha_J4_C21,  ar_alpha_J4_C21 + 1) 
  ar_beta_J4_C21   <- ifelse(  beta_J4_C21[i]==  beta_J4_C21[i-1],   ar_beta_J4_C21,   ar_beta_J4_C21 + 1) 
  ar_gamma1_J4_C21 <- ifelse(gamma1_J4_C21[i]==gamma1_J4_C21[i-1], ar_gamma1_J4_C21, ar_gamma1_J4_C21 + 1) 
  ar_gamma2_J4_C21 <- ifelse(gamma2_J4_C21[i]==gamma2_J4_C21[i-1], ar_gamma2_J4_C21, ar_gamma2_J4_C21 + 1)   
  ar_gamma3_J4_C21 <- ifelse(gamma3_J4_C21[i]==gamma3_J4_C21[i-1], ar_gamma3_J4_C21, ar_gamma3_J4_C21 + 1)   
  ar_gamma4_J4_C21 <- ifelse(gamma4_J4_C21[i]==gamma4_J4_C21[i-1], ar_gamma4_J4_C21, ar_gamma3_J4_C21 + 1)   
  
  ar_alpha_J4_C22  <- ifelse( alpha_J4_C22[i]== alpha_J4_C22[i-1], ar_alpha_J4_C22, ar_alpha_J4_C22 + 1) 
  ar_beta_J4_C22   <- ifelse(  beta_J4_C22[i]==  beta_J4_C22[i-1], ar_beta_J4_C22, ar_beta_J4_C22 + 1) 
  ar_gamma1_J4_C22 <- ifelse(gamma1_J4_C22[i]==gamma1_J4_C22[i-1], ar_gamma1_J4_C22, ar_gamma1_J4_C22 + 1) 
  ar_gamma2_J4_C22 <- ifelse(gamma2_J4_C22[i]==gamma2_J4_C22[i-1], ar_gamma2_J4_C22, ar_gamma2_J4_C22 + 1)   
  ar_gamma3_J4_C22 <- ifelse(gamma3_J4_C22[i]==gamma3_J4_C22[i-1], ar_gamma3_J4_C22, ar_gamma3_J4_C22 + 1)   
  ar_gamma4_J4_C22 <- ifelse(gamma4_J4_C22[i]==gamma4_J4_C22[i-1], ar_gamma4_J4_C22, ar_gamma4_J4_C22 + 1)   
}

ar_alpha_J4_C21
ar_alpha_J4_C22
ar_beta_J4_C21
ar_beta_J4_C22
ar_gamma1_J4_C21
ar_gamma1_J4_C22
ar_gamma2_J4_C21
ar_gamma2_J4_C22
ar_gamma3_J4_C21
ar_gamma3_J4_C22
ar_gamma4_J4_C21
ar_gamma4_J4_C22

################################# - Traceplots, thinning and posterior densities - ######################################################

# - J = 0


## - Thinned sequence

par_postJ0_gn1_th <- par_postJ0_gibbs_norm1[seq(10010, 100000, 10),]
par_postJ0_gn2_th <- par_postJ0_gibbs_norm2[seq(10010, 100000, 10),]
par_postJ0_mn1_th <- par_postJ0_metr_norm1[seq(10010, 100000, 10),]
par_postJ0_mn2_th <- par_postJ0_metr_norm2[seq(10010, 100000, 10),]

## - Density
alpha_J0_gn1 <- density(log(par_postJ0_gn1_th[,1]))
alpha_J0_gn2 <- density(log(par_postJ0_gn2_th[,1]))
alpha_J0_mn1 <- density(log(par_postJ0_mn1_th[,1]))
alpha_J0_mn2 <- density(log(par_postJ0_mn2_th[,1]))

beta_J0_gn1 <- density(par_postJ0_gn1_th[,2])
beta_J0_gn2 <- density(par_postJ0_gn2_th[,2])
beta_J0_mn1 <- density(par_postJ0_mn1_th[,2])
beta_J0_mn2 <- density(par_postJ0_mn2_th[,2])

## - Traceplot

par(mfrow=c(2,3), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, log(par_postJ0_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha",ylim=c(-3.5, -2.5)) #ylim=c(-3.2, -1), 
lines(1:100000, log(par_postJ0_gibbs_norm2[,1]), type = "l", col="red")
lines(1:100000, log(par_postJ0_metr_norm1[,1]), type = "l", col="blue")
lines(1:100000, log(par_postJ0_metr_norm2[,1]), type = "l", col="green")

plot(c(1:100000), log(par_postJ0_gn1_erg_m[,1]), type="l", ylab = "alpha", xlab="Iteration", ylim=c(-3.5,-2.5)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ0_gn2_erg_m[,1]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ0_mn1_erg_m[,1]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ0_mn2_erg_m[,1]), type="l", col="green")   #, lty=4

plot(alpha_J0_gn1, xlab = "alpha", ylim=c(0, 30), main=" ") #, main="Density of the posterior distribution of alpha"
lines(alpha_J0_gn2, col="red")
lines(alpha_J0_mn1, col="blue")
lines(alpha_J0_mn2, col="green")
abline(v=log(mean(par_postJ0_gn1_th[,1])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ0_gn2_th[,1])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ0_mn1_th[,1])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ0_mn2_th[,1])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ0_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim=c(0.08, 0.15)) #
lines(1:100000, par_postJ0_gibbs_norm2[,2], type = "l", col="red")
lines(1:100000, par_postJ0_metr_norm1[,2], type = "l", col="blue")
lines(1:100000, par_postJ0_metr_norm2[,2], type = "l", col="green")

plot(c(1:100000), par_postJ0_gn1_erg_m[,2], type="l", ylab = "beta", xlab="Iteration", ylim=c(0.08, 0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ0_gn2_erg_m[,2], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ0_mn1_erg_m[,2], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ0_mn2_erg_m[,2], type="l", col="green")   #, lty=4

plot(beta_J0_gn1, xlab = "beta", ylim=c(0, 250), main=" ") #, ylim=c(0, 10.3), 
lines(beta_J0_gn2, col="red")
lines(beta_J0_mn1, col="blue")
lines(beta_J0_mn2, col="green")
abline(v=mean(par_postJ0_gn1_th[,2]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ0_gn2_th[,2]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ0_mn1_th[,2]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ0_mn2_th[,2]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

# - J = 1 

## - Thinned sequence

par_postJ1_gn1_th <- par_postJ1_gibbs_norm1[seq(10010, 100000, 10),]
par_postJ1_gn2_th <- par_postJ1_gibbs_norm2[seq(10010, 100000, 10),]
par_postJ1_mn1_th <- par_postJ1_metr_norm1[seq(10010, 100000, 10),]
par_postJ1_mn2_th <- par_postJ1_metr_norm2[seq(10010, 100000, 10),]

## - Density
### - tau
alpha_J1_gn1 <- density(log(par_postJ1_gn1_th[,1]))
alpha_J1_gn2 <- density(log(par_postJ1_gn2_th[,1]))
alpha_J1_mn1 <- density(log(par_postJ1_mn1_th[,1]))
alpha_J1_mn2 <- density(log(par_postJ1_mn2_th[,1]))

beta_J1_gn1 <- density(par_postJ1_gn1_th[,2])
beta_J1_gn2 <- density(par_postJ1_gn2_th[,2])
beta_J1_mn1 <- density(par_postJ1_mn1_th[,2])
beta_J1_mn2 <- density(par_postJ1_mn2_th[,2])

gamma1_J1_gn1 <- density(log(par_postJ1_gn1_th[,3]))
gamma1_J1_gn2 <- density(log(par_postJ1_gn2_th[,3]))
gamma1_J1_mn1 <- density(log(par_postJ1_mn1_th[,3]))
gamma1_J1_mn2 <- density(log(par_postJ1_mn2_th[,3]))

### - zeta

zeta0_J1_gn1 <- density(par_postJ1_gn1_th[,4])
zeta0_J1_gn2 <- density(par_postJ1_gn2_th[,4])
zeta0_J1_mn1 <- density(par_postJ1_mn1_th[,4])
zeta0_J1_mn2 <- density(par_postJ1_mn2_th[,4])

zeta1_J1_gn1 <- density(par_postJ1_gn1_th[,5])
zeta1_J1_gn2 <- density(par_postJ1_gn2_th[,5])
zeta1_J1_mn1 <- density(par_postJ1_mn1_th[,5])
zeta1_J1_mn2 <- density(par_postJ1_mn2_th[,5])

zeta2_J1_gn1 <- density(par_postJ1_gn1_th[,6])
zeta2_J1_gn2 <- density(par_postJ1_gn2_th[,6])
zeta2_J1_mn1 <- density(par_postJ1_mn1_th[,6])
zeta2_J1_mn2 <- density(par_postJ1_mn2_th[,6])

zeta3_J1_gn1 <- density(par_postJ1_gn1_th[,7])
zeta3_J1_gn2 <- density(par_postJ1_gn2_th[,7])
zeta3_J1_mn1 <- density(par_postJ1_mn1_th[,7])
zeta3_J1_mn2 <- density(par_postJ1_mn2_th[,7])

zeta4_J1_gn1 <- density(par_postJ1_gn1_th[,8])
zeta4_J1_gn2 <- density(par_postJ1_gn2_th[,8])
zeta4_J1_mn1 <- density(par_postJ1_mn1_th[,8])
zeta4_J1_mn2 <- density(par_postJ1_mn2_th[,8])

zeta5_J1_gn1 <- density(par_postJ1_gn1_th[,9])
zeta5_J1_gn2 <- density(par_postJ1_gn2_th[,9])
zeta5_J1_mn1 <- density(par_postJ1_mn1_th[,9])
zeta5_J1_mn2 <- density(par_postJ1_mn2_th[,9])

zeta6_J1_gn1 <- density(par_postJ1_gn1_th[,10])
zeta6_J1_gn2 <- density(par_postJ1_gn2_th[,10])
zeta6_J1_mn1 <- density(par_postJ1_mn1_th[,10])
zeta6_J1_mn2 <- density(par_postJ1_mn2_th[,10])

zeta7_J1_gn1 <- density(par_postJ1_gn1_th[,11])
zeta7_J1_gn2 <- density(par_postJ1_gn2_th[,11])
zeta7_J1_mn1 <- density(par_postJ1_mn1_th[,11])
zeta7_J1_mn2 <- density(par_postJ1_mn2_th[,11])

zeta8_J1_gn1 <- density(par_postJ1_gn1_th[,12])
zeta8_J1_gn2 <- density(par_postJ1_gn2_th[,12])
zeta8_J1_mn1 <- density(par_postJ1_mn1_th[,12])
zeta8_J1_mn2 <- density(par_postJ1_mn2_th[,12])

zeta9_J1_gn1 <- density(par_postJ1_gn1_th[,13])
zeta9_J1_gn2 <- density(par_postJ1_gn2_th[,13])
zeta9_J1_mn1 <- density(par_postJ1_mn1_th[,13])
zeta9_J1_mn2 <- density(par_postJ1_mn2_th[,13])

zeta10_J1_gn1 <- density(par_postJ1_gn1_th[,14])
zeta10_J1_gn2 <- density(par_postJ1_gn2_th[,14])
zeta10_J1_mn1 <- density(par_postJ1_mn1_th[,14])
zeta10_J1_mn2 <- density(par_postJ1_mn2_th[,14])

zeta11_J1_gn1 <- density(par_postJ1_gn1_th[,15])
zeta11_J1_gn2 <- density(par_postJ1_gn2_th[,15])
zeta11_J1_mn1 <- density(par_postJ1_mn1_th[,15])
zeta11_J1_mn2 <- density(par_postJ1_mn2_th[,15])

par(mfrow=c(3,3), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, log(par_postJ1_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha", ylim=c(-3.5, -2)) #ylim=c(-3.2, -1), 
lines(1:100000, log(par_postJ1_gibbs_norm2[,1]), type = "l", col="red")
lines(1:100000, log(par_postJ1_metr_norm1[,1]), type = "l", col="blue")
lines(1:100000, log(par_postJ1_metr_norm2[,1]), type = "l", col="green")

plot(c(1:100000), log(par_postJ1_gn1_erg_m[,1]), type="l", ylab = "alpha", xlab="Iteration", ylim=c(-3.5,-2)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ1_gn2_erg_m[,1]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ1_mn1_erg_m[,1]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ1_mn2_erg_m[,1]), type="l", col="green")   #, lty=4

plot(alpha_J1_gn1, xlab = "alpha", main=" ", ylim=c(0, 11)) #, main="Density of the posterior distribution of alpha"
lines(alpha_J1_gn2, col="red")
lines(alpha_J1_mn1, col="blue")
lines(alpha_J1_mn2, col="green")
abline(v=log(mean(par_postJ1_gn1_th[,1])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ1_gn2_th[,1])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ1_mn1_th[,1])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ1_mn2_th[,1])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim=c(0.08, 0.15)) #, ylim=c(0.05, 0.15)
lines(1:100000, par_postJ1_gibbs_norm2[,2], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,2], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,2], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,2], type="l", ylab = "beta", xlab="Iteration", ylim=c(0.08, 0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,2], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,2], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,2], type="l", col="green")   #, lty=4

plot(beta_J1_gn1, xlab = "beta", main=" ", ylim=c(0, 230)) #, ylim=c(0, 10.3), 
lines(beta_J1_gn2, col="red")
lines(beta_J1_mn1, col="blue")
lines(beta_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,2]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,2]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,2]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,2]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ1_gibbs_norm1[,3]), type = "l", xlab="Iteration", ylab = "psi1", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ1_gibbs_norm2[,3]), type = "l", col="red")
lines(1:100000, log(par_postJ1_metr_norm1[,3]), type = "l", col="blue")
lines(1:100000, log(par_postJ1_metr_norm2[,3]), type = "l", col="green")

plot(c(1:100000), log(par_postJ1_gn1_erg_m[,3]), type="l", ylab = "psi1", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ1_gn2_erg_m[,3]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ1_mn1_erg_m[,3]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ1_mn2_erg_m[,3]), type="l", col="green")   #, lty=4

plot(gamma1_J1_gn1, xlab = "psi1", main=" ", ylim=c(0, 8)) #, main="Density of the posterior distribution of alpha"
lines(gamma1_J1_gn2, col="red")
lines(gamma1_J1_mn1, col="blue")
lines(gamma1_J1_mn2, col="green")
abline(v=log(mean(par_postJ1_gn1_th[,3])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ1_gn2_th[,3])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ1_mn1_th[,3])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ1_mn2_th[,3])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

#zeta

### - zeta

zeta03_J1_gn1 <- density(par_postJ1_gn1_th[,4])
zeta0_J1_gn2 <- density(par_postJ1_gn2_th[,4])
zeta0_J1_mn1 <- density(par_postJ1_mn1_th[,4])
zeta0_J1_mn2 <- density(par_postJ1_mn2_th[,4])

zeta1_J1_gn1 <- density(par_postJ1_gn1_th[,5])
zeta1_J1_gn2 <- density(par_postJ1_gn2_th[,5])
zeta1_J1_mn1 <- density(par_postJ1_mn1_th[,5])
zeta1_J1_mn2 <- density(par_postJ1_mn2_th[,5])

zeta2_J1_gn1 <- density(par_postJ1_gn1_th[,6])
zeta2_J1_gn2 <- density(par_postJ1_gn2_th[,6])
zeta2_J1_mn1 <- density(par_postJ1_mn1_th[,6])
zeta2_J1_mn2 <- density(par_postJ1_mn2_th[,6])

zeta3_J1_gn1 <- density(par_postJ1_gn1_th[,7])
zeta3_J1_gn2 <- density(par_postJ1_gn2_th[,7])
zeta3_J1_mn1 <- density(par_postJ1_mn1_th[,7])
zeta3_J1_mn2 <- density(par_postJ1_mn2_th[,7])

zeta4_J1_gn1 <- density(par_postJ1_gn1_th[,8])
zeta4_J1_gn2 <- density(par_postJ1_gn2_th[,8])
zeta4_J1_mn1 <- density(par_postJ1_mn1_th[,8])
zeta4_J1_mn2 <- density(par_postJ1_mn2_th[,8])

zeta5_J1_gn1 <- density(par_postJ1_gn1_th[,9])
zeta5_J1_gn2 <- density(par_postJ1_gn2_th[,9])
zeta5_J1_mn1 <- density(par_postJ1_mn1_th[,9])
zeta5_J1_mn2 <- density(par_postJ1_mn2_th[,9])

zeta6_J1_gn1 <- density(par_postJ1_gn1_th[,10])
zeta6_J1_gn2 <- density(par_postJ1_gn2_th[,10])
zeta6_J1_mn1 <- density(par_postJ1_mn1_th[,10])
zeta6_J1_mn2 <- density(par_postJ1_mn2_th[,10])

zeta7_J1_gn1 <- density(par_postJ1_gn1_th[,11])
zeta7_J1_gn2 <- density(par_postJ1_gn2_th[,11])
zeta7_J1_mn1 <- density(par_postJ1_mn1_th[,11])
zeta7_J1_mn2 <- density(par_postJ1_mn2_th[,11])

zeta8_J1_gn1 <- density(par_postJ1_gn1_th[,12])
zeta8_J1_gn2 <- density(par_postJ1_gn2_th[,12])
zeta8_J1_mn1 <- density(par_postJ1_mn1_th[,12])
zeta8_J1_mn2 <- density(par_postJ1_mn2_th[,12])

zeta9_J1_gn1 <- density(par_postJ1_gn1_th[,13])
zeta9_J1_gn2 <- density(par_postJ1_gn2_th[,13])
zeta9_J1_mn1 <- density(par_postJ1_mn1_th[,13])
zeta9_J1_mn2 <- density(par_postJ1_mn2_th[,13])

zeta10_J1_gn1 <- density(par_postJ1_gn1_th[,14])
zeta10_J1_gn2 <- density(par_postJ1_gn2_th[,14])
zeta10_J1_mn1 <- density(par_postJ1_mn1_th[,14])
zeta10_J1_mn2 <- density(par_postJ1_mn2_th[,14])

zeta11_J1_gn1 <- density(par_postJ1_gn1_th[,15])
zeta11_J1_gn2 <- density(par_postJ1_gn2_th[,15])
zeta11_J1_mn1 <- density(par_postJ1_mn1_th[,15])
zeta11_J1_mn2 <- density(par_postJ1_mn2_th[,15])


par(mfrow=c(3,3), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, par_postJ1_gibbs_norm1[,4], type = "l", xlab="Iteration", ylab = "zeta0", ylim=c(0,0.5))
lines(1:100000, par_postJ1_gibbs_norm2[,4], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,4], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,4], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,4], type="l", ylab = "zeta0", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,4], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,4], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,4], type="l", col="green")   #, lty=4

plot(zeta0_J1_gn1, xlab = "zeta0", main=" ", ylim=c(0, 45)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta0_J1_gn2, col="red")
lines(zeta0_J1_mn1, col="blue")
lines(zeta0_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,4]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,4]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,4]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,4]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,5], type = "l", xlab="Iteration", ylab = "zeta1", ylim=c(0,0.5))
lines(1:100000, par_postJ1_gibbs_norm2[,5], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,5], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,5], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,5], type="l", ylab = "zeta1", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,5], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,5], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,5], type="l", col="green")   #, lty=4

plot(zeta1_J1_gn1, xlab = "zeta1", main=" ", ylim=c(0, 14)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1_J1_gn2, col="red")
lines(zeta1_J1_mn1, col="blue")
lines(zeta1_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,5]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,5]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,5]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,5]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,6], type = "l", xlab="Iteration", ylab = "zeta2", ylim = c(0, 0.4))
lines(1:100000, par_postJ1_gibbs_norm2[,6], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,6], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,6], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,6], type="l", ylab = "zeta2", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,6], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,6], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,6], type="l", col="green")   #, lty=4

plot(zeta2_J1_gn1, xlab = "zeta2", main=" ", ylim=c(0, 45)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta2_J1_gn2, col="red")
lines(zeta2_J1_mn1, col="blue")
lines(zeta2_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,6]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,6]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,6]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,6]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,7], type = "l", xlab="Iteration", ylab = "zeta3", ylim = c(0, 0.3))
lines(1:100000, par_postJ1_gibbs_norm2[,7], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,7], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,7], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,7], type="l", ylab = "zeta3", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,7], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,7], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,7], type="l", col="green")   #, lty=4

plot(zeta3_J1_gn1, xlab = "zeta3", main=" ", ylim=c(0, 60)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta3_J1_gn2, col="red")
lines(zeta3_J1_mn1, col="blue")
lines(zeta3_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,7]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,7]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,7]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,7]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,8], type = "l", xlab="Iteration", ylab = "zeta4", ylim = c(0, 0.5))
lines(1:100000, par_postJ1_gibbs_norm2[,8], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,8], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,8], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,8], type="l", ylab = "zeta4", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,8], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,8], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,8], type="l", col="green")   #, lty=4

plot(zeta4_J1_gn1, xlab = "zeta4", main=" ", ylim = c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta4_J1_gn2, col="red")
lines(zeta4_J1_mn1, col="blue")
lines(zeta4_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,8]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,8]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,8]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,8]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,9], type = "l", xlab="Iteration", ylab = "zeta5", ylim = c(0, 0.5))
lines(1:100000, par_postJ1_gibbs_norm2[,9], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,9], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,9], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,9], type="l", ylab = "zeta5", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,9], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,9], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,9], type="l", col="green")   #, lty=4

plot(zeta5_J1_gn1, xlab = "zeta5", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta5_J1_gn2, col="red")
lines(zeta5_J1_mn1, col="blue")
lines(zeta5_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,9]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,9]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,9]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,9]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,10], type = "l", xlab="Iteration", ylab = "zeta6", ylim = c(0, 0.1))
lines(1:100000, par_postJ1_gibbs_norm2[,10], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,10], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,10], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,10], type="l", ylab = "zeta6", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,10], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,10], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,10], type="l", col="green")   #, lty=4

plot(zeta6_J1_gn1, xlab = "zeta6", main=" ", ylim=c(0, 185)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta6_J1_gn2, col="red")
lines(zeta6_J1_mn1, col="blue")
lines(zeta6_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,10]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,10]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,10]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,10]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,11], type = "l", xlab="Iteration", ylab = "zeta7", ylim = c(0, 0.2))
lines(1:100000, par_postJ1_gibbs_norm2[,11], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,11], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,11], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,11], type="l", ylab = "zeta7", xlab="Iteration", ylim=c(0, 0.2)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,11], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,11], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,11], type="l", col="green")   #, lty=4

plot(zeta7_J1_gn1, xlab = "zeta7", main=" ", ylim=c(0, 65)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta7_J1_gn2, col="red")
lines(zeta7_J1_mn1, col="blue")
lines(zeta7_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,11]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,11]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,11]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,11]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,12], type = "l", xlab="Iteration", ylab = "zeta8", ylim = c(0, 0.2))
lines(1:100000, par_postJ1_gibbs_norm2[,12], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,12], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,12], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,12], type="l", ylab = "zeta8", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,12], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,12], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,12], type="l", col="green")   #, lty=4

plot(zeta8_J1_gn1, xlab = "zeta8", main=" ", ylim=c(0, 110)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta8_J1_gn2, col="red")
lines(zeta8_J1_mn1, col="blue")
lines(zeta8_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,12]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,12]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,12]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,12]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,13], type = "l", xlab="Iteration", ylab = "zeta9", ylim = c(0, 0.2))
lines(1:100000, par_postJ1_gibbs_norm2[,13], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,13], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,13], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,13], type="l", ylab = "zeta9", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,13], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,13], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,13], type="l", col="green")   #, lty=4

plot(zeta9_J1_gn1, xlab = "zeta9", main=" ", ylim=c(0, 90)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta9_J1_gn2, col="red")
lines(zeta9_J1_mn1, col="blue")
lines(zeta9_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,13]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,13]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,13]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,13]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,14], type = "l", xlab="Iteration", ylab = "zeta10", ylim = c(0, 0.4))
lines(1:100000, par_postJ1_gibbs_norm2[,14], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,14], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,14], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,14], type="l", ylab = "zeta10", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,14], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,14], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,14], type="l", col="green")   #, lty=4

plot(zeta10_J1_gn1, xlab = "zeta10", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta10_J1_gn2, col="red")
lines(zeta10_J1_mn1, col="blue")
lines(zeta10_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,14]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,14]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,14]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,14]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ1_gibbs_norm1[,15], type = "l", xlab="Iteration", ylab = "zeta11", ylim = c(0, 0.4))
lines(1:100000, par_postJ1_gibbs_norm2[,15], type = "l", col="red")
lines(1:100000, par_postJ1_metr_norm1[,15], type = "l", col="blue")
lines(1:100000, par_postJ1_metr_norm2[,15], type = "l", col="green")

plot(c(1:100000), par_postJ1_gn1_erg_m[,15], type="l", ylab = "zeta11", xlab="Iteration", ylim=c(0.06, 0.14)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ1_gn2_erg_m[,15], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ1_mn1_erg_m[,15], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ1_mn2_erg_m[,15], type="l", col="green")   #, lty=4

plot(zeta11_J1_gn1, xlab = "zeta11", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta11_J1_gn2, col="red")
lines(zeta11_J1_mn1, col="blue")
lines(zeta11_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,15]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,15]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,15]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,15]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

###### - zeta J1 aggregated

zeta03_J1_gn1 <- density(par_postJ1_gn1_th[,4] + par_postJ1_gn1_th[,7])
zeta03_J1_gn2 <- density(par_postJ1_gn2_th[,4] + par_postJ1_gn2_th[,7])
zeta03_J1_mn1 <- density(par_postJ1_mn1_th[,4] + par_postJ1_mn1_th[,7])
zeta03_J1_mn2 <- density(par_postJ1_mn2_th[,4] + par_postJ1_mn2_th[,7])

zeta14_J1_gn1 <- density(par_postJ1_gn1_th[,5] + par_postJ1_gn1_th[,8])
zeta14_J1_gn2 <- density(par_postJ1_gn2_th[,5] + par_postJ1_gn2_th[,8])
zeta14_J1_mn1 <- density(par_postJ1_mn1_th[,5] + par_postJ1_mn1_th[,8])
zeta14_J1_mn2 <- density(par_postJ1_mn2_th[,5] + par_postJ1_mn2_th[,8])

zeta25_J1_gn1 <- density(par_postJ1_gn1_th[,6] + par_postJ1_gn1_th[,9])
zeta25_J1_gn2 <- density(par_postJ1_gn2_th[,6] + par_postJ1_gn2_th[,9])
zeta25_J1_mn1 <- density(par_postJ1_mn1_th[,6] + par_postJ1_mn1_th[,9])
zeta25_J1_mn2 <- density(par_postJ1_mn2_th[,6] + par_postJ1_mn2_th[,9])

zeta69_J1_gn1 <- density(par_postJ1_gn1_th[,10] + par_postJ1_gn1_th[,13])
zeta69_J1_gn2 <- density(par_postJ1_gn2_th[,10] + par_postJ1_gn2_th[,13])
zeta69_J1_mn1 <- density(par_postJ1_mn1_th[,10] + par_postJ1_mn1_th[,13])
zeta69_J1_mn2 <- density(par_postJ1_mn2_th[,10] + par_postJ1_mn2_th[,13])

zeta710_J1_gn1 <- density(par_postJ1_gn1_th[,11] + par_postJ1_gn1_th[,14])
zeta710_J1_gn2 <- density(par_postJ1_gn2_th[,11] + par_postJ1_gn2_th[,14])
zeta710_J1_mn1 <- density(par_postJ1_mn1_th[,11] + par_postJ1_mn1_th[,14])
zeta710_J1_mn2 <- density(par_postJ1_mn2_th[,11] + par_postJ1_mn2_th[,14])

zeta811_J1_gn1 <- density(par_postJ1_gn1_th[,12] + par_postJ1_gn1_th[,15])
zeta811_J1_gn2 <- density(par_postJ1_gn2_th[,12] + par_postJ1_gn2_th[,15])
zeta811_J1_mn1 <- density(par_postJ1_mn1_th[,12] + par_postJ1_gn1_th[,15])
zeta811_J1_mn2 <- density(par_postJ1_mn2_th[,12] + par_postJ1_mn2_th[,15])

par(mfrow=c(6,3), mai=c(0.7,0.7, 0.2,0.2))

## - BL C0
plot(1:100000, (par_postJ1_gibbs_norm1[,4] + par_postJ1_gibbs_norm1[,7]), type = "l", xlab="Iteration", ylab = "z0+z3", ylim=c(0,0.4))
lines(1:100000,(par_postJ1_gibbs_norm2[,4] + par_postJ1_gibbs_norm2[,7]), type = "l", col="red")
lines(1:100000, (par_postJ1_metr_norm1[,4] + par_postJ1_metr_norm1[,7]), type = "l", col="blue")
lines(1:100000, (par_postJ1_metr_norm2[,4] + par_postJ1_metr_norm2[,7]), type = "l", col="green")

plot(c(1:100000), (par_postJ1_gn1_erg_m[,4] + par_postJ1_gn1_erg_m[,7]), type="l", ylab = "z0+z3", xlab="Iteration", ylim=c(0,0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ1_gn2_erg_m[,4] + par_postJ1_gn2_erg_m[,7]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ1_mn1_erg_m[,4] + par_postJ1_mn1_erg_m[,7]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ1_mn2_erg_m[,4] + par_postJ1_mn2_erg_m[,7]), type="l", col="green")   #, lty=4

plot(zeta03_J1_gn1, xlab = "z0+z3", main=" ", ylim=c(0, 45)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta03_J1_mn2, col="red")
lines(zeta03_J1_gn1, col="blue")
lines(zeta03_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,4] + par_postJ1_gn1_th[,7]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,4] + par_postJ1_gn2_th[,7]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,4] + par_postJ1_mn1_th[,7]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,4] + par_postJ1_mn2_th[,7]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C0
plot(1:100000, (par_postJ1_gibbs_norm1[,10] + par_postJ1_gibbs_norm1[,13]), type = "l", xlab="Iteration", ylab = "z6+z9", ylim=c(0,0.2))
lines(1:100000,(par_postJ1_gibbs_norm2[,10] + par_postJ1_gibbs_norm2[,13]), type = "l", col="red")
lines(1:100000, (par_postJ1_metr_norm1[,10] + par_postJ1_metr_norm1[,13]), type = "l", col="blue")
lines(1:100000, (par_postJ1_metr_norm2[,10] + par_postJ1_metr_norm2[,13]), type = "l", col="green")

plot(c(1:100000), (par_postJ1_gn1_erg_m[,10] + par_postJ1_gn1_erg_m[,13]), type="l", ylab = "z6+z9", xlab="Iteration", ylim=c(0, 0.2)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ1_gn2_erg_m[,10] + par_postJ1_gn2_erg_m[,13]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ1_mn1_erg_m[,10] + par_postJ1_mn1_erg_m[,13]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ1_mn2_erg_m[,10] + par_postJ1_mn2_erg_m[,13]), type="l", col="green")   #, lty=4

plot(zeta69_J1_gn1, xlab = "z6+z9", main=" ", ylim=c(0, 50)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta69_J1_gn2, col="red")
lines(zeta69_J1_mn1, col="blue")
lines(zeta69_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,10] + par_postJ1_gn1_th[,13]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,10] + par_postJ1_gn2_th[,13]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,10] + par_postJ1_mn1_th[,13]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,10] + par_postJ1_mn2_th[,13]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BL C1
plot(1:100000, (par_postJ1_gibbs_norm1[,5] + par_postJ1_gibbs_norm1[,8]), type = "l", xlab="Iteration", ylab = "z1+z4", ylim=c(0.25,0.75))
lines(1:100000,(par_postJ1_gibbs_norm2[,5] + par_postJ1_gibbs_norm2[,8]), type = "l", col="red")
lines(1:100000, (par_postJ1_metr_norm1[,5] + par_postJ1_metr_norm1[,8]), type = "l", col="blue")
lines(1:100000, (par_postJ1_metr_norm2[,5] + par_postJ1_metr_norm2[,8]), type = "l", col="green")

plot(c(1:100000), (par_postJ1_gn1_erg_m[,5] + par_postJ1_gn1_erg_m[,8]), type="l", ylab = "z1+z4", xlab="Iteration", ylim=c(0.25,0.75)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ1_gn2_erg_m[,5] + par_postJ1_gn2_erg_m[,8]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ1_mn1_erg_m[,5] + par_postJ1_mn1_erg_m[,8]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ1_mn2_erg_m[,5] + par_postJ1_mn2_erg_m[,8]), type="l", col="green")   #, lty=4

plot(zeta14_J1_gn1, xlab = "z1+z4", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta14_J1_gn2, col="red")
lines(zeta14_J1_mn1, col="blue")
lines(zeta14_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,5] + par_postJ1_gn1_th[,8]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,5] + par_postJ1_gn2_th[,8]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,5] + par_postJ1_mn1_th[,8]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,5] + par_postJ1_mn2_th[,8]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C1
plot(1:100000, (par_postJ1_gibbs_norm1[,11] + par_postJ1_gibbs_norm1[,14]), type = "l", xlab="Iteration", ylab = "z7+z10", ylim=c(0,0.4))
lines(1:100000,(par_postJ1_gibbs_norm2[,11] + par_postJ1_gibbs_norm2[,14]), type = "l", col="red")
lines(1:100000, (par_postJ1_metr_norm1[,11] + par_postJ1_metr_norm1[,14]), type = "l", col="blue")
lines(1:100000, (par_postJ1_metr_norm2[,11] + par_postJ1_metr_norm2[,14]), type = "l", col="green")

plot(c(1:100000), (par_postJ1_gn1_erg_m[,11] + par_postJ1_gn1_erg_m[,14]), type="l", ylab = "z7+z10", xlab="Iteration", ylim=c(0, 0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ1_gn2_erg_m[,11] + par_postJ1_gn2_erg_m[,14]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ1_mn1_erg_m[,11] + par_postJ1_mn1_erg_m[,14]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ1_mn2_erg_m[,11] + par_postJ1_mn2_erg_m[,14]), type="l", col="green")   #, lty=4

plot(zeta710_J1_gn1, xlab = "z7+z10", main=" ", ylim=c(0, 22)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta710_J1_gn2, col="red")
lines(zeta710_J1_mn1, col="blue")
lines(zeta710_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,11] + par_postJ1_gn1_th[,14]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,11] + par_postJ1_gn2_th[,14]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,11] + par_postJ1_mn1_th[,14]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,11] + par_postJ1_mn2_th[,14]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BL C2
plot(1:100000, (par_postJ1_gibbs_norm1[,6] + par_postJ1_gibbs_norm1[,9]), type = "l", xlab="Iteration", ylab = "z2+z5", ylim=c(0,0.4))
lines(1:100000,(par_postJ1_gibbs_norm2[,6] + par_postJ1_gibbs_norm2[,9]), type = "l", col="red")
lines(1:100000, (par_postJ1_metr_norm1[,6] + par_postJ1_metr_norm1[,9]), type = "l", col="blue")
lines(1:100000, (par_postJ1_metr_norm2[,6] + par_postJ1_metr_norm2[,9]), type = "l", col="green")

plot(c(1:100000), (par_postJ1_gn1_erg_m[,6] + par_postJ1_gn1_erg_m[,9]), type="l", ylab = "z2+z5", xlab="Iteration", ylim=c(0,0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ1_gn2_erg_m[,6] + par_postJ1_gn2_erg_m[,9]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ1_mn1_erg_m[,6] + par_postJ1_mn1_erg_m[,9]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ1_mn2_erg_m[,6] + par_postJ1_mn2_erg_m[,9]), type="l", col="green")   #, lty=4

plot(zeta25_J1_gn1, xlab = "z2+z5", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta25_J1_gn2, col="red")
lines(zeta25_J1_mn1, col="blue")
lines(zeta25_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,6] + par_postJ1_gn1_th[,9]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,6] + par_postJ1_gn2_th[,9]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,6] + par_postJ1_mn1_th[,9]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,6] + par_postJ1_mn2_th[,9]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C2
plot(1:100000, (par_postJ1_gibbs_norm1[,12] + par_postJ1_gibbs_norm1[,15]), type = "l", xlab="Iteration", ylab = "z8+z11", ylim=c(0,0.4))
lines(1:100000,(par_postJ1_gibbs_norm2[,12] + par_postJ1_gibbs_norm2[,15]), type = "l", col="red")
lines(1:100000, (par_postJ1_metr_norm1[,12] + par_postJ1_metr_norm1[,15]), type = "l", col="blue")
lines(1:100000, (par_postJ1_metr_norm2[,12] + par_postJ1_metr_norm2[,15]), type = "l", col="green")

plot(c(1:100000), (par_postJ1_gn1_erg_m[,12] + par_postJ1_gn1_erg_m[,15]), type="l", ylab = "z8+z11", xlab="Iteration", ylim=c(0, 0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ1_gn2_erg_m[,12] + par_postJ1_gn2_erg_m[,15]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ1_mn1_erg_m[,12] + par_postJ1_mn1_erg_m[,15]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ1_mn2_erg_m[,12] + par_postJ1_mn2_erg_m[,15]), type="l", col="green")   #, lty=4

plot(zeta811_J1_gn1, xlab = "z8+z11", main=" ", ylim=c(0, 22)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta811_J1_gn2, col="red")
lines(zeta811_J1_mn1, col="blue")
lines(zeta811_J1_mn2, col="green")
abline(v=mean(par_postJ1_gn1_th[,12] + par_postJ1_gn1_th[,15]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_gn2_th[,12] + par_postJ1_gn2_th[,15]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn1_th[,12] + par_postJ1_mn1_th[,15]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ1_mn2_th[,12] + par_postJ1_mn2_th[,15]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)


# - J = 2
par(mfrow=c(4,3), mai=c(0.7,0.7, 0.2,0.2))

## - Thinned sequence

par_postJ2_gn1_th <- par_postJ2_gibbs_norm1[seq(10010, 100000, 10),]
par_postJ2_gn2_th <- par_postJ2_gibbs_norm2[seq(10010, 100000, 10),]
par_postJ2_mn1_th <- par_postJ2_metr_norm1[seq(10010, 100000, 10),]
par_postJ2_mn2_th <- par_postJ2_metr_norm2[seq(10010, 100000, 10),]

## - Density

### - tau
alpha_J2_gn1 <- density(log(par_postJ2_gn1_th[,1]))
alpha_J2_gn2 <- density(log(par_postJ2_gn2_th[,1]))
alpha_J2_mn1 <- density(log(par_postJ2_mn1_th[,1]))
alpha_J2_mn2 <- density(log(par_postJ2_mn2_th[,1]))

beta_J2_gn1 <- density(par_postJ2_gn1_th[,2])
beta_J2_gn2 <- density(par_postJ2_gn2_th[,2])
beta_J2_mn1 <- density(par_postJ2_mn1_th[,2])
beta_J2_mn2 <- density(par_postJ2_mn2_th[,2])

gamma1_J2_gn1 <- density(log(par_postJ2_gn1_th[,3]))
gamma1_J2_gn2 <- density(log(par_postJ2_gn2_th[,3]))
gamma1_J2_mn1 <- density(log(par_postJ2_mn1_th[,3]))
gamma1_J2_mn2 <- density(log(par_postJ2_mn2_th[,3]))

gamma2_J2_gn1 <- density(log(par_postJ2_gn1_th[,4]))
gamma2_J2_gn2 <- density(log(par_postJ2_gn2_th[,4]))
gamma2_J2_mn1 <- density(log(par_postJ2_mn1_th[,4]))
gamma2_J2_mn2 <- density(log(par_postJ2_mn2_th[,4]))

### - zeta

zeta0_J2_gn1 <- density(par_postJ2_gn1_th[,5])
zeta0_J2_gn2 <- density(par_postJ2_gn2_th[,5])
zeta0_J2_mn1 <- density(par_postJ2_mn1_th[,5])
zeta0_J2_mn2 <- density(par_postJ2_mn2_th[,5])

zeta1_J2_gn1 <- density(par_postJ2_gn1_th[,6])
zeta1_J2_gn2 <- density(par_postJ2_gn2_th[,6])
zeta1_J2_mn1 <- density(par_postJ2_mn1_th[,6])
zeta1_J2_mn2 <- density(par_postJ2_mn2_th[,6])

zeta2_J2_gn1 <- density(par_postJ2_gn1_th[,7])
zeta2_J2_gn2 <- density(par_postJ2_gn2_th[,7])
zeta2_J2_mn1 <- density(par_postJ2_mn1_th[,7])
zeta2_J2_mn2 <- density(par_postJ2_mn2_th[,7])

zeta3_J2_gn1 <- density(par_postJ2_gn1_th[,8])
zeta3_J2_gn2 <- density(par_postJ2_gn2_th[,8])
zeta3_J2_mn1 <- density(par_postJ2_mn1_th[,8])
zeta3_J2_mn2 <- density(par_postJ2_mn2_th[,8])

zeta4_J2_gn1 <- density(par_postJ2_gn1_th[,9])
zeta4_J2_gn2 <- density(par_postJ2_gn2_th[,9])
zeta4_J2_mn1 <- density(par_postJ2_mn1_th[,9])
zeta4_J2_mn2 <- density(par_postJ2_mn2_th[,9])

zeta5_J2_gn1 <- density(par_postJ2_gn1_th[,10])
zeta5_J2_gn2 <- density(par_postJ2_gn2_th[,10])
zeta5_J2_mn1 <- density(par_postJ2_mn1_th[,10])
zeta5_J2_mn2 <- density(par_postJ2_mn2_th[,10])

zeta6_J2_gn1 <- density(par_postJ2_gn1_th[,11])
zeta6_J2_gn2 <- density(par_postJ2_gn2_th[,11])
zeta6_J2_mn1 <- density(par_postJ2_mn1_th[,11])
zeta6_J2_mn2 <- density(par_postJ2_mn2_th[,11])

zeta7_J2_gn1 <- density(par_postJ2_gn1_th[,12])
zeta7_J2_gn2 <- density(par_postJ2_gn2_th[,12])
zeta7_J2_mn1 <- density(par_postJ2_mn1_th[,12])
zeta7_J2_mn2 <- density(par_postJ2_mn2_th[,12])

zeta8_J2_gn1 <- density(par_postJ2_gn1_th[,13])
zeta8_J2_gn2 <- density(par_postJ2_gn2_th[,13])
zeta8_J2_mn1 <- density(par_postJ2_mn1_th[,13])
zeta8_J2_mn2 <- density(par_postJ2_mn2_th[,13])

zeta9_J2_gn1 <- density(par_postJ2_gn1_th[,14])
zeta9_J2_gn2 <- density(par_postJ2_gn2_th[,14])
zeta9_J2_mn1 <- density(par_postJ2_mn1_th[,14])
zeta9_J2_mn2 <- density(par_postJ2_mn2_th[,14])

zeta10_J2_gn1 <- density(par_postJ2_gn1_th[,15])
zeta10_J2_gn2 <- density(par_postJ2_gn2_th[,15])
zeta10_J2_mn1 <- density(par_postJ2_mn1_th[,15])
zeta10_J2_mn2 <- density(par_postJ2_mn2_th[,15])

zeta11_J2_gn1 <- density(par_postJ2_gn1_th[,16])
zeta11_J2_gn2 <- density(par_postJ2_gn2_th[,16])
zeta11_J2_mn1 <- density(par_postJ2_mn1_th[,16])
zeta11_J2_mn2 <- density(par_postJ2_mn2_th[,16])

zeta12_J2_gn1 <- density(par_postJ2_gn1_th[,17])
zeta12_J2_gn2 <- density(par_postJ2_gn2_th[,17])
zeta12_J2_mn1 <- density(par_postJ2_mn1_th[,17])
zeta12_J2_mn2 <- density(par_postJ2_mn2_th[,17])

zeta13_J2_gn1 <- density(par_postJ2_gn1_th[,18])
zeta13_J2_gn2 <- density(par_postJ2_gn2_th[,18])
zeta13_J2_mn1 <- density(par_postJ2_mn1_th[,18])
zeta13_J2_mn2 <- density(par_postJ2_mn2_th[,18])

zeta14_J2_gn1 <- density(par_postJ2_gn1_th[,19])
zeta14_J2_gn2 <- density(par_postJ2_gn2_th[,19])
zeta14_J2_mn1 <- density(par_postJ2_mn1_th[,19])
zeta14_J2_mn2 <- density(par_postJ2_mn2_th[,19])

zeta15_J2_gn1 <- density(par_postJ2_gn1_th[,20])
zeta15_J2_gn2 <- density(par_postJ2_gn2_th[,20])
zeta15_J2_mn1 <- density(par_postJ2_mn1_th[,20])
zeta15_J2_mn2 <- density(par_postJ2_mn2_th[,20])

zeta16_J2_gn1 <- density(par_postJ2_gn1_th[,21])
zeta16_J2_gn2 <- density(par_postJ2_gn2_th[,21])
zeta16_J2_mn1 <- density(par_postJ2_mn1_th[,21])
zeta16_J2_mn2 <- density(par_postJ2_mn2_th[,21])

zeta17_J2_gn1 <- density(par_postJ2_gn1_th[,22])
zeta17_J2_gn2 <- density(par_postJ2_gn2_th[,22])
zeta17_J2_mn1 <- density(par_postJ2_mn1_th[,22])
zeta17_J2_mn2 <- density(par_postJ2_mn2_th[,22])

plot(1:100000, log(par_postJ2_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha", ylim=c(-3.5, -2.5)) #ylim=c(-3.2, -1), 
lines(1:100000, log(par_postJ2_gibbs_norm2[,1]), type = "l", col="red")
lines(1:100000, log(par_postJ2_metr_norm1[,1]), type = "l", col="blue")
lines(1:100000, log(par_postJ2_metr_norm2[,1]), type = "l", col="green")
#abline(h=log(phi0_p), lty=2)

plot(c(1:100000), log(par_postJ2_gn1_erg_m[,1]), type="l", ylab = "alpha", xlab="Iteration", ylim=c(-3.5,-2.5)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ2_gn2_erg_m[,1]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ2_mn1_erg_m[,1]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ2_mn2_erg_m[,1]), type="l", col="green")   #, lty=4

plot(alpha_J2_gn1, xlab = "alpha", main=" ", ylim=c(0, 8.5)) #, main="Density of the posterior distribution of alpha"
lines(alpha_J2_gn2, col="red")
lines(alpha_J2_mn1, col="blue")
lines(alpha_J2_mn2, col="green")
abline(v=log(mean(par_postJ2_gn1_th[,1])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_gn2_th[,1])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_mn1_th[,1])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_mn2_th[,1])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim = c(0.06, 0.15)) #, ylim=c(0.05, 0.15)
lines(1:100000, par_postJ2_gibbs_norm2[,2], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,2], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,2], type = "l", col="green")
#abline(h=beta_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,2], type="l", ylab = "beta", xlab="Iteration", ylim=c(0.06,0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,2], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,2], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,2], type="l", col="green")   #, lty=4

plot(beta_J2_gn1, xlab = "beta", main=" ", ylim=c(0, 225)) #, ylim=c(0, 10.3), 
lines(beta_J2_gn2, col="red")
lines(beta_J2_mn1, col="blue")
lines(beta_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,2]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,2]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,2]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,2]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ2_gibbs_norm1[,3]), type = "l", xlab="Iteration", ylab = "psi1", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ2_gibbs_norm2[,3]), type = "l", col="red")
lines(1:100000, log(par_postJ2_metr_norm1[,3]), type = "l", col="blue")
lines(1:100000, log(par_postJ2_metr_norm2[,3]), type = "l", col="green")
#abline(h=log(phi1_p), lty=2)

plot(c(1:100000), log(par_postJ2_gn1_erg_m[,3]), type="l", ylab = "psi1", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ2_gn2_erg_m[,3]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ2_mn1_erg_m[,3]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ2_mn2_erg_m[,3]), type="l", col="green")   #, lty=4

plot(gamma1_J2_gn1, xlab = "psi1", main=" ", ylim=c(0, 5)) #, main="Density of the posterior distribution of alpha"
lines(gamma1_J2_gn2, col="red")
lines(gamma1_J2_mn1, col="blue")
lines(gamma1_J2_mn2, col="green")
abline(v=log(mean(par_postJ2_gn1_th[,3])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_gn2_th[,3])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_mn1_th[,3])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_mn2_th[,3])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ2_gibbs_norm1[,4]), type = "l", xlab="Iteration", ylab = "psi2", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ2_gibbs_norm2[,4]), type = "l", col="red")
lines(1:100000, log(par_postJ2_metr_norm1[,4]), type = "l", col="blue")
lines(1:100000, log(par_postJ2_metr_norm2[,4]), type = "l", col="green")
#abline(h=log(phi2_p), lty=2)

plot(c(1:100000), log(par_postJ2_gn1_erg_m[,4]), type="l", ylab = "psi2", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ2_gn2_erg_m[,4]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ2_mn1_erg_m[,4]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ2_mn2_erg_m[,4]), type="l", col="green")   #, lty=4

plot(gamma2_J2_gn1, xlab = "psi2", main=" ", ylim=c(0, 4)) #, main="Density of the posterior distribution of alpha"
lines(gamma2_J2_gn2, col="red")
lines(gamma2_J2_mn1, col="blue")
lines(gamma2_J2_mn2, col="green")
abline(v=log(mean(par_postJ2_gn1_th[,4])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_gn2_th[,4])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_mn1_th[,4])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ2_mn2_th[,4])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)


#zeta - save plots 9/8 - same for tau
par(mfrow=c(6,3), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, par_postJ2_gibbs_norm1[,5], type = "l", xlab="Iteration", ylab = "zeta0", ylim = c(0, 0.5))
lines(1:100000, par_postJ2_gibbs_norm2[,5], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,5], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,5], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,5], type="l", ylab = "zeta0", xlab="Iteration", ylim=c(0,0.5)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,5], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,5], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,5], type="l", col="green")   #, lty=4

plot(zeta0_J2_gn1, xlab = "zeta0", main=" ", ylim=c(0, 27)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta0_J2_gn2, col="red")
lines(zeta0_J2_mn1, col="blue")
lines(zeta0_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,5]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,5]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,5]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,5]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,6], type = "l", xlab="Iteration", ylab = "zeta1", ylim = c(0, 0.5))
lines(1:100000, par_postJ2_gibbs_norm2[,6], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,6], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,6], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,6], type="l", ylab = "zeta1", xlab="Iteration", ylim=c(0,0.5)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,6], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,6], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,6], type="l", col="green")   #, lty=4

plot(zeta1_J2_gn1, xlab = "zeta1", main=" ", ylim=c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1_J2_gn2, col="red")
lines(zeta1_J2_mn1, col="blue")
lines(zeta1_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,6]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,6]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,6]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,6]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,7], type = "l", xlab="Iteration", ylab = "zeta2", ylim = c(0, 0.15))
lines(1:100000, par_postJ2_gibbs_norm2[,7], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,7], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,7], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,7], type="l", ylab = "zeta2", xlab="Iteration", ylim=c(0,0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,7], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,7], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,7], type="l", col="green")   #, lty=4

plot(zeta2_J2_gn1, xlab = "zeta2", main=" ", ylim=c(0, 80)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta2_J2_gn2, col="red")
lines(zeta2_J2_mn1, col="blue")
lines(zeta2_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,7]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,7]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,7]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,7]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,8], type = "l", xlab="Iteration", ylab = "zeta3", ylim = c(0, 0.15))
lines(1:100000, par_postJ2_gibbs_norm2[,8], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,8], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,8], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,8], type="l", ylab = "zeta3", xlab="Iteration", ylim=c(0,0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,8], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,8], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,8], type="l", col="green")   #, lty=4

plot(zeta3_J2_gn1, xlab = "zeta3", main=" ", ylim=c(0, 35)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta3_J2_gn2, col="red")
lines(zeta3_J2_mn1, col="blue")
lines(zeta3_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,8]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,8]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,8]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,8]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,9], type = "l", xlab="Iteration", ylab = "zeta4", ylim = c(0, 0.5))
lines(1:100000, par_postJ2_gibbs_norm2[,9], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,9], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,9], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,9], type="l", ylab = "zeta4", xlab="Iteration", ylim=c(0,0.5)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,9], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,9], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,9], type="l", col="green")   #, lty=4

plot(zeta4_J2_gn1, xlab = "zeta4", main=" ", ylim=c(0, 12)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta4_J2_gn2, col="red")
lines(zeta4_J2_mn1, col="blue")
lines(zeta4_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,9]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,9]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,9]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,9]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,10], type = "l", xlab="Iteration", ylab = "zeta5", ylim=c(0,0.2))
lines(1:100000, par_postJ2_gibbs_norm2[,10], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,10], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,10], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,10], type="l", ylab = "zeta5", xlab="Iteration", ylim=c(0,0.2)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,10], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,10], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,10], type="l", col="green")   #, lty=4

plot(zeta5_J2_gn1, xlab = "zeta5", main=" ", ylim=c(0, 23)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta5_J2_gn2, col="red")
lines(zeta5_J2_mn1, col="blue")
lines(zeta5_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,10]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,10]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,10]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,10]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,11], type = "l", xlab="Iteration", ylab = "zeta6", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,11], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,11], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,11], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,11], type="l", ylab = "zeta6", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,11], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,11], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,11], type="l", col="green")   #, lty=4

plot(zeta6_J2_gn1, xlab = "zeta6", main=" ", ylim=c(0, 75)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta6_J2_gn2, col="red")
lines(zeta6_J2_mn1, col="blue")
lines(zeta6_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,11]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,11]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,11]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,11]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,12], type = "l", xlab="Iteration", ylab = "zeta7", ylim=c(0,0.3))
lines(1:100000, par_postJ2_gibbs_norm2[,12], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,12], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,12], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,12], type="l", ylab = "zeta7", xlab="Iteration", ylim=c(0,0.3)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,12], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,12], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,12], type="l", col="green")   #, lty=4

plot(zeta7_J2_gn1, xlab = "zeta7", main=" ", ylim=c(0, 18)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta7_J2_gn2, col="red")
lines(zeta7_J2_mn1, col="blue")
lines(zeta7_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,12]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,12]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,12]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,12]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,13], type = "l", xlab="Iteration", ylab = "zeta8", ylim=c(0,0.15))
lines(1:100000, par_postJ2_gibbs_norm2[,13], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,13], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,13], type = "l", col="green")

plot(c(1:100000), par_postJ2_gn1_erg_m[,13], type="l", ylab = "zeta8", xlab="Iteration", ylim=c(0,0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,13], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,13], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,13], type="l", col="green")   #, lty=4

plot(zeta8_J2_gn1, xlab = "zeta8", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta8_J2_gn2, col="red")
lines(zeta8_J2_mn1, col="blue")
lines(zeta8_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,13]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,13]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,13]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,13]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,14], type = "l", xlab="Iteration", ylab = "zeta9", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,14], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,14], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,14], type = "l", col="green")
#abline(h=z9_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,14], type="l", ylab = "zeta9", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,14], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,14], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,14], type="l", col="green")   #, lty=4

plot(zeta9_J2_gn1, xlab = "zeta9", main=" ", ylim=c(0, 155)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta9_J2_gn2, col="red")
lines(zeta9_J2_mn1, col="blue")
lines(zeta9_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,14]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,14]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,14]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,14]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,15], type = "l", xlab="Iteration", ylab = "zeta10", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,15], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,15], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,15], type = "l", col="green")
#abline(h=z10_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,15], type="l", ylab = "zeta10", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,15], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,15], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,15], type="l", col="green")   #, lty=4

plot(zeta10_J2_gn1, xlab = "zeta10", main=" ", ylim=c(0, 80)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta10_J2_gn2, col="red")
lines(zeta10_J2_mn1, col="blue")
lines(zeta10_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,15]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,15]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,15]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,15]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,16], type = "l", xlab="Iteration", ylab = "zeta11", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,16], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,16], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,16], type = "l", col="green")
#abline(h=z11_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,16], type="l", ylab = "zeta11", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,16], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,16], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,16], type="l", col="green")   #, lty=4

plot(zeta11_J2_gn1, xlab = "zeta11", main=" ", ylim=c(0, 250)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta11_J2_gn2, col="red")
lines(zeta11_J2_mn1, col="blue")
lines(zeta11_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,16]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,16]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,16]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,16]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,17], type = "l", xlab="Iteration", ylab = "zeta12", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,17], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,17], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,17], type = "l", col="green")
#abline(h=z12_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,17], type="l", ylab = "zeta12", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,17], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,17], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,17], type="l", col="green")   #, lty=4

plot(zeta12_J2_gn1, xlab = "zeta12", main=" ", ylim=c(0, 140)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta12_J2_gn2, col="red")
lines(zeta12_J2_mn1, col="blue")
lines(zeta12_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,17]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,17]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,17]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,17]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,18], type = "l", xlab="Iteration", ylab = "zeta13", ylim=c(0,0.2))
lines(1:100000, par_postJ2_gibbs_norm2[,18], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,18], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,18], type = "l", col="green")
#abline(h=z13_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,18], type="l", ylab = "zeta13", xlab="Iteration", ylim=c(0,0.2)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,18], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,18], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,18], type="l", col="green")   #, lty=4

plot(zeta13_J2_gn1, xlab = "zeta13", main=" ", ylim=c(0, 22)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta13_J2_gn2, col="red")
lines(zeta13_J2_mn1, col="blue")
lines(zeta13_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,18]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,18]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,18]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,18]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,19], type = "l", xlab="Iteration", ylab = "zeta14", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,19], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,19], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,19], type = "l", col="green")
#abline(h=z14_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,19], type="l", ylab = "zeta14", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,19], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,19], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,19], type="l", col="green")   #, lty=4

plot(zeta14_J2_gn1, xlab = "zeta14", main=" ", ylim=c(0, 120)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta14_J2_gn2, col="red")
lines(zeta14_J2_mn1, col="blue")
lines(zeta14_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,19]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,19]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,19]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,19]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,20], type = "l", xlab="Iteration", ylab = "zeta15", ylim=c(0,0.1))
lines(1:100000, par_postJ2_gibbs_norm2[,20], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,20], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,20], type = "l", col="green")
#abline(h=z15_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,20], type="l", ylab = "zeta15", xlab="Iteration", ylim=c(0,0.1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,20], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,20], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,20], type="l", col="green")   #, lty=4

plot(zeta15_J2_gn1, xlab = "zeta15", main=" ", ylim=c(0, 75)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta15_J2_gn2, col="red")
lines(zeta15_J2_mn1, col="blue")
lines(zeta15_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,20]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,20]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,20]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,20]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,21], type = "l", xlab="Iteration", ylab = "zeta16", ylim=c(0,0.2))
lines(1:100000, par_postJ2_gibbs_norm2[,21], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,21], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,21], type = "l", col="green")
#abline(h=z16_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,21], type="l", ylab = "zeta16", xlab="Iteration", ylim=c(0,0.2)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,21], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,21], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,21], type="l", col="green")   #, lty=4

plot(zeta16_J2_gn1, xlab = "zeta16", main=" ", ylim=c(0, 30)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta16_J2_gn2, col="red")
lines(zeta16_J2_mn1, col="blue")
lines(zeta16_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,21]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,21]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,21]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,21]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,22], type = "l", xlab="Iteration", ylab = "zeta17", ylim=c(0,0.25))
lines(1:100000, par_postJ2_gibbs_norm2[,22], type = "l", col="red")
lines(1:100000, par_postJ2_metr_norm1[,22], type = "l", col="blue")
lines(1:100000, par_postJ2_metr_norm2[,22], type = "l", col="green")
#abline(h=z17_p, lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,22], type="l", ylab = "zeta17", xlab="Iteration", ylim=c(0,0.25)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn2_erg_m[,22], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ2_mn1_erg_m[,22], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ2_mn2_erg_m[,22], type="l", col="green")   #, lty=4

plot(zeta17_J2_gn1, xlab = "zeta17", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta17_J2_gn2, col="red")
lines(zeta17_J2_mn1, col="blue")
lines(zeta17_J2_mn2, col="green")
abline(v=mean(par_postJ2_gn1_th[,22]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_gn2_th[,22]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn1_th[,22]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ2_mn2_th[,22]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

# - J = 3

## - Thinned sequence

par_postJ3_gn1_th <- par_postJ3_gibbs_norm1[seq(10010, 100000, 10),]
par_postJ3_gn2_th <- par_postJ3_gibbs_norm2[seq(10010, 100000, 10),]
par_postJ3_mn1_th <- par_postJ3_metr_norm1[seq(10010, 100000, 10),]
par_postJ3_mn2_th <- par_postJ3_metr_norm2[seq(10010, 100000, 10),]

## - Density

### - tau
alpha_J3_gn1 <- density(log(par_postJ3_gn1_th[,1]))
alpha_J3_gn2 <- density(log(par_postJ3_gn2_th[,1]))
alpha_J3_mn1 <- density(log(par_postJ3_mn1_th[,1]))
alpha_J3_mn2 <- density(log(par_postJ3_mn2_th[,1]))

beta_J3_gn1 <- density(par_postJ3_gn1_th[,2])
beta_J3_gn2 <- density(par_postJ3_gn2_th[,2])
beta_J3_mn1 <- density(par_postJ3_mn1_th[,2])
beta_J3_mn2 <- density(par_postJ3_mn2_th[,2])

gamma1_J3_gn1 <- density(log(par_postJ3_gn1_th[,3]))
gamma1_J3_gn2 <- density(log(par_postJ3_gn2_th[,3]))
gamma1_J3_mn1 <- density(log(par_postJ3_mn1_th[,3]))
gamma1_J3_mn2 <- density(log(par_postJ3_mn2_th[,3]))

gamma2_J3_gn1 <- density(log(par_postJ3_gn1_th[,4]))
gamma2_J3_gn2 <- density(log(par_postJ3_gn2_th[,4]))
gamma2_J3_mn1 <- density(log(par_postJ3_mn1_th[,4]))
gamma2_J3_mn2 <- density(log(par_postJ3_mn2_th[,4]))

gamma3_J3_gn1 <- density(log(par_postJ3_gn1_th[,5]))
gamma3_J3_gn2 <- density(log(par_postJ3_gn2_th[,5]))
gamma3_J3_mn1 <- density(log(par_postJ3_mn1_th[,5]))
gamma3_J3_mn2 <- density(log(par_postJ3_mn2_th[,5]))

### - zeta

zeta0_J3_gn1 <- density(par_postJ3_gn1_th[,6])
zeta0_J3_gn2 <- density(par_postJ3_gn2_th[,6])
zeta0_J3_mn1 <- density(par_postJ3_mn1_th[,6])
zeta0_J3_mn2 <- density(par_postJ3_mn2_th[,6])

zeta1_J3_gn1 <- density(par_postJ3_gn1_th[,7])
zeta1_J3_gn2 <- density(par_postJ3_gn2_th[,7])
zeta1_J3_mn1 <- density(par_postJ3_mn1_th[,7])
zeta1_J3_mn2 <- density(par_postJ3_mn2_th[,7])

zeta2_J3_gn1 <- density(par_postJ3_gn1_th[,8])
zeta2_J3_gn2 <- density(par_postJ3_gn2_th[,8])
zeta2_J3_mn1 <- density(par_postJ3_mn1_th[,8])
zeta2_J3_mn2 <- density(par_postJ3_mn2_th[,8])

zeta3_J3_gn1 <- density(par_postJ3_gn1_th[,9])
zeta3_J3_gn2 <- density(par_postJ3_gn2_th[,9])
zeta3_J3_mn1 <- density(par_postJ3_mn1_th[,9])
zeta3_J3_mn2 <- density(par_postJ3_mn2_th[,9])

zeta4_J3_gn1 <- density(par_postJ3_gn1_th[,10])
zeta4_J3_gn2 <- density(par_postJ3_gn2_th[,10])
zeta4_J3_mn1 <- density(par_postJ3_mn1_th[,10])
zeta4_J3_mn2 <- density(par_postJ3_mn2_th[,10])

zeta5_J3_gn1 <- density(par_postJ3_gn1_th[,11])
zeta5_J3_gn2 <- density(par_postJ3_gn2_th[,11])
zeta5_J3_mn1 <- density(par_postJ3_mn1_th[,11])
zeta5_J3_mn2 <- density(par_postJ3_mn2_th[,11])

zeta6_J3_gn1 <- density(par_postJ3_gn1_th[,12])
zeta6_J3_gn2 <- density(par_postJ3_gn2_th[,12])
zeta6_J3_mn1 <- density(par_postJ3_mn1_th[,12])
zeta6_J3_mn2 <- density(par_postJ3_mn2_th[,12])

zeta7_J3_gn1 <- density(par_postJ3_gn1_th[,13])
zeta7_J3_gn2 <- density(par_postJ3_gn2_th[,13])
zeta7_J3_mn1 <- density(par_postJ3_mn1_th[,13])
zeta7_J3_mn2 <- density(par_postJ3_mn2_th[,13])

zeta8_J3_gn1 <- density(par_postJ3_gn1_th[,14])
zeta8_J3_gn2 <- density(par_postJ3_gn2_th[,14])
zeta8_J3_mn1 <- density(par_postJ3_mn1_th[,14])
zeta8_J3_mn2 <- density(par_postJ3_mn2_th[,14])

zeta9_J3_gn1 <- density(par_postJ3_gn1_th[,15])
zeta9_J3_gn2 <- density(par_postJ3_gn2_th[,15])
zeta9_J3_mn1 <- density(par_postJ3_mn1_th[,15])
zeta9_J3_mn2 <- density(par_postJ3_mn2_th[,15])

zeta10_J3_gn1 <- density(par_postJ3_gn1_th[,16])
zeta10_J3_gn2 <- density(par_postJ3_gn2_th[,16])
zeta10_J3_mn1 <- density(par_postJ3_mn1_th[,16])
zeta10_J3_mn2 <- density(par_postJ3_mn2_th[,16])

zeta11_J3_gn1 <- density(par_postJ3_gn1_th[,17])
zeta11_J3_gn2 <- density(par_postJ3_gn2_th[,17])
zeta11_J3_mn1 <- density(par_postJ3_mn1_th[,17])
zeta11_J3_mn2 <- density(par_postJ3_mn2_th[,17])

zeta12_J3_gn1 <- density(par_postJ3_gn1_th[,18])
zeta12_J3_gn2 <- density(par_postJ3_gn2_th[,18])
zeta12_J3_mn1 <- density(par_postJ3_mn1_th[,18])
zeta12_J3_mn2 <- density(par_postJ3_mn2_th[,18])

zeta13_J3_gn1 <- density(par_postJ3_gn1_th[,19])
zeta13_J3_gn2 <- density(par_postJ3_gn2_th[,19])
zeta13_J3_mn1 <- density(par_postJ3_mn1_th[,19])
zeta13_J3_mn2 <- density(par_postJ3_mn2_th[,19])

zeta14_J3_gn1 <- density(par_postJ3_gn1_th[,20])
zeta14_J3_gn2 <- density(par_postJ3_gn2_th[,20])
zeta14_J3_mn1 <- density(par_postJ3_mn1_th[,20])
zeta14_J3_mn2 <- density(par_postJ3_mn2_th[,20])

zeta15_J3_gn1 <- density(par_postJ3_gn1_th[,21])
zeta15_J3_gn2 <- density(par_postJ3_gn2_th[,21])
zeta15_J3_mn1 <- density(par_postJ3_mn1_th[,21])
zeta15_J3_mn2 <- density(par_postJ3_mn2_th[,21])

zeta16_J3_gn1 <- density(par_postJ3_gn1_th[,22])
zeta16_J3_gn2 <- density(par_postJ3_gn2_th[,22])
zeta16_J3_mn1 <- density(par_postJ3_mn1_th[,22])
zeta16_J3_mn2 <- density(par_postJ3_mn2_th[,22])

zeta17_J3_gn1 <- density(par_postJ3_gn1_th[,23])
zeta17_J3_gn2 <- density(par_postJ3_gn2_th[,23])
zeta17_J3_mn1 <- density(par_postJ3_mn1_th[,23])
zeta17_J3_mn2 <- density(par_postJ3_mn2_th[,23])

zeta18_J3_gn1 <- density(par_postJ3_gn1_th[,24])
zeta18_J3_gn2 <- density(par_postJ3_gn2_th[,24])
zeta18_J3_mn1 <- density(par_postJ3_mn1_th[,24])
zeta18_J3_mn2 <- density(par_postJ3_mn2_th[,24])

zeta19_J3_gn1 <- density(par_postJ3_gn1_th[,25])
zeta19_J3_gn2 <- density(par_postJ3_gn2_th[,25])
zeta19_J3_mn1 <- density(par_postJ3_mn1_th[,25])
zeta19_J3_mn2 <- density(par_postJ3_mn2_th[,25])

zeta20_J3_gn1 <- density(par_postJ3_gn1_th[,26])
zeta20_J3_gn2 <- density(par_postJ3_gn2_th[,26])
zeta20_J3_mn1 <- density(par_postJ3_mn1_th[,26])
zeta20_J3_mn2 <- density(par_postJ3_mn2_th[,26])

zeta21_J3_gn1 <- density(par_postJ3_gn1_th[,27])
zeta21_J3_gn2 <- density(par_postJ3_gn2_th[,27])
zeta21_J3_mn1 <- density(par_postJ3_mn1_th[,27])
zeta21_J3_mn2 <- density(par_postJ3_mn2_th[,27])

zeta22_J3_gn1 <- density(par_postJ3_gn1_th[,28])
zeta22_J3_gn2 <- density(par_postJ3_gn2_th[,28])
zeta22_J3_mn1 <- density(par_postJ3_mn1_th[,28])
zeta22_J3_mn2 <- density(par_postJ3_mn2_th[,28])

zeta23_J3_gn1 <- density(par_postJ3_gn1_th[,29])
zeta23_J3_gn2 <- density(par_postJ3_gn2_th[,29])
zeta23_J3_mn1 <- density(par_postJ3_mn1_th[,29])
zeta23_J3_mn2 <- density(par_postJ3_mn2_th[,29])

# - tau

par(mfrow=c(5,3), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, log(par_postJ3_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha", ylim=c(-3.5, -2.5)) #ylim=c(-3.2, -1), 
lines(1:100000, log(par_postJ3_gibbs_norm2[,1]), type = "l", col="red")
lines(1:100000, log(par_postJ3_metr_norm1[,1]), type = "l", col="blue")
lines(1:100000, log(par_postJ3_metr_norm2[,1]), type = "l", col="green")
#abline(h=log(phi0_p), lty=2)

plot(c(1:100000), log(par_postJ3_gn1_erg_m[,1]), type="l", ylab = "alpha", xlab="Iteration", ylim=c(-3.5,-2.5)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ3_gn2_erg_m[,1]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ3_mn1_erg_m[,1]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ3_mn2_erg_m[,1]), type="l", col="green")   #, lty=4

plot(alpha_J3_gn1, xlab = "alpha", main=" ", ylim=c(0, 8)) #, main="Density of the posterior distribution of alpha"
lines(alpha_J3_gn2, col="red")
lines(alpha_J3_mn1, col="blue")
lines(alpha_J3_mn2, col="green")
abline(v=log(mean(par_postJ3_gn1_th[,1])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_gn2_th[,1])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn1_th[,1])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn2_th[,1])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ3_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim = c(0.06, 0.15)) #, ylim=c(0.05, 0.15)
lines(1:100000, par_postJ3_gibbs_norm2[,2], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,2], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,2], type = "l", col="green")
#abline(h=beta_p, lty=2)

plot(c(1:100000), par_postJ3_gn1_erg_m[,2], type="l", ylab = "beta", xlab="Iteration", ylim=c(0.06,0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ3_gn2_erg_m[,2], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ3_mn1_erg_m[,2], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ3_mn2_erg_m[,2], type="l", col="green")   #, lty=4

plot(beta_J3_gn1, xlab = "beta", main=" ", ylim=c(0, 225)) #, ylim=c(0, 10.3), 
lines(beta_J3_gn2, col="red")
lines(beta_J3_mn1, col="blue")
lines(beta_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,2]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,2]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,2]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,2]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ3_gibbs_norm1[,3]), type = "l", xlab="Iteration", ylab = "psi1", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ3_gibbs_norm2[,3]), type = "l", col="red")
lines(1:100000, log(par_postJ3_metr_norm1[,3]), type = "l", col="blue")
lines(1:100000, log(par_postJ3_metr_norm2[,3]), type = "l", col="green")
#abline(h=log(phi1_p), lty=2)

plot(c(1:100000), log(par_postJ3_gn1_erg_m[,3]), type="l", ylab = "psi1", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ3_gn2_erg_m[,3]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ3_mn1_erg_m[,3]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ3_mn2_erg_m[,3]), type="l", col="green")   #, lty=4

plot(gamma1_J3_gn1, xlab = "psi1", main=" ", ylim=c(0, 5)) #, main="Density of the posterior distribution of alpha"
lines(gamma1_J3_gn2, col="red")
lines(gamma1_J3_mn1, col="blue")
lines(gamma1_J3_mn2, col="green")
abline(v=log(mean(par_postJ3_gn1_th[,3])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_gn2_th[,3])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn1_th[,3])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn2_th[,3])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ3_gibbs_norm1[,4]), type = "l", xlab="Iteration", ylab = "psi2", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ3_gibbs_norm2[,4]), type = "l", col="red")
lines(1:100000, log(par_postJ3_metr_norm1[,4]), type = "l", col="blue")
lines(1:100000, log(par_postJ3_metr_norm2[,4]), type = "l", col="green")
#abline(h=log(phi2_p), lty=2)

plot(c(1:100000), log(par_postJ3_gn1_erg_m[,4]), type="l", ylab = "psi2", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ3_gn2_erg_m[,4]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ3_mn1_erg_m[,4]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ3_mn2_erg_m[,4]), type="l", col="green")   #, lty=4

plot(gamma2_J3_gn1, xlab = "psi2", main=" ", ylim=c(0, 4)) #, main="Density of the posterior distribution of alpha"
lines(gamma2_J3_gn2, col="red")
lines(gamma2_J3_mn1, col="blue")
lines(gamma2_J3_mn2, col="green")
abline(v=log(mean(par_postJ3_gn1_th[,4])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_gn2_th[,4])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn1_th[,4])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn2_th[,4])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ3_gibbs_norm1[,5]), type = "l", xlab="Iteration", ylab = "psi3", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ3_gibbs_norm2[,5]), type = "l", col="red")
lines(1:100000, log(par_postJ3_metr_norm1[,5]), type = "l", col="blue")
lines(1:100000, log(par_postJ3_metr_norm2[,5]), type = "l", col="green")
#abline(h=log(phi2_p), lty=2)

plot(c(1:100000), log(par_postJ3_gn1_erg_m[,5]), type="l", ylab = "psi3", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ3_gn2_erg_m[,5]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ3_mn1_erg_m[,5]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ3_mn2_erg_m[,5]), type="l", col="green")   #, lty=4

plot(gamma3_J3_gn1, xlab = "psi3", main=" ", ylim=c(0, 3.5)) #, main="Density of the posterior distribution of alpha"
lines(gamma3_J3_gn2, col="red")
lines(gamma3_J3_mn1, col="blue")
lines(gamma3_J3_mn2, col="green")
abline(v=log(mean(par_postJ3_gn1_th[,5])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_gn2_th[,5])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn1_th[,5])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ3_mn2_th[,5])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)


#zeta

plot(1:100000, par_postJ3_gibbs_norm1[,6], type = "l", xlab="Iteration", ylab = "zeta0")
lines(1:100000, par_postJ3_gibbs_norm2[,6], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,6], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,6], type = "l", col="green")

plot(zeta0_J3_gn1, xlab = "zeta0", main=" ", ylim=c(0, 28)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta0_J3_gn2, col="red")
lines(zeta0_J3_mn1, col="blue")
lines(zeta0_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,7], type = "l", xlab="Iteration", ylab = "zeta1")
lines(1:100000, par_postJ3_gibbs_norm2[,7], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,7], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,7], type = "l", col="green")

plot(zeta1_J3_gn1, xlab = "zeta1", main=" ", ylim=c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1_J3_gn2, col="red")
lines(zeta1_J3_mn1, col="blue")
lines(zeta1_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,8], type = "l", xlab="Iteration", ylab = "zeta2")
lines(1:100000, par_postJ3_gibbs_norm2[,8], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,8], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,8], type = "l", col="green")

plot(zeta2_J3_gn1, xlab = "zeta2", main=" ", ylim=c(0, 63)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta2_J3_gn2, col="red")
lines(zeta2_J3_mn1, col="blue")
lines(zeta2_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,9], type = "l", xlab="Iteration", ylab = "zeta3")
lines(1:100000, par_postJ3_gibbs_norm2[,9], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,9], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,9], type = "l", col="green")

plot(zeta3_J3_gn1, xlab = "zeta3", main=" ", ylim=c(0, 37)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta3_J3_gn2, col="red")
lines(zeta3_J3_mn1, col="blue")
lines(zeta3_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,10], type = "l", xlab="Iteration", ylab = "zeta4")
lines(1:100000, par_postJ3_gibbs_norm2[,10], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,10], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,10], type = "l", col="green")

plot(zeta4_J3_gn1, xlab = "zeta4", main=" ", ylim=c(0, 12)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta4_J3_gn2, col="red")
lines(zeta4_J3_mn1, col="blue")
lines(zeta4_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,11], type = "l", xlab="Iteration", ylab = "zeta5")
lines(1:100000, par_postJ3_gibbs_norm2[,11], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,11], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,11], type = "l", col="green")

plot(zeta5_J3_gn1, xlab = "zeta5", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta5_J3_gn2, col="red")
lines(zeta5_J3_mn1, col="blue")
lines(zeta5_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,12], type = "l", xlab="Iteration", ylab = "zeta6")
lines(1:100000, par_postJ3_gibbs_norm2[,12], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,12], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,12], type = "l", col="green")

plot(zeta6_J3_gn1, xlab = "zeta6", main=" ", ylim=c(0, 160)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta6_J3_gn2, col="red")
lines(zeta6_J3_mn1, col="blue")
lines(zeta6_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,13], type = "l", xlab="Iteration", ylab = "zeta7")
lines(1:100000, par_postJ3_gibbs_norm2[,13], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,13], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,13], type = "l", col="green")

plot(zeta7_J3_gn1, xlab = "zeta7", main=" ", ylim=c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta7_J3_gn2, col="red")
lines(zeta7_J3_mn1, col="blue")
lines(zeta7_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,14], type = "l", xlab="Iteration", ylab = "zeta8")
lines(1:100000, par_postJ3_gibbs_norm2[,14], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,14], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,14], type = "l", col="green")

plot(zeta8_J3_gn1, xlab = "zeta8", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta8_J3_gn2, col="red")
lines(zeta8_J3_mn1, col="blue")
lines(zeta8_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,15], type = "l", xlab="Iteration", ylab = "zeta9")
lines(1:100000, par_postJ3_gibbs_norm2[,15], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,15], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,15], type = "l", col="green")

plot(zeta9_J3_gn1, xlab = "zeta9", main=" ", ylim=c(0, 200)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta9_J3_gn2, col="red")
lines(zeta9_J3_mn1, col="blue")
lines(zeta9_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,16], type = "l", xlab="Iteration", ylab = "zeta10")
lines(1:100000, par_postJ3_gibbs_norm2[,16], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,16], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,16], type = "l", col="green")

plot(zeta10_J3_gn1, xlab = "zeta10", main=" ", ylim=c(0, 17)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta10_J3_gn2, col="red")
lines(zeta10_J3_mn1, col="blue")
lines(zeta10_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,17], type = "l", xlab="Iteration", ylab = "zeta11")
lines(1:100000, par_postJ3_gibbs_norm2[,17], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,17], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,17], type = "l", col="green")

plot(zeta11_J3_gn1, xlab = "zeta11", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta11_J3_gn2, col="red")
lines(zeta11_J3_mn1, col="blue")
lines(zeta11_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,18], type = "l", xlab="Iteration", ylab = "zeta12")
lines(1:100000, par_postJ3_gibbs_norm2[,18], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,18], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,18], type = "l", col="green")

plot(zeta12_J3_gn1, xlab = "zeta12", main=" ", ylim=c(0, 210)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta12_J3_gn2, col="red")
lines(zeta12_J3_mn1, col="blue")
lines(zeta12_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,19], type = "l", xlab="Iteration", ylab = "zeta13")
lines(1:100000, par_postJ3_gibbs_norm2[,19], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,19], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,19], type = "l", col="green")

plot(zeta13_J3_gn1, xlab = "zeta13", main=" ", ylim=c(0, 210)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta13_J3_gn2, col="red")
lines(zeta13_J3_mn1, col="blue")
lines(zeta13_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,20], type = "l", xlab="Iteration", ylab = "zeta14")
lines(1:100000, par_postJ3_gibbs_norm2[,20], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,20], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,20], type = "l", col="green")

plot(zeta14_J3_gn1, xlab = "zeta14", main=" ", ylim=c(0, 275)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta14_J3_gn2, col="red")
lines(zeta14_J3_mn1, col="blue")
lines(zeta14_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,21], type = "l", xlab="Iteration", ylab = "zeta15")
lines(1:100000, par_postJ3_gibbs_norm2[,21], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,21], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,21], type = "l", col="green")

plot(zeta15_J3_gn1, xlab = "zeta15", main=" ", ylim=c(0, 205)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta15_J3_gn2, col="red")
lines(zeta15_J3_mn1, col="blue")
lines(zeta15_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,22], type = "l", xlab="Iteration", ylab = "zeta16")
lines(1:100000, par_postJ3_gibbs_norm2[,22], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,22], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,22], type = "l", col="green")

plot(zeta16_J3_gn1, xlab = "zeta16", main=" ", ylim=c(0, 31)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta16_J3_gn2, col="red")
lines(zeta16_J3_mn1, col="blue")
lines(zeta16_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,23], type = "l", xlab="Iteration", ylab = "zeta17")
lines(1:100000, par_postJ3_gibbs_norm2[,23], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,23], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,23], type = "l", col="green")

plot(zeta17_J3_gn1, xlab = "zeta17", main=" ", ylim=c(0, 135)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta17_J3_gn2, col="red")
lines(zeta17_J3_mn1, col="blue")
lines(zeta17_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,24], type = "l", xlab="Iteration", ylab = "zeta18")
lines(1:100000, par_postJ3_gibbs_norm2[,24], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,24], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,24], type = "l", col="green")

plot(zeta18_J3_gn1, xlab = "zeta18", main=" ", ylim=c(0, 190)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta18_J3_gn2, col="red")
lines(zeta18_J3_mn1, col="blue")
lines(zeta18_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,25], type = "l", xlab="Iteration", ylab = "zeta19")
lines(1:100000, par_postJ3_gibbs_norm2[,25], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,25], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,25], type = "l", col="green")

plot(zeta19_J3_gn1, xlab = "zeta19", main=" ", ylim=c(0, 30)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta19_J3_gn2, col="red")
lines(zeta19_J3_mn1, col="blue")
lines(zeta19_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,26], type = "l", xlab="Iteration", ylab = "zeta20")
lines(1:100000, par_postJ3_gibbs_norm2[,26], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,26], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,26], type = "l", col="green")

plot(zeta20_J3_gn1, xlab = "zeta17", main=" ", ylim=c(0, 85)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta20_J3_gn2, col="red")
lines(zeta20_J3_mn1, col="blue")
lines(zeta20_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,27], type = "l", xlab="Iteration", ylab = "zeta21")
lines(1:100000, par_postJ3_gibbs_norm2[,27], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,27], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,27], type = "l", col="green")

plot(zeta21_J3_gn1, xlab = "zeta21", main=" ", ylim=c(0, 205)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta21_J3_gn2, col="red")
lines(zeta21_J3_mn1, col="blue")
lines(zeta21_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,28], type = "l", xlab="Iteration", ylab = "zeta22")
lines(1:100000, par_postJ3_gibbs_norm2[,28], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,28], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,28], type = "l", col="green")

plot(zeta22_J3_gn1, xlab = "zeta22", main=" ", ylim=c(0, 38)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta22_J3_gn2, col="red")
lines(zeta22_J3_mn1, col="blue")
lines(zeta22_J3_mn2, col="green")

plot(1:100000, par_postJ3_gibbs_norm1[,29], type = "l", xlab="Iteration", ylab = "zeta23")
lines(1:100000, par_postJ3_gibbs_norm2[,29], type = "l", col="red")
lines(1:100000, par_postJ3_metr_norm1[,29], type = "l", col="blue")
lines(1:100000, par_postJ3_metr_norm2[,29], type = "l", col="green")

plot(zeta23_J3_gn1, xlab = "zeta23", main=" ", ylim=c(0, 23)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta23_J3_gn2, col="red")
lines(zeta23_J3_mn1, col="blue")
lines(zeta23_J3_mn2, col="green")


###### - zeta J3 aggregated
zeta0369_J3_gn1 <- density(par_postJ3_gn1_th[,6] + par_postJ3_gn1_th[,9] + par_postJ3_gn1_th[,12] + par_postJ3_gn1_th[,15])
zeta0369_J3_gn2 <- density(par_postJ3_gn2_th[,6] + par_postJ3_gn2_th[,9] + par_postJ3_gn2_th[,12] + par_postJ3_gn2_th[,15])
zeta0369_J3_mn1 <- density(par_postJ3_mn1_th[,6] + par_postJ3_mn1_th[,9] + par_postJ3_mn1_th[,12] + par_postJ3_mn1_th[,15])
zeta0369_J3_mn2 <- density(par_postJ3_mn2_th[,6] + par_postJ3_mn2_th[,9] + par_postJ3_mn2_th[,12] + par_postJ3_mn2_th[,15])

zeta14710_J3_gn1 <- density(par_postJ3_gn1_th[,7] + par_postJ3_gn1_th[,10] + par_postJ3_gn1_th[,13] + par_postJ3_gn1_th[,16])
zeta14710_J3_gn2 <- density(par_postJ3_gn2_th[,7] + par_postJ3_gn2_th[,10] + par_postJ3_gn2_th[,13] + par_postJ3_gn2_th[,16])
zeta14710_J3_mn1 <- density(par_postJ3_mn1_th[,7] + par_postJ3_mn1_th[,10] + par_postJ3_mn1_th[,13] + par_postJ3_mn1_th[,16])
zeta14710_J3_mn2 <- density(par_postJ3_mn2_th[,7] + par_postJ3_mn2_th[,10] + par_postJ3_mn2_th[,13] + par_postJ3_mn2_th[,16])

zeta25811_J3_gn1 <- density(par_postJ3_gn1_th[,8] + par_postJ3_gn1_th[,11] + par_postJ3_gn1_th[,14] + par_postJ3_gn1_th[,17])
zeta25811_J3_gn2 <- density(par_postJ3_gn2_th[,8] + par_postJ3_gn2_th[,11] + par_postJ3_gn2_th[,14] + par_postJ3_gn2_th[,17])
zeta25811_J3_mn1 <- density(par_postJ3_mn1_th[,8] + par_postJ3_mn1_th[,11] + par_postJ3_mn1_th[,14] + par_postJ3_mn1_th[,17])
zeta25811_J3_mn2 <- density(par_postJ3_mn2_th[,8] + par_postJ3_mn2_th[,11] + par_postJ3_mn2_th[,14] + par_postJ3_mn2_th[,17])

zeta12151821_J3_gn1 <- density(par_postJ3_gn1_th[,18] + par_postJ3_gn1_th[,21] + par_postJ3_gn1_th[,24] + par_postJ3_gn1_th[,27])
zeta12151821_J3_gn2 <- density(par_postJ3_gn2_th[,18] + par_postJ3_gn2_th[,21] + par_postJ3_gn2_th[,24] + par_postJ3_gn2_th[,27])
zeta12151821_J3_mn1 <- density(par_postJ3_mn1_th[,18] + par_postJ3_mn1_th[,21] + par_postJ3_mn1_th[,24] + par_postJ3_mn1_th[,27])
zeta12151821_J3_mn2 <- density(par_postJ3_mn2_th[,18] + par_postJ3_mn2_th[,21] + par_postJ3_mn2_th[,24] + par_postJ3_mn2_th[,27])

zeta13161922_J3_gn1 <- density(par_postJ3_gn1_th[,19] + par_postJ3_gn1_th[,22] + par_postJ3_gn1_th[,25] + par_postJ3_gn1_th[,28])
zeta13161922_J3_gn2 <- density(par_postJ3_gn2_th[,19] + par_postJ3_gn2_th[,22] + par_postJ3_gn2_th[,25] + par_postJ3_gn2_th[,28])
zeta13161922_J3_mn1 <- density(par_postJ3_mn1_th[,19] + par_postJ3_mn1_th[,22] + par_postJ3_mn1_th[,25] + par_postJ3_mn1_th[,28])
zeta13161922_J3_mn2 <- density(par_postJ3_mn2_th[,19] + par_postJ3_mn2_th[,22] + par_postJ3_mn2_th[,25] + par_postJ3_mn2_th[,28])

zeta14172023_J3_gn1 <- density(par_postJ3_gn1_th[,20] + par_postJ3_gn1_th[,23] + par_postJ3_gn1_th[,26] + par_postJ3_gn1_th[,29])
zeta14172023_J3_gn2 <- density(par_postJ3_gn2_th[,20] + par_postJ3_gn2_th[,23] + par_postJ3_gn2_th[,26] + par_postJ3_gn2_th[,29])
zeta14172023_J3_mn1 <- density(par_postJ3_mn1_th[,20] + par_postJ3_mn1_th[,23] + par_postJ3_mn1_th[,26] + par_postJ3_mn1_th[,29])
zeta14172023_J3_mn2 <- density(par_postJ3_mn2_th[,20] + par_postJ3_mn2_th[,23] + par_postJ3_mn2_th[,26] + par_postJ3_mn2_th[,29])

par(mfrow=c(6,3), mai=c(0.7,0.7, 0.2,0.2))

## - BL C0
plot(1:100000, (par_postJ3_gibbs_norm1[,6] + par_postJ3_gibbs_norm1[,9] + par_postJ3_gibbs_norm1[,12] + par_postJ3_gibbs_norm1[,15]), type = "l", xlab="Iteration", ylab = "z0+z3+z6+z9", ylim=c(0,0.4))
lines(1:100000,(par_postJ3_gibbs_norm2[,6] + par_postJ3_gibbs_norm2[,9] + par_postJ3_gibbs_norm2[,12] + par_postJ3_gibbs_norm2[,15]), type = "l", col="red")
lines(1:100000, (par_postJ3_metr_norm1[,6] + par_postJ3_metr_norm1[,9] + par_postJ3_metr_norm1[,12] + par_postJ3_metr_norm1[,15]), type = "l", col="blue")
lines(1:100000, (par_postJ3_metr_norm2[,6] + par_postJ3_metr_norm2[,9] + par_postJ3_metr_norm2[,12] + par_postJ3_metr_norm2[,15]), type = "l", col="green")

plot(c(1:100000), (par_postJ3_gn1_erg_m[,6] + par_postJ3_gn1_erg_m[,9] + par_postJ3_gn1_erg_m[,12] + par_postJ3_gn1_erg_m[,15]), type="l", ylab = "z0+z3+z6+z9", xlab="Iteration", ylim=c(0, 0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ3_gn2_erg_m[,6] + par_postJ3_gn2_erg_m[,9] + par_postJ3_gn2_erg_m[,12] + par_postJ3_gn2_erg_m[,15]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ3_mn1_erg_m[,6] + par_postJ3_mn1_erg_m[,9] + par_postJ3_mn1_erg_m[,12] + par_postJ3_mn1_erg_m[,15]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ3_mn2_erg_m[,6] + par_postJ3_mn2_erg_m[,9] + par_postJ3_mn2_erg_m[,12] + par_postJ3_mn2_erg_m[,15]), type="l", col="green")   #, lty=4

plot(zeta0369_J3_gn1, xlab = "z0+z3+z6+z9", main=" ", ylim=c(0, 45)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta0369_J3_gn2, col="red")
lines(zeta0369_J3_mn1, col="blue")
lines(zeta0369_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,6] + par_postJ3_gn1_th[,9] + par_postJ3_gn1_th[,12] + par_postJ3_gn1_th[,15]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,6] + par_postJ3_gn2_th[,9] + par_postJ3_gn2_th[,12] + par_postJ3_gn2_th[,15]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,6] + par_postJ3_mn1_th[,9] + par_postJ3_mn1_th[,12] + par_postJ3_mn1_th[,15]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,6] + par_postJ3_mn2_th[,9] + par_postJ3_mn2_th[,12] + par_postJ3_mn2_th[,15]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C0
plot(1:100000, (par_postJ3_gibbs_norm1[,18] + par_postJ3_gibbs_norm1[,21] + par_postJ3_gibbs_norm1[,24] + par_postJ3_gibbs_norm1[,27]), type = "l", xlab="Iteration", ylab = "z12+z15+z18+z21", ylim=c(0,0.2))
lines(1:100000,(par_postJ3_gibbs_norm2[,18] + par_postJ3_gibbs_norm2[,21] + par_postJ3_gibbs_norm2[,24] + par_postJ3_gibbs_norm2[,27]), type = "l", col="red")
lines(1:100000, (par_postJ3_metr_norm1[,18] + par_postJ3_metr_norm1[,21] + par_postJ3_metr_norm1[,24] + par_postJ3_metr_norm1[,27]), type = "l", col="blue")
lines(1:100000, (par_postJ3_metr_norm2[,18] + par_postJ3_metr_norm2[,21] + par_postJ3_metr_norm2[,24] + par_postJ3_metr_norm2[,27]), type = "l", col="green")

plot(c(1:100000), (par_postJ3_gn1_erg_m[,18] + par_postJ3_gn1_erg_m[,21] + par_postJ3_gn1_erg_m[,24] + par_postJ3_gn1_erg_m[,27]), type="l", ylab = "z12+z15+z18+z21", xlab="Iteration", ylim=c(0,0.2)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ3_gn2_erg_m[,18] + par_postJ3_gn2_erg_m[,21] + par_postJ3_gn2_erg_m[,24] + par_postJ3_gn2_erg_m[,27]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ3_mn1_erg_m[,18] + par_postJ3_mn1_erg_m[,21] + par_postJ3_mn1_erg_m[,24] + par_postJ3_mn1_erg_m[,27]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ3_mn2_erg_m[,18] + par_postJ3_mn2_erg_m[,21] + par_postJ3_mn2_erg_m[,24] + par_postJ3_mn2_erg_m[,27]), type="l", col="green")   #, lty=4

plot(zeta12151821_J3_gn1, xlab = "z12+z15+z18+z21", main=" ", ylim=c(0, 50)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta12151821_J3_gn2, col="red")
lines(zeta12151821_J3_mn1, col="blue")
lines(zeta12151821_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,18] + par_postJ3_gn1_th[,21] + par_postJ3_gn1_th[,24] + par_postJ3_gn1_th[,27]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,18] + par_postJ3_gn2_th[,21] + par_postJ3_gn2_th[,24] + par_postJ3_gn2_th[,27]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,18] + par_postJ3_mn1_th[,21] + par_postJ3_mn1_th[,24] + par_postJ3_mn1_th[,27]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,18] + par_postJ3_mn2_th[,21] + par_postJ3_mn2_th[,24] + par_postJ3_mn2_th[,27]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BL C1
plot(1:100000, (par_postJ3_gibbs_norm1[,7] + par_postJ3_gibbs_norm1[,10] + par_postJ3_gibbs_norm1[,13] + par_postJ3_gibbs_norm1[,16]), type = "l", xlab="Iteration", ylab = "z1+z4+z7+z10", ylim=c(0.25,0.75))
lines(1:100000,(par_postJ3_gibbs_norm2[,7] + par_postJ3_gibbs_norm2[,10] + par_postJ3_gibbs_norm2[,13] + par_postJ3_gibbs_norm2[,16]), type = "l", col="red")
lines(1:100000, (par_postJ3_metr_norm1[,7] + par_postJ3_metr_norm1[,10] + par_postJ3_metr_norm1[,13] + par_postJ3_metr_norm1[,16]), type = "l", col="blue")
lines(1:100000, (par_postJ3_metr_norm2[,7] + par_postJ3_metr_norm2[,10] + par_postJ3_metr_norm2[,13] + par_postJ3_metr_norm2[,16]), type = "l", col="green")

plot(c(1:100000), (par_postJ3_gn1_erg_m[,7] + par_postJ3_gn1_erg_m[,10] + par_postJ3_gn1_erg_m[,13] + par_postJ3_gn1_erg_m[,16]), type="l", ylab = "z1+z4+z7+z10", xlab="Iteration", ylim=c(0.25,0.75)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ3_gn2_erg_m[,7] + par_postJ3_gn2_erg_m[,10] + par_postJ3_gn2_erg_m[,13] + par_postJ3_gn2_erg_m[,16]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ3_mn1_erg_m[,7] + par_postJ3_mn1_erg_m[,10] + par_postJ3_mn1_erg_m[,13] + par_postJ3_mn1_erg_m[,16]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ3_mn2_erg_m[,7] + par_postJ3_mn2_erg_m[,10] + par_postJ3_mn2_erg_m[,13] + par_postJ3_mn2_erg_m[,16]), type="l", col="green")   #, lty=4

plot(zeta14710_J3_gn1, xlab = "z1+z4+z7+z10", main=" ", ylim=c(0, 22)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta14710_J3_gn2, col="red")
lines(zeta14710_J3_mn1, col="blue")
lines(zeta14710_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,7] + par_postJ3_gn1_th[,10] + par_postJ3_gn1_th[,13] + par_postJ3_gn1_th[,16]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,7] + par_postJ3_gn2_th[,10] + par_postJ3_gn2_th[,13] + par_postJ3_gn2_th[,16]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,7] + par_postJ3_mn1_th[,10] + par_postJ3_mn1_th[,13] + par_postJ3_mn1_th[,16]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,7] + par_postJ3_mn2_th[,10] + par_postJ3_mn2_th[,13] + par_postJ3_mn2_th[,16]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C1
plot(1:100000, (par_postJ3_gibbs_norm1[,19] + par_postJ3_gibbs_norm1[,22] + par_postJ3_gibbs_norm1[,25] + par_postJ3_gibbs_norm1[,28]), type = "l", xlab="Iteration", ylab = "z13+z16+z19+z22", ylim=c(0,0.4))
lines(1:100000,(par_postJ3_gibbs_norm2[,19] + par_postJ3_gibbs_norm2[,22] + par_postJ3_gibbs_norm2[,25] + par_postJ3_gibbs_norm2[,28]), type = "l", col="red")
lines(1:100000, (par_postJ3_metr_norm1[,19] + par_postJ3_metr_norm1[,22] + par_postJ3_metr_norm1[,25] + par_postJ3_metr_norm1[,28]), type = "l", col="blue")
lines(1:100000, (par_postJ3_metr_norm2[,19] + par_postJ3_metr_norm2[,22] + par_postJ3_metr_norm2[,25] + par_postJ3_metr_norm2[,28]), type = "l", col="green")

plot(c(1:100000), (par_postJ3_gn1_erg_m[,19] + par_postJ3_gn1_erg_m[,22] + par_postJ3_gn1_erg_m[,25] + par_postJ3_gn1_erg_m[,28]), type="l", ylab = "z13+z16+z19+z22", xlab="Iteration", ylim=c(0,0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ3_gn2_erg_m[,19] + par_postJ3_gn2_erg_m[,22] + par_postJ3_gn2_erg_m[,25] + par_postJ3_gn2_erg_m[,28]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ3_mn1_erg_m[,19] + par_postJ3_mn1_erg_m[,22] + par_postJ3_mn1_erg_m[,25] + par_postJ3_mn1_erg_m[,28]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ3_mn2_erg_m[,19] + par_postJ3_mn2_erg_m[,22] + par_postJ3_mn2_erg_m[,25] + par_postJ3_mn2_erg_m[,28]), type="l", col="green")   #, lty=4

plot(zeta13161922_J3_gn1, xlab = "z13+z16+z19+z22", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta13161922_J3_gn2, col="red")
lines(zeta13161922_J3_mn1, col="blue")
lines(zeta13161922_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,19] + par_postJ3_gn1_th[,22] + par_postJ3_gn1_th[,25] + par_postJ3_gn1_th[,28]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,19] + par_postJ3_gn2_th[,22] + par_postJ3_gn2_th[,25] + par_postJ3_gn2_th[,28]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,19] + par_postJ3_mn1_th[,22] + par_postJ3_mn1_th[,25] + par_postJ3_mn1_th[,28]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,19] + par_postJ3_mn2_th[,22] + par_postJ3_mn2_th[,25] + par_postJ3_mn2_th[,28]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BL C2
plot(1:100000, (par_postJ3_gibbs_norm1[,8] + par_postJ3_gibbs_norm1[,11] + par_postJ3_gibbs_norm1[,14] + par_postJ3_gibbs_norm1[,17]), type = "l", xlab="Iteration", ylab = "z2+z5+z8+z11", ylim=c(0,0.4))
lines(1:100000,(par_postJ3_gibbs_norm2[,8] + par_postJ3_gibbs_norm2[,11] + par_postJ3_gibbs_norm2[,14] + par_postJ3_gibbs_norm2[,17]), type = "l", col="red")
lines(1:100000, (par_postJ3_metr_norm1[,8] + par_postJ3_metr_norm1[,11] + par_postJ3_metr_norm1[,14] + par_postJ3_metr_norm1[,17]), type = "l", col="blue")
lines(1:100000, (par_postJ3_metr_norm2[,8] + par_postJ3_metr_norm2[,11] + par_postJ3_metr_norm2[,14] + par_postJ3_metr_norm2[,17]), type = "l", col="green")

plot(c(1:100000), (par_postJ3_gn1_erg_m[,8] + par_postJ3_gn1_erg_m[,11] + par_postJ3_gn1_erg_m[,14] + par_postJ3_gn1_erg_m[,17]), type="l", ylab = "z2+z5+z8+z11", xlab="Iteration", ylim=c(0,0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ3_gn2_erg_m[,8] + par_postJ3_gn2_erg_m[,11] + par_postJ3_gn2_erg_m[,14] + par_postJ3_gn2_erg_m[,17]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ3_mn1_erg_m[,8] + par_postJ3_mn1_erg_m[,11] + par_postJ3_mn1_erg_m[,14] + par_postJ3_mn1_erg_m[,17]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ3_mn2_erg_m[,8] + par_postJ3_mn2_erg_m[,11] + par_postJ3_mn2_erg_m[,14] + par_postJ3_mn2_erg_m[,17]), type="l", col="green")   #, lty=4

plot(zeta25811_J3_gn1, xlab = "z2+z5+z8+z11", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta25811_J3_gn2, col="red")
lines(zeta25811_J3_mn1, col="blue")
lines(zeta25811_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,8] + par_postJ3_gn1_th[,11] + par_postJ3_gn1_th[,14] + par_postJ3_gn1_th[,17]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,8] + par_postJ3_gn2_th[,11] + par_postJ3_gn2_th[,14] + par_postJ3_gn2_th[,17]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,8] + par_postJ3_mn1_th[,11] + par_postJ3_mn1_th[,14] + par_postJ3_mn1_th[,17]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,8] + par_postJ3_mn2_th[,11] + par_postJ3_mn2_th[,14] + par_postJ3_mn2_th[,17]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C2
plot(1:100000, (par_postJ3_gibbs_norm1[,20] + par_postJ3_gibbs_norm1[,23] + par_postJ3_gibbs_norm1[,26] + par_postJ3_gibbs_norm1[,29]), type = "l", xlab="Iteration", ylab = "z14+z17+z20+z23", ylim=c(0,0.4))
lines(1:100000,(par_postJ3_gibbs_norm2[,20] + par_postJ3_gibbs_norm2[,23] + par_postJ3_gibbs_norm2[,26] + par_postJ3_gibbs_norm2[,29]), type = "l", col="red")
lines(1:100000, (par_postJ3_metr_norm1[,20] + par_postJ3_metr_norm1[,23] + par_postJ3_metr_norm1[,26] + par_postJ3_metr_norm1[,29]), type = "l", col="blue")
lines(1:100000, (par_postJ3_metr_norm2[,20] + par_postJ3_metr_norm2[,23] + par_postJ3_metr_norm2[,26] + par_postJ3_metr_norm2[,29]), type = "l", col="green")

plot(c(1:100000), (par_postJ3_gn1_erg_m[,20] + par_postJ3_gn1_erg_m[,23] + par_postJ3_gn1_erg_m[,26] + par_postJ3_gn1_erg_m[,29]), type="l", ylab = "z14+z17+z20+z23", xlab="Iteration", ylim=c(0, 0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ3_gn2_erg_m[,20] + par_postJ3_gn2_erg_m[,23] + par_postJ3_gn2_erg_m[,26] + par_postJ3_gn2_erg_m[,29]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ3_mn1_erg_m[,20] + par_postJ3_mn1_erg_m[,23] + par_postJ3_mn1_erg_m[,26] + par_postJ3_mn1_erg_m[,29]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ3_mn2_erg_m[,20] + par_postJ3_mn2_erg_m[,23] + par_postJ3_mn2_erg_m[,26] + par_postJ3_mn2_erg_m[,29]), type="l", col="green")   #, lty=4

plot(zeta14172023_J3_gn1, xlab = "z14+z17+z20+z23", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta14172023_J3_gn2, col="red")
lines(zeta14172023_J3_mn1, col="blue")
lines(zeta14172023_J3_mn2, col="green")
abline(v=mean(par_postJ3_gn1_th[,20] + par_postJ3_gn1_th[,23] + par_postJ3_gn1_th[,26] + par_postJ3_gn1_th[,29]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_gn2_th[,20] + par_postJ3_gn2_th[,23] + par_postJ3_gn2_th[,26] + par_postJ3_gn2_th[,29]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn1_th[,20] + par_postJ3_mn1_th[,23] + par_postJ3_mn1_th[,26] + par_postJ3_mn1_th[,29]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ3_mn2_th[,20] + par_postJ3_mn2_th[,23] + par_postJ3_mn2_th[,26] + par_postJ3_mn2_th[,29]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

# - J = 4

## - Thinned sequence

par_postJ4_gn1_th <- par_postJ4_gibbs_norm1[seq(10010, 100000, 10),]
par_postJ4_gn2_th <- par_postJ4_gibbs_norm2[seq(10010, 100000, 10),]
par_postJ4_mn1_th <- par_postJ4_metr_norm1[seq(10010, 100000, 10),]
par_postJ4_mn2_th <- par_postJ4_metr_norm2[seq(10010, 100000, 10),]

#### - Posterior mean
mean_post_J4_gn1_th <- colMeans(par_postJ4_gn1_th)
mean_post_J4_gn1_th[c(1,3,4,5,6)] <- log(mean_post_J4_gn1_th[c(1,3,4,5,6)])
mean_post_J4_gn2_th <- colMeans(par_postJ4_gn2_th)
mean_post_J4_gn2_th[c(1,3,4,5,6)] <- log(mean_post_J4_gn2_th[c(1,3,4,5,6)])
mean_post_J4_mn1_th <- colMeans(par_postJ4_mn1_th)
mean_post_J4_mn1_th[c(1,3,4,5,6)] <- log(mean_post_J4_mn1_th[c(1,3,4,5,6)])
mean_post_J4_mn2_th <- colMeans(par_postJ4_mn2_th)
mean_post_J4_mn2_th[c(1,3,4,5,6)] <- log(mean_post_J4_mn2_th[c(1,3,4,5,6)])

### - tau
alpha_J4_gn1 <- density(log(par_postJ4_gn1_th[,1]))
alpha_J4_gn2 <- density(log(par_postJ4_gn2_th[,1]))
alpha_J4_mn1 <- density(log(par_postJ4_mn1_th[,1]))
alpha_J4_mn2 <- density(log(par_postJ4_mn2_th[,1]))

beta_J4_gn1 <- density(par_postJ4_gn1_th[,2])
beta_J4_gn2 <- density(par_postJ4_gn2_th[,2])
beta_J4_mn1 <- density(par_postJ4_mn1_th[,2])
beta_J4_mn2 <- density(par_postJ4_mn2_th[,2])

gamma1_J4_gn1 <- density(log(par_postJ4_gn1_th[,3]))
gamma1_J4_gn2 <- density(log(par_postJ4_gn2_th[,3]))
gamma1_J4_mn1 <- density(log(par_postJ4_mn1_th[,3]))
gamma1_J4_mn2 <- density(log(par_postJ4_mn2_th[,3]))

gamma2_J4_gn1 <- density(log(par_postJ4_gn1_th[,4]))
gamma2_J4_gn2 <- density(log(par_postJ4_gn2_th[,4]))
gamma2_J4_mn1 <- density(log(par_postJ4_mn1_th[,4]))
gamma2_J4_mn2 <- density(log(par_postJ4_mn2_th[,4]))

gamma3_J4_gn1 <- density(log(par_postJ4_gn1_th[,5]))
gamma3_J4_gn2 <- density(log(par_postJ4_gn2_th[,5]))
gamma3_J4_mn1 <- density(log(par_postJ4_mn1_th[,5]))
gamma3_J4_mn2 <- density(log(par_postJ4_mn2_th[,5]))

gamma4_J4_gn1 <- density(log(par_postJ4_gn1_th[,6]))
gamma4_J4_gn2 <- density(log(par_postJ4_gn2_th[,6]))
gamma4_J4_mn1 <- density(log(par_postJ4_mn1_th[,6]))
gamma4_J4_mn2 <- density(log(par_postJ4_mn2_th[,6]))

### - zeta

zeta0_J4_gn1 <- density(par_postJ4_gn1_th[,7])
zeta0_J4_gn2 <- density(par_postJ4_gn2_th[,7])
zeta0_J4_mn1 <- density(par_postJ4_mn1_th[,7])
zeta0_J4_mn2 <- density(par_postJ4_mn2_th[,7])

zeta1_J4_gn1 <- density(par_postJ4_gn1_th[,8])
zeta1_J4_gn2 <- density(par_postJ4_gn2_th[,8])
zeta1_J4_mn1 <- density(par_postJ4_mn1_th[,8])
zeta1_J4_mn2 <- density(par_postJ4_mn2_th[,8])

zeta2_J4_gn1 <- density(par_postJ4_gn1_th[,9])
zeta2_J4_gn2 <- density(par_postJ4_gn2_th[,9])
zeta2_J4_mn1 <- density(par_postJ4_mn1_th[,9])
zeta2_J4_mn2 <- density(par_postJ4_mn2_th[,9])

zeta3_J4_gn1 <- density(par_postJ4_gn1_th[,10])
zeta3_J4_gn2 <- density(par_postJ4_gn2_th[,10])
zeta3_J4_mn1 <- density(par_postJ4_mn1_th[,10])
zeta3_J4_mn2 <- density(par_postJ4_mn2_th[,10])

zeta4_J4_gn1 <- density(par_postJ4_gn1_th[,11])
zeta4_J4_gn2 <- density(par_postJ4_gn2_th[,11])
zeta4_J4_mn1 <- density(par_postJ4_mn1_th[,11])
zeta4_J4_mn2 <- density(par_postJ4_mn2_th[,11])

zeta5_J4_gn1 <- density(par_postJ4_gn1_th[,12])
zeta5_J4_gn2 <- density(par_postJ4_gn2_th[,12])
zeta5_J4_mn1 <- density(par_postJ4_mn1_th[,12])
zeta5_J4_mn2 <- density(par_postJ4_mn2_th[,12])

zeta6_J4_gn1 <- density(par_postJ4_gn1_th[,13])
zeta6_J4_gn2 <- density(par_postJ4_gn2_th[,13])
zeta6_J4_mn1 <- density(par_postJ4_mn1_th[,13])
zeta6_J4_mn2 <- density(par_postJ4_mn2_th[,13])

zeta7_J4_gn1 <- density(par_postJ4_gn1_th[,14])
zeta7_J4_gn2 <- density(par_postJ4_gn2_th[,14])
zeta7_J4_mn1 <- density(par_postJ4_mn1_th[,14])
zeta7_J4_mn2 <- density(par_postJ4_mn2_th[,14])

zeta8_J4_gn1 <- density(par_postJ4_gn1_th[,15])
zeta8_J4_gn2 <- density(par_postJ4_gn2_th[,15])
zeta8_J4_mn1 <- density(par_postJ4_mn1_th[,15])
zeta8_J4_mn2 <- density(par_postJ4_mn2_th[,15])

zeta9_J4_gn1 <- density(par_postJ4_gn1_th[,16])
zeta9_J4_gn2 <- density(par_postJ4_gn2_th[,16])
zeta9_J4_mn1 <- density(par_postJ4_mn1_th[,16])
zeta9_J4_mn2 <- density(par_postJ4_mn2_th[,16])

zeta10_J4_gn1 <- density(par_postJ4_gn1_th[,17])
zeta10_J4_gn2 <- density(par_postJ4_gn2_th[,17])
zeta10_J4_mn1 <- density(par_postJ4_mn1_th[,17])
zeta10_J4_mn2 <- density(par_postJ4_mn2_th[,17])

zeta11_J4_gn1 <- density(par_postJ4_gn1_th[,18])
zeta11_J4_gn2 <- density(par_postJ4_gn2_th[,18])
zeta11_J4_mn1 <- density(par_postJ4_mn1_th[,18])
zeta11_J4_mn2 <- density(par_postJ4_mn2_th[,18])

zeta12_J4_gn1 <- density(par_postJ4_gn1_th[,19])
zeta12_J4_gn2 <- density(par_postJ4_gn2_th[,19])
zeta12_J4_mn1 <- density(par_postJ4_mn1_th[,19])
zeta12_J4_mn2 <- density(par_postJ4_mn2_th[,19])

zeta13_J4_gn1 <- density(par_postJ4_gn1_th[,20])
zeta13_J4_gn2 <- density(par_postJ4_gn2_th[,20])
zeta13_J4_mn1 <- density(par_postJ4_mn1_th[,20])
zeta13_J4_mn2 <- density(par_postJ4_mn2_th[,20])

zeta14_J4_gn1 <- density(par_postJ4_gn1_th[,21])
zeta14_J4_gn2 <- density(par_postJ4_gn2_th[,21])
zeta14_J4_mn1 <- density(par_postJ4_mn1_th[,21])
zeta14_J4_mn2 <- density(par_postJ4_mn2_th[,21])

zeta15_J4_gn1 <- density(par_postJ4_gn1_th[,22])
zeta15_J4_gn2 <- density(par_postJ4_gn2_th[,22])
zeta15_J4_mn1 <- density(par_postJ4_mn1_th[,22])
zeta15_J4_mn2 <- density(par_postJ4_mn2_th[,22])

zeta16_J4_gn1 <- density(par_postJ4_gn1_th[,23])
zeta16_J4_gn2 <- density(par_postJ4_gn2_th[,23])
zeta16_J4_mn1 <- density(par_postJ4_mn1_th[,23])
zeta16_J4_mn2 <- density(par_postJ4_mn2_th[,23])

zeta17_J4_gn1 <- density(par_postJ4_gn1_th[,24])
zeta17_J4_gn2 <- density(par_postJ4_gn2_th[,24])
zeta17_J4_mn1 <- density(par_postJ4_mn1_th[,24])
zeta17_J4_mn2 <- density(par_postJ4_mn2_th[,24])

zeta18_J4_gn1 <- density(par_postJ4_gn1_th[,25])
zeta18_J4_gn2 <- density(par_postJ4_gn2_th[,25])
zeta18_J4_mn1 <- density(par_postJ4_mn1_th[,25])
zeta18_J4_mn2 <- density(par_postJ4_mn2_th[,25])

zeta19_J4_gn1 <- density(par_postJ4_gn1_th[,26])
zeta19_J4_gn2 <- density(par_postJ4_gn2_th[,26])
zeta19_J4_mn1 <- density(par_postJ4_mn1_th[,26])
zeta19_J4_mn2 <- density(par_postJ4_mn2_th[,26])

zeta20_J4_gn1 <- density(par_postJ4_gn1_th[,27])
zeta20_J4_gn2 <- density(par_postJ4_gn2_th[,27])
zeta20_J4_mn1 <- density(par_postJ4_mn1_th[,27])
zeta20_J4_mn2 <- density(par_postJ4_mn2_th[,27])

zeta21_J4_gn1 <- density(par_postJ4_gn1_th[,28])
zeta21_J4_gn2 <- density(par_postJ4_gn2_th[,28])
zeta21_J4_mn1 <- density(par_postJ4_mn1_th[,28])
zeta21_J4_mn2 <- density(par_postJ4_mn2_th[,28])

zeta22_J4_gn1 <- density(par_postJ4_gn1_th[,29])
zeta22_J4_gn2 <- density(par_postJ4_gn2_th[,29])
zeta22_J4_mn1 <- density(par_postJ4_mn1_th[,29])
zeta22_J4_mn2 <- density(par_postJ4_mn2_th[,29])

zeta23_J4_gn1 <- density(par_postJ4_gn1_th[,30])
zeta23_J4_gn2 <- density(par_postJ4_gn2_th[,30])
zeta23_J4_mn1 <- density(par_postJ4_mn1_th[,30])
zeta23_J4_mn2 <- density(par_postJ4_mn2_th[,30])

zeta24_J4_gn1 <- density(par_postJ4_gn1_th[,31])
zeta24_J4_gn2 <- density(par_postJ4_gn2_th[,31])
zeta24_J4_mn1 <- density(par_postJ4_mn1_th[,31])
zeta24_J4_mn2 <- density(par_postJ4_mn2_th[,31])

zeta25_J4_gn1 <- density(par_postJ4_gn1_th[,32])
zeta25_J4_gn2 <- density(par_postJ4_gn2_th[,32])
zeta25_J4_mn1 <- density(par_postJ4_mn1_th[,32])
zeta25_J4_mn2 <- density(par_postJ4_mn2_th[,32])

zeta26_J4_gn1 <- density(par_postJ4_gn1_th[,33])
zeta26_J4_gn2 <- density(par_postJ4_gn2_th[,33])
zeta26_J4_mn1 <- density(par_postJ4_mn1_th[,33])
zeta26_J4_mn2 <- density(par_postJ4_mn2_th[,33])

zeta27_J4_gn1 <- density(par_postJ4_gn1_th[,34])
zeta27_J4_gn2 <- density(par_postJ4_gn2_th[,34])
zeta27_J4_mn1 <- density(par_postJ4_mn1_th[,34])
zeta27_J4_mn2 <- density(par_postJ4_mn2_th[,34])

zeta28_J4_gn1 <- density(par_postJ4_gn1_th[,35])
zeta28_J4_gn2 <- density(par_postJ4_gn2_th[,35])
zeta28_J4_mn1 <- density(par_postJ4_mn1_th[,35])
zeta28_J4_mn2 <- density(par_postJ4_mn2_th[,35])

zeta29_J4_gn1 <- density(par_postJ4_gn1_th[,36])
zeta29_J4_gn2 <- density(par_postJ4_gn2_th[,36])
zeta29_J4_mn1 <- density(par_postJ4_mn1_th[,36])
zeta29_J4_mn2 <- density(par_postJ4_mn2_th[,36])

par(mfrow=c(6,3), mai=c(0.7,0.7, 0.2,0.2))

#par(mfrow=c(1,1), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, log(par_postJ4_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha", ylim=c(-3.5, -2.5)) #ylim=c(-3.2, -1), 
lines(1:100000, log(par_postJ4_gibbs_norm2[,1]), type = "l", col="red")
lines(1:100000, log(par_postJ4_metr_norm1[,1]), type = "l", col="blue")
lines(1:100000, log(par_postJ4_metr_norm2[,1]), type = "l", col="green")
#abline(h=log(phi0_p), lty=2)

plot(c(1:100000), log(par_postJ4_gn1_erg_m[,1]), type="l", ylab = "alpha", xlab="Iteration", ylim=c(-3.5,-2.5)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ4_gn2_erg_m[,1]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ4_mn1_erg_m[,1]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ4_mn2_erg_m[,1]), type="l", col="green")   #, lty=4

plot(alpha_J4_gn1, xlab = "alpha", main=" ", ylim=c(0, 8)) #, main="Density of the posterior distribution of alpha"
lines(alpha_J4_gn2, col="red")
lines(alpha_J4_mn1, col="blue")
lines(alpha_J4_mn2, col="green")
abline(v=log(mean(par_postJ4_gn1_th[,1])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_gn2_th[,1])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn1_th[,1])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn2_th[,1])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim = c(0.06, 0.15)) #, ylim=c(0.05, 0.15)
lines(1:100000, par_postJ4_gibbs_norm2[,2], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,2], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,2], type = "l", col="green")
#abline(h=beta_p, lty=2)

plot(c(1:100000), par_postJ4_gn1_erg_m[,2], type="l", ylab = "beta", xlab="Iteration", ylim=c(0.06,0.15)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ4_gn2_erg_m[,2], type="l", col="red") #, lty=2
lines(c(1:100000), par_postJ4_mn1_erg_m[,2], type="l", col="blue")   #, lty=3
lines(c(1:100000), par_postJ4_mn2_erg_m[,2], type="l", col="green")   #, lty=4

plot(beta_J4_gn1, xlab = "beta", main=" ", ylim=c(0, 230)) #, ylim=c(0, 10.3), 
lines(beta_J4_gn2, col="red")
lines(beta_J4_mn1, col="blue")
lines(beta_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,2]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,2]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,2]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,2]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ4_gibbs_norm1[,3]), type = "l", xlab="Iteration", ylab = "psi1", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ4_gibbs_norm2[,3]), type = "l", col="red")
lines(1:100000, log(par_postJ4_metr_norm1[,3]), type = "l", col="blue")
lines(1:100000, log(par_postJ4_metr_norm2[,3]), type = "l", col="green")
#abline(h=log(phi1_p), lty=2)

plot(c(1:100000), log(par_postJ4_gn1_erg_m[,3]), type="l", ylab = "psi1", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ4_gn2_erg_m[,3]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ4_mn1_erg_m[,3]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ4_mn2_erg_m[,3]), type="l", col="green")   #, lty=4

plot(gamma1_J4_gn1, xlab = "psi1", main=" ", ylim=c(0, 4)) #, main="Density of the posterior distribution of alpha"
lines(gamma1_J4_gn2, col="red")
lines(gamma1_J4_mn1, col="blue")
lines(gamma1_J4_mn2, col="green")
abline(v=log(mean(par_postJ4_gn1_th[,3])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_gn2_th[,3])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn1_th[,3])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn2_th[,3])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ4_gibbs_norm1[,4]), type = "l", xlab="Iteration", ylab = "psi2", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ4_gibbs_norm2[,4]), type = "l", col="red")
lines(1:100000, log(par_postJ4_metr_norm1[,4]), type = "l", col="blue")
lines(1:100000, log(par_postJ4_metr_norm2[,4]), type = "l", col="green")
#abline(h=log(phi2_p), lty=2)

plot(c(1:100000), log(par_postJ4_gn1_erg_m[,4]), type="l", ylab = "psi2", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ4_gn2_erg_m[,4]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ4_mn1_erg_m[,4]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ4_mn2_erg_m[,4]), type="l", col="green")   #, lty=4

plot(gamma2_J4_gn1, xlab = "psi2", main=" ", ylim=c(0, 240), xlim=c(-0.7,0)) #, main="Density of the posterior distribution of alpha"
lines(gamma2_J4_gn2, col="red")
lines(gamma2_J4_mn1, col="blue")
lines(gamma2_J4_mn2, col="green")
abline(v=log(mean(par_postJ4_gn1_th[,4])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_gn2_th[,4])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn1_th[,4])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn2_th[,4])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ4_gibbs_norm1[,5]), type = "l", xlab="Iteration", ylab = "psi3", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ4_gibbs_norm2[,5]), type = "l", col="red")
lines(1:100000, log(par_postJ4_metr_norm1[,5]), type = "l", col="blue")
lines(1:100000, log(par_postJ4_metr_norm2[,5]), type = "l", col="green")
#abline(h=log(phi2_p), lty=2)

plot(c(1:100000), log(par_postJ4_gn1_erg_m[,5]), type="l", ylab = "psi3", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ4_gn2_erg_m[,5]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ4_mn1_erg_m[,5]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ4_mn2_erg_m[,5]), type="l", col="green")   #, lty=4

plot(gamma3_J4_gn1, xlab = "psi3", main=" ", ylim=c(0, 4)) #, main="Density of the posterior distribution of alpha"
lines(gamma3_J4_gn2, col="red")
lines(gamma3_J4_mn1, col="blue")
lines(gamma3_J4_mn2, col="green")
abline(v=log(mean(par_postJ4_gn1_th[,5])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_gn2_th[,5])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn1_th[,5])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn2_th[,5])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, log(par_postJ4_gibbs_norm1[,6]), type = "l", xlab="Iteration", ylab = "psi4", ylim=c(-1, 0)) #
lines(1:100000, log(par_postJ4_gibbs_norm2[,6]), type = "l", col="red")
lines(1:100000, log(par_postJ4_metr_norm1[,6]), type = "l", col="blue")
lines(1:100000, log(par_postJ4_metr_norm2[,6]), type = "l", col="green")
#abline(h=log(phi2_p), lty=2)

plot(c(1:100000), log(par_postJ4_gn1_erg_m[,6]), type="l", ylab = "psi4", xlab="Iteration", ylim=c(-1, 0)) #, ylim=c(0,1)
lines(c(1:100000), log(par_postJ4_gn2_erg_m[,6]), type="l", col="red") #, lty=2
lines(c(1:100000), log(par_postJ4_mn1_erg_m[,6]), type="l", col="blue")   #, lty=3
lines(c(1:100000), log(par_postJ4_mn2_erg_m[,6]), type="l", col="green")   #, lty=4

plot(gamma4_J4_gn1, xlab = "psi4", main=" ", ylim=c(0, 3)) #, main="Density of the posterior distribution of alpha"
lines(gamma4_J4_gn2, col="red")
lines(gamma4_J4_mn1, col="blue")
lines(gamma4_J4_mn2, col="green")
abline(v=log(mean(par_postJ4_gn1_th[,6])), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_gn2_th[,6])), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn1_th[,6])), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=log(mean(par_postJ4_mn2_th[,6])), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

#zeta

par(mfrow=c(5,3), mai=c(0.7,0.7, 0.2,0.2))

plot(1:100000, par_postJ4_gibbs_norm1[,7], type = "l", xlab="Iteration", ylab = "zeta0")
lines(1:100000, par_postJ4_gibbs_norm2[,7], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,7], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,7], type = "l", col="green")

plot(zeta0_J4_gn1, xlab = "zeta0", main=" ", ylim=c(0, 22)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta0_J4_gn2, col="red")
lines(zeta0_J4_mn1, col="blue")
lines(zeta0_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,7]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,7]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,7]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,7]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,8], type = "l", xlab="Iteration", ylab = "zeta1")
lines(1:100000, par_postJ4_gibbs_norm2[,8], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,8], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,8], type = "l", col="green")

plot(zeta1_J4_gn1, xlab = "zeta1", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1_J4_gn2, col="red")
lines(zeta1_J4_mn1, col="blue")
lines(zeta1_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,8]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,8]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,8]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,8]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,9], type = "l", xlab="Iteration", ylab = "zeta2")
lines(1:100000, par_postJ4_gibbs_norm2[,9], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,9], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,9], type = "l", col="green")

plot(zeta2_J4_gn1, xlab = "zeta2", main=" ", ylim=c(0, 190)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta2_J4_gn2, col="red")
lines(zeta2_J4_mn1, col="blue")
lines(zeta2_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,9]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,9]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,9]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,9]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,10], type = "l", xlab="Iteration", ylab = "zeta3")
lines(1:100000, par_postJ4_gibbs_norm2[,10], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,10], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,10], type = "l", col="green")

plot(zeta3_J4_gn1, xlab = "zeta3", main=" ", ylim=c(0, 32)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta3_J4_gn2, col="red")
lines(zeta3_J4_mn1, col="blue")
lines(zeta3_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,10]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,10]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,10]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,10]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,11], type = "l", xlab="Iteration", ylab = "zeta4")
lines(1:100000, par_postJ4_gibbs_norm2[,11], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,11], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,11], type = "l", col="green")

plot(zeta4_J4_gn1, xlab = "zeta4", main=" ", ylim=c(0, 14)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta4_J4_gn2, col="red")
lines(zeta4_J4_mn1, col="blue")
lines(zeta4_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,11]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,11]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,11]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,11]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,12], type = "l", xlab="Iteration", ylab = "zeta5")
lines(1:100000, par_postJ4_gibbs_norm2[,12], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,12], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,12], type = "l", col="green")

plot(zeta5_J4_gn1, xlab = "zeta5", main=" ", ylim=c(0, 52)) #, ylim=c(0, 25), main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta5_J4_gn2, col="red")
lines(zeta5_J4_mn1, col="blue")
lines(zeta5_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,12]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,12]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,12]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,12]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,13], type = "l", xlab="Iteration", ylab = "zeta6")
lines(1:100000, par_postJ4_gibbs_norm2[,13], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,13], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,13], type = "l", col="green")

plot(zeta6_J4_gn1, xlab = "zeta6", main=" ", ylim=c(0, 35)) #, ylim=c(0, 160), main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta6_J4_gn2, col="red")
lines(zeta6_J4_mn1, col="blue")
lines(zeta6_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,13]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,13]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,13]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,13]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,14], type = "l", xlab="Iteration", ylab = "zeta7")
lines(1:100000, par_postJ4_gibbs_norm2[,14], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,14], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,14], type = "l", col="green")

plot(zeta7_J4_gn1, xlab = "zeta7", main=" ", ylim=c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta7_J4_gn2, col="red")
lines(zeta7_J4_mn1, col="blue")
lines(zeta7_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,14]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,14]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,14]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,14]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,15], type = "l", xlab="Iteration", ylab = "zeta8")
lines(1:100000, par_postJ4_gibbs_norm2[,15], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,15], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,15], type = "l", col="green")

plot(zeta8_J4_gn1, xlab = "zeta8", main=" ", ylim=c(0, 27)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta8_J4_gn2, col="red")
lines(zeta8_J4_mn1, col="blue")
lines(zeta8_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,15]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,15]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,15]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,15]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,16], type = "l", xlab="Iteration", ylab = "zeta9")
lines(1:100000, par_postJ4_gibbs_norm2[,16], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,16], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,16], type = "l", col="green")

plot(zeta9_J4_gn1, xlab = "zeta9", main=" ", ylim=c(0, 350)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta9_J4_gn2, col="red")
lines(zeta9_J4_mn1, col="blue")
lines(zeta9_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,16]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,16]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,16]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,16]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,17], type = "l", xlab="Iteration", ylab = "zeta10")
lines(1:100000, par_postJ4_gibbs_norm2[,17], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,17], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,17], type = "l", col="green")

plot(zeta10_J4_gn1, xlab = "zeta10", main=" ", ylim=c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta10_J4_gn2, col="red")
lines(zeta10_J4_mn1, col="blue")
lines(zeta10_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,17]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,17]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,17]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,17]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,18], type = "l", xlab="Iteration", ylab = "zeta11")
lines(1:100000, par_postJ4_gibbs_norm2[,18], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,18], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,18], type = "l", col="green")

plot(zeta11_J4_gn1, xlab = "zeta11", main=" ", ylim=c(0, 22)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta11_J4_gn2, col="red")
lines(zeta11_J4_mn1, col="blue")
lines(zeta11_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,18]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,18]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,18]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,18]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,19], type = "l", xlab="Iteration", ylab = "zeta12")
lines(1:100000, par_postJ4_gibbs_norm2[,19], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,19], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,19], type = "l", col="green")

plot(zeta12_J4_gn1, xlab = "zeta12", main=" ", ylim=c(0, 350)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta12_J4_gn2, col="red")
lines(zeta12_J4_mn1, col="blue")
lines(zeta12_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,19]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,19]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,19]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,19]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,20], type = "l", xlab="Iteration", ylab = "zeta13")
lines(1:100000, par_postJ4_gibbs_norm2[,20], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,20], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,20], type = "l", col="green")

plot(zeta13_J4_gn1, xlab = "zeta13", main=" ", ylim=c(0, 25)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta13_J4_gn2, col="red")
lines(zeta13_J4_mn1, col="blue")
lines(zeta13_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,20]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,20]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,20]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,20]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,21], type = "l", xlab="Iteration", ylab = "zeta14")
lines(1:100000, par_postJ4_gibbs_norm2[,21], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,21], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,21], type = "l", col="green")

plot(zeta14_J4_gn1, xlab = "zeta14", main=" ", ylim=c(0, 5000), xlim=c(0,0.003)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta14_J4_gn2, col="red")
lines(zeta14_J4_mn1, col="blue")
lines(zeta14_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,21]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,21]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,21]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,21]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,22], type = "l", xlab="Iteration", ylab = "zeta15")
lines(1:100000, par_postJ4_gibbs_norm2[,22], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,22], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,22], type = "l", col="green")

plot(zeta15_J4_gn1, xlab = "zeta15", main=" ", ylim=c(0, 160)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta15_J4_gn2, col="red")
lines(zeta15_J4_mn1, col="blue")
lines(zeta15_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,22]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,22]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,22]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,22]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,23], type = "l", xlab="Iteration", ylab = "zeta16")
lines(1:100000, par_postJ4_gibbs_norm2[,23], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,23], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,23], type = "l", col="green")

plot(zeta16_J4_gn1, xlab = "zeta16", main=" ", ylim=c(0, 80)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta16_J4_gn2, col="red")
lines(zeta16_J4_mn1, col="blue")
lines(zeta16_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,23]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,23]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,23]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,23]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,24], type = "l", xlab="Iteration", ylab = "zeta17")
lines(1:100000, par_postJ4_gibbs_norm2[,24], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,24], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,24], type = "l", col="green")

plot(zeta17_J4_gn1, xlab = "zeta17", main=" ", ylim=c(0, 720)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta17_J4_gn2, col="red")
lines(zeta17_J4_mn1, col="blue")
lines(zeta17_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,24]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,24]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,24]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,24]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,25], type = "l", xlab="Iteration", ylab = "zeta18")
lines(1:100000, par_postJ4_gibbs_norm2[,25], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,25], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,25], type = "l", col="green")

plot(zeta18_J4_gn1, xlab = "zeta18", main=" ", ylim=c(0, 150)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta18_J4_gn2, col="red")
lines(zeta18_J4_mn1, col="blue")
lines(zeta18_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,25]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,25]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,25]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,25]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,26], type = "l", xlab="Iteration", ylab = "zeta19")
lines(1:100000, par_postJ4_gibbs_norm2[,26], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,26], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,26], type = "l", col="green")

plot(zeta19_J4_gn1, xlab = "zeta19", main=" ", ylim=c(0, 55)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta19_J4_gn2, col="red")
lines(zeta19_J4_mn1, col="blue")
lines(zeta19_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,26]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,26]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,26]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,26]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,27], type = "l", xlab="Iteration", ylab = "zeta20")
lines(1:100000, par_postJ4_gibbs_norm2[,27], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,27], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,27], type = "l", col="green")

plot(zeta20_J4_gn1, xlab = "zeta20", main=" ", ylim=c(0, 480)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta20_J4_gn2, col="red")
lines(zeta20_J4_mn1, col="blue")
lines(zeta20_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,27]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,27]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,27]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,27]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,28], type = "l", xlab="Iteration", ylab = "zeta21")
lines(1:100000, par_postJ4_gibbs_norm2[,28], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,28], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,28], type = "l", col="green")

plot(zeta21_J4_gn1, xlab = "zeta21", main=" ", ylim=c(0, 150)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta21_J4_gn2, col="red")
lines(zeta21_J4_mn1, col="blue")
lines(zeta21_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,28]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,28]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,28]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,28]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,29], type = "l", xlab="Iteration", ylab = "zeta22")
lines(1:100000, par_postJ4_gibbs_norm2[,29], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,29], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,29], type = "l", col="green")

plot(zeta22_J4_gn1, xlab = "zeta22", main=" ", ylim=c(0, 35)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta22_J4_gn2, col="red")
lines(zeta22_J4_mn1, col="blue")
lines(zeta22_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,29]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,29]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,29]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,29]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,30], type = "l", xlab="Iteration", ylab = "zeta23")
lines(1:100000, par_postJ4_gibbs_norm2[,30], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,30], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,30], type = "l", col="green")

plot(zeta23_J4_gn1, xlab = "zeta23", main=" ", ylim=c(0, 75)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta23_J4_gn2, col="red")
lines(zeta23_J4_mn1, col="blue")
lines(zeta23_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,30]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,30]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,30]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,30]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,31], type = "l", xlab="Iteration", ylab = "zeta24")
lines(1:100000, par_postJ4_gibbs_norm2[,31], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,31], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,31], type = "l", col="green")

plot(zeta24_J4_gn1, xlab = "zeta24", main=" ", ylim=c(0, 400)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta24_J4_gn2, col="red")
lines(zeta24_J4_mn1, col="blue")
lines(zeta24_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,31]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,31]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,31]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,31]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,32], type = "l", xlab="Iteration", ylab = "zeta25")
lines(1:100000, par_postJ4_gibbs_norm2[,32], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,32], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,32], type = "l", col="green")

plot(zeta25_J4_gn1, xlab = "zeta25", main=" ", ylim=c(0, 35)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta25_J4_gn2, col="red")
lines(zeta25_J4_mn1, col="blue")
lines(zeta25_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,32]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,32]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,32]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,32]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,33], type = "l", xlab="Iteration", ylab = "zeta26")
lines(1:100000, par_postJ4_gibbs_norm2[,33], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,33], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,33], type = "l", col="green")

plot(zeta26_J4_gn1, xlab = "zeta26", main=" ", ylim=c(0, 30)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta26_J4_gn2, col="red")
lines(zeta26_J4_mn1, col="blue")
lines(zeta26_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,33]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,33]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,33]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,33]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,34], type = "l", xlab="Iteration", ylab = "zeta27")
lines(1:100000, par_postJ4_gibbs_norm2[,34], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,34], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,34], type = "l", col="green")

plot(zeta27_J4_gn1, xlab = "zeta27", main=" ", ylim=c(0, 400)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta27_J4_gn2, col="red")
lines(zeta27_J4_mn1, col="blue")
lines(zeta27_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,34]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,34]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,34]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,34]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,35], type = "l", xlab="Iteration", ylab = "zeta28")
lines(1:100000, par_postJ4_gibbs_norm2[,35], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,35], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,35], type = "l", col="green")

plot(zeta28_J4_gn1, xlab = "zeta28", main=" ", ylim=c(0, 35)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta28_J4_gn2, col="red")
lines(zeta28_J4_mn1, col="blue")
lines(zeta28_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,35]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,35]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,35]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,35]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

plot(1:100000, par_postJ4_gibbs_norm1[,36], type = "l", xlab="Iteration", ylab = "zeta29")
lines(1:100000, par_postJ4_gibbs_norm2[,36], type = "l", col="red")
lines(1:100000, par_postJ4_metr_norm1[,36], type = "l", col="blue")
lines(1:100000, par_postJ4_metr_norm2[,36], type = "l", col="green")

plot(zeta29_J4_gn1, xlab = "zeta29", main=" ", ylim=c(0, 15)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta29_J4_gn2, col="red")
lines(zeta29_J4_mn1, col="blue")
lines(zeta29_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,36]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,36]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,36]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,36]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

###### - zeta J4 aggregated

zeta036912_J4_gn1 <- density(par_postJ4_gn1_th[,7] + par_postJ4_gn1_th[,10] + par_postJ4_gn1_th[,13] + par_postJ4_gn1_th[,16] + par_postJ4_gn1_th[,19])
zeta036912_J4_gn2 <- density(par_postJ4_gn2_th[,7] + par_postJ4_gn2_th[,10] + par_postJ4_gn2_th[,13] + par_postJ4_gn2_th[,16] + par_postJ4_gn2_th[,19])
zeta036912_J4_mn1 <- density(par_postJ4_mn1_th[,7] + par_postJ4_mn1_th[,10] + par_postJ4_mn1_th[,13] + par_postJ4_mn1_th[,16] + par_postJ4_mn1_th[,19])
zeta036912_J4_mn2 <- density(par_postJ4_mn2_th[,7] + par_postJ4_mn2_th[,10] + par_postJ4_mn2_th[,13] + par_postJ4_mn2_th[,16] + par_postJ4_mn2_th[,19])

zeta1471013_J4_gn1 <- density(par_postJ4_gn1_th[,8] + par_postJ4_gn1_th[,11] + par_postJ4_gn1_th[,14] + par_postJ4_gn1_th[,17] + par_postJ4_gn1_th[,20])
zeta1471013_J4_gn2 <- density(par_postJ4_gn2_th[,8] + par_postJ4_gn2_th[,11] + par_postJ4_gn2_th[,14] + par_postJ4_gn2_th[,17] + par_postJ4_gn2_th[,20])
zeta1471013_J4_mn1 <- density(par_postJ4_mn1_th[,8] + par_postJ4_mn1_th[,11] + par_postJ4_mn1_th[,14] + par_postJ4_mn1_th[,17] + par_postJ4_mn1_th[,20])
zeta1471013_J4_mn2 <- density(par_postJ4_mn2_th[,8] + par_postJ4_mn2_th[,11] + par_postJ4_mn2_th[,14] + par_postJ4_mn2_th[,17] + par_postJ4_mn2_th[,20])

zeta2581114_J4_gn1 <- density(par_postJ4_gn1_th[,9] + par_postJ4_gn1_th[,12] + par_postJ4_gn1_th[,15] + par_postJ4_gn1_th[,18] + par_postJ4_gn1_th[,21])
zeta2581114_J4_gn2 <- density(par_postJ4_gn2_th[,9] + par_postJ4_gn2_th[,12] + par_postJ4_gn2_th[,15] + par_postJ4_gn2_th[,18] + par_postJ4_gn2_th[,21])
zeta2581114_J4_mn1 <- density(par_postJ4_mn1_th[,9] + par_postJ4_mn1_th[,12] + par_postJ4_mn1_th[,15] + par_postJ4_mn1_th[,18] + par_postJ4_mn1_th[,21])
zeta2581114_J4_mn2 <- density(par_postJ4_mn2_th[,9] + par_postJ4_mn2_th[,12] + par_postJ4_mn2_th[,15] + par_postJ4_mn2_th[,18] + par_postJ4_mn2_th[,21])

zeta1518212427_J4_gn1 <- density(par_postJ4_gn1_th[,22] + par_postJ4_gn1_th[,25] + par_postJ4_gn1_th[,28] + par_postJ4_gn1_th[,31] + par_postJ4_gn1_th[,34])
zeta1518212427_J4_gn2 <- density(par_postJ4_gn2_th[,22] + par_postJ4_gn2_th[,25] + par_postJ4_gn2_th[,28] + par_postJ4_gn2_th[,31] + par_postJ4_gn2_th[,34])
zeta1518212427_J4_mn1 <- density(par_postJ4_mn1_th[,22] + par_postJ4_mn1_th[,25] + par_postJ4_mn1_th[,28] + par_postJ4_mn1_th[,31] + par_postJ4_mn1_th[,34])
zeta1518212427_J4_mn2 <- density(par_postJ4_mn2_th[,22] + par_postJ4_mn2_th[,25] + par_postJ4_mn2_th[,28] + par_postJ4_mn2_th[,31] + par_postJ4_mn2_th[,34])

zeta1619222528_J4_gn1 <- density(par_postJ4_gn1_th[,23] + par_postJ4_gn1_th[,26] + par_postJ4_gn1_th[,29] + par_postJ4_gn1_th[,32] + par_postJ4_gn1_th[,35])
zeta1619222528_J4_gn2 <- density(par_postJ4_gn2_th[,23] + par_postJ4_gn2_th[,26] + par_postJ4_gn2_th[,29] + par_postJ4_gn2_th[,32] + par_postJ4_gn2_th[,35])
zeta1619222528_J4_mn1 <- density(par_postJ4_mn1_th[,23] + par_postJ4_mn1_th[,26] + par_postJ4_mn1_th[,29] + par_postJ4_mn1_th[,32] + par_postJ4_mn1_th[,35])
zeta1619222528_J4_mn2 <- density(par_postJ4_mn1_th[,23] + par_postJ4_mn1_th[,26] + par_postJ4_mn1_th[,29] + par_postJ4_mn1_th[,32] + par_postJ4_mn1_th[,35])

zeta1720232629_J4_gn1 <- density(par_postJ4_gn1_th[,24] + par_postJ4_gn1_th[,27] + par_postJ4_gn1_th[,30] + par_postJ4_gn1_th[,33] + par_postJ4_gn1_th[,36])
zeta1720232629_J4_gn2 <- density(par_postJ4_gn2_th[,24] + par_postJ4_gn2_th[,27] + par_postJ4_gn2_th[,30] + par_postJ4_gn2_th[,33] + par_postJ4_gn2_th[,36])
zeta1720232629_J4_mn1 <- density(par_postJ4_mn1_th[,24] + par_postJ4_mn1_th[,27] + par_postJ4_mn1_th[,30] + par_postJ4_mn1_th[,33] + par_postJ4_mn1_th[,36])
zeta1720232629_J4_mn2 <- density(par_postJ4_mn2_th[,24] + par_postJ4_mn2_th[,27] + par_postJ4_mn2_th[,30] + par_postJ4_mn2_th[,33] + par_postJ4_mn2_th[,36])

par(mfrow=c(6,3), mai=c(0.7,0.7, 0.2,0.2))

## - BL C0
plot(1:100000, (par_postJ4_gibbs_norm1[,7] + par_postJ4_gibbs_norm1[,10] + par_postJ4_gibbs_norm1[,13] + par_postJ4_gibbs_norm1[,16] + par_postJ4_gibbs_norm1[,19]), type = "l", xlab="Iteration", ylab = "z0+z3+z6+z9+z12", ylim=c(0,0.4))
lines(1:100000,(par_postJ4_gibbs_norm2[,7] + par_postJ4_gibbs_norm2[,10] + par_postJ4_gibbs_norm2[,13] + par_postJ4_gibbs_norm2[,16] + par_postJ4_gibbs_norm2[,19]), type = "l", col="red")
lines(1:100000, (par_postJ4_metr_norm1[,7] + par_postJ4_metr_norm1[,10] + par_postJ4_metr_norm1[,13] + par_postJ4_metr_norm1[,16] + par_postJ4_metr_norm1[,19]), type = "l", col="blue")
lines(1:100000, (par_postJ4_metr_norm2[,7] + par_postJ4_metr_norm2[,10] + par_postJ4_metr_norm2[,13] + par_postJ4_metr_norm2[,16] + par_postJ4_metr_norm2[,19]), type = "l", col="green")

plot(c(1:100000), (par_postJ4_gn1_erg_m[,7] + par_postJ4_gn1_erg_m[,10] + par_postJ4_gn1_erg_m[,13] + par_postJ4_gn1_erg_m[,16] + par_postJ4_gn1_erg_m[,19]), type="l", ylab = "z0+z3+z6+z9+z12", xlab="Iteration", ylim=c(0, 0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ4_gn2_erg_m[,7] + par_postJ4_gn2_erg_m[,10] + par_postJ4_gn2_erg_m[,13] + par_postJ4_gn2_erg_m[,16] + par_postJ4_gn2_erg_m[,19]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ4_mn1_erg_m[,7] + par_postJ4_mn1_erg_m[,10] + par_postJ4_mn1_erg_m[,13] + par_postJ4_mn1_erg_m[,16] + par_postJ4_mn1_erg_m[,19]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ4_mn2_erg_m[,7] + par_postJ4_mn2_erg_m[,10] + par_postJ4_mn2_erg_m[,13] + par_postJ4_mn2_erg_m[,16] + par_postJ4_mn2_erg_m[,19]), type="l", col="green")   #, lty=4

plot(zeta036912_J4_gn1, xlab = "z0+z3+z6+z9+z12", main=" ", ylim=c(0, 45)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta036912_J4_gn2, col="red")
lines(zeta036912_J4_mn1, col="blue")
lines(zeta036912_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,7] + par_postJ4_gn1_th[,10] + par_postJ4_gn1_th[,13] + par_postJ4_gn1_th[,16] + par_postJ4_gn1_th[,19]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,7] + par_postJ4_gn2_th[,10] + par_postJ4_gn2_th[,13] + par_postJ4_gn2_th[,16] + par_postJ4_gn2_th[,19]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,7] + par_postJ4_mn1_th[,10] + par_postJ4_mn1_th[,13] + par_postJ4_mn1_th[,16] + par_postJ4_mn1_th[,19]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,7] + par_postJ4_mn2_th[,10] + par_postJ4_mn2_th[,13] + par_postJ4_mn2_th[,16] + par_postJ4_mn2_th[,19]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C0
plot(1:100000, (par_postJ4_gibbs_norm1[,22] + par_postJ4_gibbs_norm1[,25] + par_postJ4_gibbs_norm1[,28] + par_postJ4_gibbs_norm1[,31] + par_postJ4_gibbs_norm1[,34]), type = "l", xlab="Iteration", ylab = "z15+z18+z21+z24+z27", ylim=c(0,0.2))
lines(1:100000,(par_postJ4_gibbs_norm2[,22] + par_postJ4_gibbs_norm2[,25] + par_postJ4_gibbs_norm2[,28] + par_postJ4_gibbs_norm2[,31] + par_postJ4_gibbs_norm2[,34]), type = "l", col="red")
lines(1:100000, (par_postJ4_metr_norm1[,22] + par_postJ4_metr_norm1[,25] + par_postJ4_metr_norm1[,28] + par_postJ4_metr_norm1[,31] + par_postJ4_metr_norm1[,34]), type = "l", col="blue")
lines(1:100000, (par_postJ4_metr_norm2[,22] + par_postJ4_metr_norm2[,25] + par_postJ4_metr_norm2[,28] + par_postJ4_metr_norm2[,31] + par_postJ4_metr_norm2[,34]), type = "l", col="green")

plot(c(1:100000), (par_postJ4_gn1_erg_m[,22] + par_postJ4_gn1_erg_m[,25] + par_postJ4_gn1_erg_m[,28] + par_postJ4_gn1_erg_m[,31] + par_postJ4_gn1_erg_m[,34]), type="l", ylab = "z15+z18+z21+z24+z27", xlab="Iteration", ylim=c(0,0.2)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ4_gn2_erg_m[,22] + par_postJ4_gn2_erg_m[,25] + par_postJ4_gn2_erg_m[,28] + par_postJ4_gn2_erg_m[,31] + par_postJ4_gn2_erg_m[,34]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ4_mn1_erg_m[,22] + par_postJ4_mn1_erg_m[,25] + par_postJ4_mn1_erg_m[,28] + par_postJ4_mn1_erg_m[,31] + par_postJ4_mn1_erg_m[,34]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ4_mn2_erg_m[,22] + par_postJ4_mn2_erg_m[,25] + par_postJ4_mn2_erg_m[,28] + par_postJ4_mn2_erg_m[,31] + par_postJ4_mn2_erg_m[,34]), type="l", col="green")   #, lty=4

plot(zeta1518212427_J4_gn1, xlab = "z15+z18+z21+z24+z27", main=" ", ylim=c(0, 45)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1518212427_J4_gn2, col="red")
lines(zeta1518212427_J4_mn1, col="blue")
lines(zeta1518212427_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,22] + par_postJ4_gn1_th[,25] + par_postJ4_gn1_th[,28] + par_postJ4_gn1_th[,31] + par_postJ4_gn1_th[,34]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,22] + par_postJ4_gn2_th[,25] + par_postJ4_gn2_th[,28] + par_postJ4_gn2_th[,31] + par_postJ4_gn2_th[,34]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,22] + par_postJ4_mn1_th[,25] + par_postJ4_mn1_th[,28] + par_postJ4_mn1_th[,31] + par_postJ4_mn1_th[,34]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,22] + par_postJ4_mn2_th[,25] + par_postJ4_mn2_th[,28] + par_postJ4_mn2_th[,31] + par_postJ4_mn2_th[,34]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BL C1
plot(1:100000, (par_postJ4_gibbs_norm1[,8] + par_postJ4_gibbs_norm1[,11] + par_postJ4_gibbs_norm1[,14] + par_postJ4_gibbs_norm1[,17] + par_postJ4_gibbs_norm1[,20]), type = "l", xlab="Iteration", ylab = "z1+z4+z7+z10+z13", ylim=c(0.25,0.75))
lines(1:100000,(par_postJ4_gibbs_norm2[,8] + par_postJ4_gibbs_norm2[,11] + par_postJ4_gibbs_norm2[,14] + par_postJ4_gibbs_norm2[,17] + par_postJ4_gibbs_norm2[,20]), type = "l", col="red")
lines(1:100000, (par_postJ4_metr_norm1[,8] + par_postJ4_metr_norm1[,11] + par_postJ4_metr_norm1[,14] + par_postJ4_metr_norm1[,17] + par_postJ4_metr_norm1[,20]), type = "l", col="blue")
lines(1:100000, (par_postJ4_metr_norm2[,8] + par_postJ4_metr_norm2[,11] + par_postJ4_metr_norm2[,14] + par_postJ4_metr_norm2[,17] + par_postJ4_metr_norm2[,20]), type = "l", col="green")

plot(c(1:100000), (par_postJ4_gn1_erg_m[,8] + par_postJ4_gn1_erg_m[,11] + par_postJ4_gn1_erg_m[,14] + par_postJ4_gn1_erg_m[,17] + par_postJ4_gn1_erg_m[,20]), type="l", ylab = "z1+z4+z7+z10+z13", xlab="Iteration", ylim=c(0.25,0.75)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ4_gn2_erg_m[,8] + par_postJ4_gn2_erg_m[,11] + par_postJ4_gn2_erg_m[,14] + par_postJ4_gn2_erg_m[,17] + par_postJ4_gn2_erg_m[,20]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ4_mn1_erg_m[,8] + par_postJ4_mn1_erg_m[,11] + par_postJ4_mn1_erg_m[,14] + par_postJ4_mn1_erg_m[,17] + par_postJ4_mn1_erg_m[,20]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ4_mn2_erg_m[,8] + par_postJ4_mn2_erg_m[,11] + par_postJ4_mn2_erg_m[,14] + par_postJ4_mn2_erg_m[,17] + par_postJ4_mn2_erg_m[,20]), type="l", col="green")   #, lty=4

plot(zeta1471013_J4_gn1, xlab = "z1+z4+z7+z10+z13", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1471013_J4_gn2, col="red")
lines(zeta1471013_J4_mn1, col="blue")
lines(zeta1471013_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,8] + par_postJ4_gn1_th[,11] + par_postJ4_gn1_th[,14] + par_postJ4_gn1_th[,17] + par_postJ4_gn1_th[,20]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,8] + par_postJ4_gn2_th[,11] + par_postJ4_gn2_th[,14] + par_postJ4_gn2_th[,17] + par_postJ4_gn2_th[,20]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,8] + par_postJ4_mn1_th[,11] + par_postJ4_mn1_th[,14] + par_postJ4_mn1_th[,17] + par_postJ4_mn1_th[,20]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,8] + par_postJ4_mn2_th[,11] + par_postJ4_mn2_th[,14] + par_postJ4_mn2_th[,17] + par_postJ4_mn2_th[,20]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C1
plot(1:100000, (par_postJ4_gibbs_norm1[,23] + par_postJ4_gibbs_norm1[,26] + par_postJ4_gibbs_norm1[,29] + par_postJ4_gibbs_norm1[,32] + par_postJ4_gibbs_norm1[,35]), type = "l", xlab="Iteration", ylab = "z16+z19+z22+z25+z28", ylim=c(0,0.4))
lines(1:100000,(par_postJ4_gibbs_norm2[,23] + par_postJ4_gibbs_norm2[,26] + par_postJ4_gibbs_norm2[,29] + par_postJ4_gibbs_norm2[,32] + par_postJ4_gibbs_norm2[,35]), type = "l", col="red")
lines(1:100000, (par_postJ4_metr_norm1[,23] + par_postJ4_metr_norm1[,26] + par_postJ4_metr_norm1[,29] + par_postJ4_metr_norm1[,32] + par_postJ4_metr_norm1[,35]), type = "l", col="blue")
lines(1:100000, (par_postJ4_metr_norm2[,23] + par_postJ4_metr_norm2[,26] + par_postJ4_metr_norm2[,29] + par_postJ4_metr_norm2[,32] + par_postJ4_metr_norm2[,35]), type = "l", col="green")

plot(c(1:100000), (par_postJ4_gn1_erg_m[,23] + par_postJ4_gn1_erg_m[,26] + par_postJ4_gn1_erg_m[,29] + par_postJ4_gn1_erg_m[,32] + par_postJ4_gn1_erg_m[,35]), type="l", ylab = "z16+z19+z22+z25+z28", xlab="Iteration", ylim=c(0,0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ4_gn2_erg_m[,23] + par_postJ4_gn2_erg_m[,26] + par_postJ4_gn2_erg_m[,29] + par_postJ4_gn2_erg_m[,32] + par_postJ4_gn2_erg_m[,35]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ4_mn1_erg_m[,23] + par_postJ4_mn1_erg_m[,26] + par_postJ4_mn1_erg_m[,29] + par_postJ4_mn1_erg_m[,32] + par_postJ4_mn1_erg_m[,35]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ4_mn2_erg_m[,23] + par_postJ4_mn2_erg_m[,26] + par_postJ4_mn2_erg_m[,29] + par_postJ4_mn2_erg_m[,32] + par_postJ4_mn2_erg_m[,35]), type="l", col="green")   #, lty=4

plot(zeta1619222528_J4_gn1, xlab = "z16+z19+z22+z25+z28", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1619222528_J4_gn2, col="red")
lines(zeta1619222528_J4_mn1, col="blue")
lines(zeta1619222528_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,23] + par_postJ4_gn1_th[,26] + par_postJ4_gn1_th[,29] + par_postJ4_gn1_th[,32] + par_postJ4_gn1_th[,35]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,23] + par_postJ4_gn2_th[,26] + par_postJ4_gn2_th[,29] + par_postJ4_gn2_th[,32] + par_postJ4_gn2_th[,35]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,23] + par_postJ4_mn1_th[,26] + par_postJ4_mn1_th[,29] + par_postJ4_mn1_th[,32] + par_postJ4_mn1_th[,35]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,23] + par_postJ4_mn2_th[,26] + par_postJ4_mn2_th[,29] + par_postJ4_mn2_th[,32] + par_postJ4_mn2_th[,35]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BL C2
plot(1:100000, (par_postJ4_gibbs_norm1[,9] + par_postJ4_gibbs_norm1[,12] + par_postJ4_gibbs_norm1[,15] + par_postJ4_gibbs_norm1[,18] + par_postJ4_gibbs_norm1[,21]), type = "l", xlab="Iteration", ylab = "z2+z5+z8+z11+z14", ylim=c(0,0.4))
lines(1:100000,(par_postJ4_gibbs_norm2[,9] + par_postJ4_gibbs_norm2[,12] + par_postJ4_gibbs_norm2[,15] + par_postJ4_gibbs_norm2[,18] + par_postJ4_gibbs_norm2[,21]), type = "l", col="red")
lines(1:100000, (par_postJ4_metr_norm1[,9] + par_postJ4_metr_norm1[,12] + par_postJ4_metr_norm1[,15] + par_postJ4_metr_norm1[,18] + par_postJ4_metr_norm1[,21]), type = "l", col="blue")
lines(1:100000, (par_postJ4_metr_norm2[,9] + par_postJ4_metr_norm2[,12] + par_postJ4_metr_norm2[,15] + par_postJ4_metr_norm2[,18] + par_postJ4_metr_norm2[,21]), type = "l", col="green")

plot(c(1:100000), (par_postJ4_gn1_erg_m[,9] + par_postJ4_gn1_erg_m[,12] + par_postJ4_gn1_erg_m[,15] + par_postJ4_gn1_erg_m[,18] + par_postJ4_gn1_erg_m[,21]), type="l", ylab = "z2+z5+z8+z11+z14", xlab="Iteration", ylim=c(0,0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ4_gn2_erg_m[,9] + par_postJ4_gn2_erg_m[,12] + par_postJ4_gn2_erg_m[,15] + par_postJ4_gn2_erg_m[,18] + par_postJ4_gn2_erg_m[,21]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ4_mn1_erg_m[,9] + par_postJ4_mn1_erg_m[,12] + par_postJ4_mn1_erg_m[,15] + par_postJ4_mn1_erg_m[,18] + par_postJ4_mn1_erg_m[,21]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ4_mn2_erg_m[,9] + par_postJ4_mn2_erg_m[,12] + par_postJ4_mn2_erg_m[,15] + par_postJ4_mn2_erg_m[,18] + par_postJ4_mn2_erg_m[,21]), type="l", col="green")   #, lty=4

plot(zeta2581114_J4_gn1, xlab = "z2+z5+z8+z11+z14", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta2581114_J4_gn2, col="red")
lines(zeta2581114_J4_mn1, col="blue")
lines(zeta2581114_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,9] + par_postJ4_gn1_th[,12] + par_postJ4_gn1_th[,15] + par_postJ4_gn1_th[,18] + par_postJ4_gn1_th[,21]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,9] + par_postJ4_gn2_th[,12] + par_postJ4_gn2_th[,15] + par_postJ4_gn2_th[,18] + par_postJ4_gn2_th[,21]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,9] + par_postJ4_mn1_th[,12] + par_postJ4_mn1_th[,15] + par_postJ4_mn1_th[,18] + par_postJ4_mn1_th[,21]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,9] + par_postJ4_mn2_th[,12] + par_postJ4_mn2_th[,15] + par_postJ4_mn2_th[,18] + par_postJ4_mn2_th[,21]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)

## - BH C2
plot(1:100000, (par_postJ4_gibbs_norm1[,24] + par_postJ4_gibbs_norm1[,27] + par_postJ4_gibbs_norm1[,30] + par_postJ4_gibbs_norm1[,33] + par_postJ4_gibbs_norm1[,36]), type = "l", xlab="Iteration", ylab = "z17+z20+z23+z26+z29", ylim=c(0,0.4))
lines(1:100000,(par_postJ4_gibbs_norm2[,24] + par_postJ4_gibbs_norm2[,27] + par_postJ4_gibbs_norm2[,30] + par_postJ4_gibbs_norm2[,33] + par_postJ4_gibbs_norm2[,36]), type = "l", col="red")
lines(1:100000, (par_postJ4_metr_norm1[,24] + par_postJ4_metr_norm1[,27] + par_postJ4_metr_norm1[,30] + par_postJ4_metr_norm1[,33] + par_postJ4_metr_norm1[,36]), type = "l", col="blue")
lines(1:100000, (par_postJ4_metr_norm2[,24] + par_postJ4_metr_norm2[,27] + par_postJ4_metr_norm2[,30] + par_postJ4_metr_norm2[,33] + par_postJ4_metr_norm2[,36]), type = "l", col="green")

plot(c(1:100000), (par_postJ4_gn1_erg_m[,24] + par_postJ4_gn1_erg_m[,27] + par_postJ4_gn1_erg_m[,30] + par_postJ4_gn1_erg_m[,33] + par_postJ4_gn1_erg_m[,36]), type="l", ylab = "z17+z20+z23+z26+z29", xlab="Iteration", ylim=c(0, 0.4)) #, ylim=c(0,1)
lines(c(1:100000),(par_postJ4_gn2_erg_m[,24] + par_postJ4_gn2_erg_m[,27] + par_postJ4_gn2_erg_m[,30] + par_postJ4_gn2_erg_m[,33] + par_postJ4_gn2_erg_m[,36]), type="l", col="red") #, lty=2
lines(c(1:100000),(par_postJ4_mn1_erg_m[,24] + par_postJ4_mn1_erg_m[,27] + par_postJ4_mn1_erg_m[,30] + par_postJ4_mn1_erg_m[,33] + par_postJ4_mn1_erg_m[,36]), type="l", col="blue")   #, lty=3
lines(c(1:100000),(par_postJ4_mn2_erg_m[,24] + par_postJ4_mn2_erg_m[,27] + par_postJ4_mn2_erg_m[,30] + par_postJ4_mn2_erg_m[,33] + par_postJ4_mn2_erg_m[,36]), type="l", col="green")   #, lty=4

plot(zeta1720232629_J4_gn1, xlab = "z17+z20+z23+z26+z29", main=" ", ylim=c(0, 20)) #, main="Density of the posterior distribution of alpha", ylim=c(0, 7)
lines(zeta1720232629_J4_gn2, col="red")
lines(zeta1720232629_J4_mn1, col="blue")
lines(zeta1720232629_J4_mn2, col="green")
abline(v=mean(par_postJ4_gn1_th[,24] + par_postJ4_gn1_th[,27] + par_postJ4_gn1_th[,30] + par_postJ4_gn1_th[,33] + par_postJ4_gn1_th[,36]), lty=2) # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_gn2_th[,24] + par_postJ4_gn2_th[,27] + par_postJ4_gn2_th[,30] + par_postJ4_gn2_th[,33] + par_postJ4_gn2_th[,36]), lty=2, col="red") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn1_th[,24] + par_postJ4_mn1_th[,27] + par_postJ4_mn1_th[,30] + par_postJ4_mn1_th[,33] + par_postJ4_mn1_th[,36]), lty=2, col="blue") # abline(v=log(phi0_p), lty=2)
abline(v=mean(par_postJ4_mn2_th[,24] + par_postJ4_mn2_th[,27] + par_postJ4_mn2_th[,30] + par_postJ4_mn2_th[,33] + par_postJ4_mn2_th[,36]), lty=2, col="green") # abline(v=log(phi0_p), lty=2)


################################################################ - WAIC and WBIC - ###########################################################

# - J=0
S_sample_J0_gn1 <- par_postJ0_gn1_th#[sample(nrow(par_postJ0_gn1_th), 1000),]
S_sample_J0_gn2 <- par_postJ0_gn2_th#[sample(nrow(par_postJ0_gn2_th), 1000),]
S_sample_J0_mn1 <- par_postJ0_mn1_th#[sample(nrow(par_postJ0_mn1_th), 1000),]
S_sample_J0_mn2 <- par_postJ0_mn2_th#[sample(nrow(par_postJ0_mn2_th), 1000),]

log_avg_J0_gn1 <- rep(0, nrow(DS_C))
avg_log_J0_gn1 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_gn1[,1]
b <- S_sample_J0_gn1[,2]

for (i in 1:nrow(DS_C)){ 
  log_avg_J0_gn1[i] <- log(mean((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))  
}

s_log_avg_J0_gn1 <- sum(log_avg_J0_gn1)

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_gn1[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

s_avg_log_J0_gn1 <- sum(avg_log_J0_gn1)

p0_W_gn1 <- 2 * (s_log_avg_J0_gn1 - s_avg_log_J0_gn1)

WAIC_J0_gn1 <- - s_log_avg_J0_gn1 + p0_W_gn1

#######

log_avg_J0_gn2 <- rep(0, nrow(DS_C))
avg_log_J0_gn2 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_gn2[,1]
b <- S_sample_J0_gn2[,2]

for (i in 1:nrow(DS_C)){ 
  log_avg_J0_gn2[i] <- log(mean((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))  
}

s_log_avg_J0_gn2 <- sum(log_avg_J0_gn2)

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_gn2[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

s_avg_log_J0_gn2 <- sum(avg_log_J0_gn2)

p0_W_gn2 <- 2 * (s_log_avg_J0_gn2 - s_avg_log_J0_gn2)

WAIC_J0_gn2 <- - s_log_avg_J0_gn2 + p0_W_gn2


###### - MH

log_avg_J0_mn1 <- rep(0, nrow(DS_C))
avg_log_J0_mn1 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_mn1[,1]
b <- S_sample_J0_mn1[,2]


for (i in 1:nrow(DS_C)){ 
  log_avg_J0_mn1[i] <- log(mean((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))  
}

s_log_avg_J0_mn1 <- sum(log_avg_J0_mn1)

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_mn1[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

s_avg_log_J0_mn1 <- sum(avg_log_J0_mn1)

p0_W_mn1 <- 2 * (s_log_avg_J0_mn1 - s_avg_log_J0_mn1)

WAIC_J0_mn1 <- - s_log_avg_J0_mn1 + p0_W_mn1


#

log_avg_J0_mn2 <- rep(0, nrow(DS_C))
avg_log_J0_mn2 <- rep(0, nrow(DS_C))

phi0 <- S_sample_J0_mn2[,1]
b <- S_sample_J0_mn2[,2]


for (i in 1:nrow(DS_C)){ 
  log_avg_J0_mn2[i] <- log(mean((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))  
}

s_log_avg_J0_mn2 <- sum(log_avg_J0_mn2)

for (i in 1:nrow(DS_C)){ 
  avg_log_J0_mn2[i] <- mean(log((exp( -((exp( b * t[i]) - 1) /  b ) * phi0 * exp( b * (x[i] - 77.5) ) ) * ((phi0 * exp( b * (x[i] - 77.5 + t[i])  ))^d[i]) )))
}

s_avg_log_J0_mn2 <- sum(avg_log_J0_mn2)

p0_W_mn2 <- 2 * (s_log_avg_J0_mn2 - s_avg_log_J0_mn2)

WAIC_J0_mn2 <- - s_log_avg_J0_mn2 + p0_W_mn2


# - J=1
S_sample_J1_gn1 <- par_postJ1_gn1_th#[sample(nrow(par_postJ1_gn1_th), 1000),]
S_sample_J1_gn2 <- par_postJ1_gn2_th#[sample(nrow(par_postJ1_gn2_th), 1000),]
S_sample_J1_mn1 <- par_postJ1_mn1_th#[sample(nrow(par_postJ1_mn1_th), 1000),]
S_sample_J1_mn2 <- par_postJ1_mn2_th#[sample(nrow(par_postJ1_mn2_th), 1000),]


log_avg_J1_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J1_P1_gn1 <- rep(0, nrow(DS1))

log_avg_J1_P2_gn1 <- rep(0, nrow(DS2))
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
  log_avg_J1_P1_gn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7  + z8 ) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J1_P2_gn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_log_avg_J1_gn1 <- sum(log_avg_J1_P1_gn1) + sum(log_avg_J1_P2_gn1)
s_avg_log_J1_gn1 <- sum(avg_log_J1_P1_gn1) + sum(avg_log_J1_P2_gn1)

p1_W_gn1 <- 2 * (s_log_avg_J1_gn1 - s_avg_log_J1_gn1)

WAIC_J1_gn1 <- - s_log_avg_J1_gn1 + p1_W_gn1


############################################

log_avg_J1_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J1_P1_gn2 <- rep(0, nrow(DS1))

log_avg_J1_P2_gn2 <- rep(0, nrow(DS2))
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
  log_avg_J1_P1_gn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J1_P2_gn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_log_avg_J1_gn2 <- sum(log_avg_J1_P1_gn2) + sum(log_avg_J1_P2_gn2)

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_avg_log_J1_gn2 <- sum(avg_log_J1_P1_gn2) + sum(avg_log_J1_P2_gn2)

p1_W_gn2 <- 2 * (s_log_avg_J1_gn2 - s_avg_log_J1_gn2)

WAIC_J1_gn2 <- - s_log_avg_J1_gn2 + p1_W_gn2

###### - MH

log_avg_J1_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J1_P1_mn1 <- rep(0, nrow(DS1))

log_avg_J1_P2_mn1 <- rep(0, nrow(DS2))
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
  log_avg_J1_P1_mn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J1_P2_mn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_log_avg_J1_mn1 <- sum(log_avg_J1_P1_mn1) + sum(log_avg_J1_P2_mn1)

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_avg_log_J1_mn1 <- sum(avg_log_J1_P1_mn1) + sum(avg_log_J1_P2_mn1)

p1_W_mn1 <- 2 * (s_log_avg_J1_mn1 - s_avg_log_J1_mn1)

WAIC_J1_mn1 <- - s_log_avg_J1_mn1 + p1_W_mn1


################

log_avg_J1_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J1_P1_mn2 <- rep(0, nrow(DS1))

log_avg_J1_P2_mn2 <- rep(0, nrow(DS2))
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
  log_avg_J1_P1_mn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J1_P2_mn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_log_avg_J1_mn2 <- sum(log_avg_J1_P1_mn2) + sum(log_avg_J1_P2_mn2)

for (i in 1:nrow(DS1)){ 
  avg_log_J1_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0        * exp(b * (x1[i] - 77.5)) ) * ((phi0        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z6 + z7 + z8) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z9 + z10 + z11) * bH[i])) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5) + bH[i] * (z6 + z7 + z8 + z9 + z10 + z11))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J1_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0        * exp(b * (x2[i] - 77.5)) ) * ((phi0        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0 + z6) * c0[i] + (z1 + z7 ) * c1[i] + (z2 + z8 ) * c2[i]) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3 + z9) * c0[i] + (z4 + z10) * c1[i] + (z5 + z11) * c2[i])) / (c0[i] * (z0 + z3 + z6 + z9) + c1[i] * (z1 + z4 + z7 + z10) + c2[i] * (z2 + z5 + z8 + z11))  )  )
}

s_avg_log_J1_mn2 <- sum(avg_log_J1_P1_mn2) + sum(avg_log_J1_P2_mn2)

p1_W_mn2 <- 2 * (s_log_avg_J1_mn2 - s_avg_log_J1_mn2)

WAIC_J1_mn2 <- - s_log_avg_J1_mn2 + p1_W_mn2



# - J=2

S_sample_J2_gn1 <- par_postJ2_gn1_th#[sample(nrow(par_postJ2_gn1_th), 1000),]
S_sample_J2_gn2 <- par_postJ2_gn2_th#[sample(nrow(par_postJ2_gn2_th), 1000),]
S_sample_J2_mn1 <- par_postJ2_mn1_th#[sample(nrow(par_postJ2_mn1_th), 1000),]
S_sample_J2_mn2 <- par_postJ2_mn2_th#[sample(nrow(par_postJ2_mn2_th), 1000),]

log_avg_J2_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J2_P1_gn1 <- rep(0, nrow(DS1))

log_avg_J2_P2_gn1 <- rep(0, nrow(DS2))
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
  log_avg_J2_P1_gn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J2_P2_gn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}


s_log_avg_J2_gn1 <- sum(log_avg_J2_P1_gn1) + sum(log_avg_J2_P2_gn1)
s_avg_log_J2_gn1 <- sum(avg_log_J2_P1_gn1) + sum(avg_log_J2_P2_gn1)

p2_W_gn1 <- 2 * (s_log_avg_J2_gn1 - s_avg_log_J2_gn1)

WAIC_J2_gn1 <- - s_log_avg_J2_gn1 + p2_W_gn1


##########

log_avg_J2_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J2_P1_gn2 <- rep(0, nrow(DS1))

log_avg_J2_P2_gn2 <- rep(0, nrow(DS2))
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
  log_avg_J2_P1_gn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J2_P2_gn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}


s_log_avg_J2_gn2 <- sum(log_avg_J2_P1_gn2) + sum(log_avg_J2_P2_gn2)
s_avg_log_J2_gn2 <- sum(avg_log_J2_P1_gn2) + sum(avg_log_J2_P2_gn2)

p2_W_gn2 <- 2 * (s_log_avg_J2_gn2 - s_avg_log_J2_gn2)

WAIC_J2_gn2 <- - s_log_avg_J2_gn2 + p2_W_gn2


###########

log_avg_J2_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J2_P1_mn1 <- rep(0, nrow(DS1))

log_avg_J2_P2_mn1 <- rep(0, nrow(DS2))
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
  log_avg_J2_P1_mn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J2_P2_mn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}


s_log_avg_J2_mn1 <- sum(log_avg_J2_P1_mn1) + sum(log_avg_J2_P2_mn1)
s_avg_log_J2_mn1 <- sum(avg_log_J2_P1_mn1) + sum(avg_log_J2_P2_mn1)

p2_W_mn1 <- 2 * (s_log_avg_J2_mn1 - s_avg_log_J2_mn1)

WAIC_J2_mn1 <- - s_log_avg_J2_mn1 + p2_W_mn1


###########

log_avg_J2_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J2_P1_mn2 <- rep(0, nrow(DS1))

log_avg_J2_P2_mn2 <- rep(0, nrow(DS2))
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
  log_avg_J2_P1_mn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J2_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0               * exp(b * (x1[i] - 77.5)) ) * ((phi0               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0 + z1 + z2) * bL[i] + (z9  + z10 + z11) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3 + z4 + z5) * bL[i] + (z12 + z13 + z14) * bH[i]) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6 + z7 + z8) * bL[i] + (z15 + z16 + z17) * bH[i]) ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) + bH[i] * (z9 + z10 + z11 + z12 + z13 + z14 + z15 + z16 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J2_P2_mn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J2_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0               * exp(b * (x2[i] - 77.5)) ) * ((phi0               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z9 ) + c1[i] * (z1 + z10) + c2[i] * (z2 + z11)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z12) + c1[i] * (z4 + z13) + c2[i] * (z5 + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z15) + c1[i] * (z7 + z16) + c2[i] * (z8 + z17)) ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17))  )  )
}


s_log_avg_J2_mn2 <- sum(log_avg_J2_P1_mn2) + sum(log_avg_J2_P2_mn2)
s_avg_log_J2_mn2 <- sum(avg_log_J2_P1_mn2) + sum(avg_log_J2_P2_mn2)

p2_W_mn2 <- 2 * (s_log_avg_J2_mn2 - s_avg_log_J2_mn2)

WAIC_J2_mn2 <- - s_log_avg_J2_mn2 + p2_W_mn2


# - J=3

S_sample_J3_gn1 <- par_postJ3_gn1_th#[sample(nrow(par_postJ3_gn1_th), 1000),]
S_sample_J3_gn2 <- par_postJ3_gn2_th#[sample(nrow(par_postJ3_gn2_th), 1000),]
S_sample_J3_mn1 <- par_postJ3_mn1_th#[sample(nrow(par_postJ3_mn1_th), 1000),]
S_sample_J3_mn2 <- par_postJ3_mn2_th#[sample(nrow(par_postJ3_mn2_th), 1000),]

log_avg_J3_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J3_P1_gn1 <- rep(0, nrow(DS1))

log_avg_J3_P2_gn1 <- rep(0, nrow(DS2))
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
  log_avg_J3_P1_gn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J3_P2_gn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

s_log_avg_J3_gn1 <- sum(log_avg_J3_P1_gn1) + sum(log_avg_J3_P2_gn1)
s_avg_log_J3_gn1 <- sum(avg_log_J3_P1_gn1) + sum(avg_log_J3_P2_gn1)

p3_W_gn1 <- 2 * (s_log_avg_J3_gn1 - s_avg_log_J3_gn1)

WAIC_J3_gn1 <- - s_log_avg_J3_gn1 + p3_W_gn1


####################################################################

log_avg_J3_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J3_P1_gn2 <- rep(0, nrow(DS1))

log_avg_J3_P2_gn2 <- rep(0, nrow(DS2))
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
  log_avg_J3_P1_gn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J3_P2_gn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

s_log_avg_J3_gn2 <- sum(log_avg_J3_P1_gn2) + sum(log_avg_J3_P2_gn2)
s_avg_log_J3_gn2 <- sum(avg_log_J3_P1_gn2) + sum(avg_log_J3_P2_gn2)

p3_W_gn2 <- 2 * (s_log_avg_J3_gn2 - s_avg_log_J3_gn2)

WAIC_J3_gn2 <- - s_log_avg_J3_gn2 + p3_W_gn2


####################################################################

log_avg_J3_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J3_P1_mn1 <- rep(0, nrow(DS1))

log_avg_J3_P2_mn1 <- rep(0, nrow(DS2))
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
  log_avg_J3_P1_mn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J3_P2_mn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

s_log_avg_J3_mn1 <- sum(log_avg_J3_P1_mn1) + sum(log_avg_J3_P2_mn1)
s_avg_log_J3_mn1 <- sum(avg_log_J3_P1_mn1) + sum(avg_log_J3_P2_mn1)

p3_W_mn1 <- 2 * (s_log_avg_J3_mn1 - s_avg_log_J3_mn1)

WAIC_J3_mn1 <- - s_log_avg_J3_mn1 + p3_W_mn1


####################################################################

log_avg_J3_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J3_P1_mn2 <- rep(0, nrow(DS1))

log_avg_J3_P2_mn2 <- rep(0, nrow(DS2))
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
  log_avg_J3_P1_mn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J3_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                      * exp(b * (x1[i] - 77.5)) ) * ((phi0                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0 + z1  + z2 ) + bH[i] * (z12 + z13 + z14)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3 + z4  + z5 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6 + z7  + z8 ) + bH[i] * (z18 + z19 + z20)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9 + z10 + z11) + bH[i] * (z21 + z22 + z23))  ) / (bL[i] * (z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + z9 + z10 + z11) + bH[i] * (z12 + z13 + z14 + z15 + z16 + z17 + z18 + z19 + z20 + z21 + z22 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J3_P2_mn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J3_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                      * exp(b * (x2[i] - 77.5)) ) * ((phi0                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0 + z12) + c1[i] * (z1  + z13) + c2[i] * (z2  + z14)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3 + z15) + c1[i] * (z4  + z16) + c2[i] * (z5  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6 + z18) + c1[i] * (z7  + z19) + c2[i] * (z8  + z20)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9 + z21) + c1[i] * (z10 + z22) + c2[i] * (z11 + z23))  ) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23))  )  )
}

s_log_avg_J3_mn2 <- sum(log_avg_J3_P1_mn2) + sum(log_avg_J3_P2_mn2)
s_avg_log_J3_mn2 <- sum(avg_log_J3_P1_mn2) + sum(avg_log_J3_P2_mn2)

p3_W_mn2 <- 2 * (s_log_avg_J3_mn2 - s_avg_log_J3_mn2)

WAIC_J3_mn2 <- - s_log_avg_J3_mn2 + p3_W_mn2


# - J=4

S_sample_J4_gn1 <- par_postJ4_gn1_th#[sample(nrow(par_postJ4_gn1_th), 1000),]
S_sample_J4_gn2 <- par_postJ4_gn2_th#[sample(nrow(par_postJ4_gn2_th), 1000),]
S_sample_J4_mn1 <- par_postJ4_mn1_th#[sample(nrow(par_postJ4_mn1_th), 1000),]
S_sample_J4_mn2 <- par_postJ4_mn2_th#[sample(nrow(par_postJ4_mn2_th), 1000),]

log_avg_J4_P1_gn1 <- rep(0, nrow(DS1))
avg_log_J4_P1_gn1 <- rep(0, nrow(DS1))

log_avg_J4_P2_gn1 <- rep(0, nrow(DS2))
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
  log_avg_J4_P1_gn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_gn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                    exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J4_P2_gn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_gn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                    exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

s_log_avg_J4_gn1 <- sum(log_avg_J4_P1_gn1) + sum(log_avg_J4_P2_gn1)
s_avg_log_J4_gn1 <- sum(avg_log_J4_P1_gn1) + sum(avg_log_J4_P2_gn1)

p4_W_gn1 <- 2 * (s_log_avg_J4_gn1 - s_avg_log_J4_gn1)

WAIC_J4_gn1 <- - s_log_avg_J4_gn1 + p4_W_gn1


####

log_avg_J4_P1_gn2 <- rep(0, nrow(DS1))
avg_log_J4_P1_gn2 <- rep(0, nrow(DS1))

log_avg_J4_P2_gn2 <- rep(0, nrow(DS2))
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
  log_avg_J4_P1_gn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_gn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J4_P2_gn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_gn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

s_log_avg_J4_gn2 <- sum(log_avg_J4_P1_gn2) + sum(log_avg_J4_P2_gn2)
s_avg_log_J4_gn2 <- sum(avg_log_J4_P1_gn2) + sum(avg_log_J4_P2_gn2)

p4_W_gn2 <- 2 * (s_log_avg_J4_gn2 - s_avg_log_J4_gn2)

WAIC_J4_gn2 <- - s_log_avg_J4_gn2 + p4_W_gn2



## - Metropolis

log_avg_J4_P1_mn1 <- rep(0, nrow(DS1))
avg_log_J4_P1_mn1 <- rep(0, nrow(DS1))

log_avg_J4_P2_mn1 <- rep(0, nrow(DS2))
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
  log_avg_J4_P1_mn1[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_mn1[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J4_P2_mn1[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_mn1[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

s_log_avg_J4_mn1 <- sum(log_avg_J4_P1_mn1) + sum(log_avg_J4_P2_mn1)
s_avg_log_J4_mn1 <- sum(avg_log_J4_P1_mn1) + sum(avg_log_J4_P2_mn1)

p4_W_mn1 <- 2 * (s_log_avg_J4_mn1 - s_avg_log_J4_mn1)

WAIC_J4_mn1 <- - s_log_avg_J4_mn1 + p4_W_mn1


####

log_avg_J4_P1_mn2 <- rep(0, nrow(DS1))
avg_log_J4_P1_mn2 <- rep(0, nrow(DS1))

log_avg_J4_P2_mn2 <- rep(0, nrow(DS2))
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
  log_avg_J4_P1_mn2[i] <- log(mean((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS1)){ 
  avg_log_J4_P1_mn2[i] <- mean(log((exp( -((exp(b * t1[i]) - 1) / b) * phi0                             * exp(b * (x1[i] - 77.5)) ) * ((phi0                             * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0  + z1  + z2 ) + bH[i] * (z15 + z16 + z17)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3  + z4  + z5 ) + bH[i] * (z18 + z19 + z20)) + 
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6  + z7  + z8 ) + bH[i] * (z21 + z22 + z23)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9  + z10 + z11) + bH[i] * (z24 + z25 + z26)) +
                                      exp( -((exp(b * t1[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12 + z13 + z14) + bH[i] * (z27 + z28 + z29))) / (bL[i] * (z0 + z1 +z2 + z3 + z4 + z5 + z6 +z7 + z8 + z9 + z10 + z11 +z12 + z13 + z14) + bH[i] * (z15 + z16 +z17 + z18 + z19 + z20 + z21 +z22 + z23 + z24 + z25 + z26 +z27 + z28 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  log_avg_J4_P2_mn2[i] <- log(mean((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

for (i in 1:nrow(DS2)){ 
  avg_log_J4_P2_mn2[i] <- mean(log((exp( -((exp(b * t2[i]) - 1) / b) * phi0                             * exp(b * (x2[i] - 77.5)) ) * ((phi0                             * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0  + z15) + c1[i] * (z1  + z16) + c2[i] * (z2  + z17)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1                      * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1                      * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3  + z18) + c1[i] * (z4  + z19) + c2[i] * (z5  + z20)) + 
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2               * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6  + z21) + c1[i] * (z7  + z22) + c2[i] * (z8  + z23)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3        * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9  + z24) + c1[i] * (z10 + z25) + c2[i] * (z11 + z26)) +
                                      exp( -((exp(b * t2[i]) - 1) / b) * phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5)) ) * ((phi0 * phi1 * phi2 * phi3 * phi4 * exp(b * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12 + z27) + c1[i] * (z13 + z28) + c2[i] * (z14 + z29))) / (c0[i] * (z0 + z3 + z6 + z9 + z12 + z15 + z18 + z21 + z24 + z27) + c1[i] * (z1 + z4 + z7 + z10 + z13 + z16 + z19 + z22 + z25 + z28) + c2[i] * (z2 + z5 + z8 + z11 + z14 + z17 + z20 + z23 + z26 + z29))  )  )
}

s_log_avg_J4_mn2 <- sum(log_avg_J4_P1_mn2) + sum(log_avg_J4_P2_mn2)
s_avg_log_J4_mn2 <- sum(avg_log_J4_P1_mn2) + sum(avg_log_J4_P2_mn2)

p4_W_mn2 <- 2 * (s_log_avg_J4_mn2 - s_avg_log_J4_mn2)

WAIC_J4_mn2 <- - s_log_avg_J4_mn2 + p4_W_mn2

# - Info Crit Table

Info_Crit_table <- matrix(NA, 12, 5)

## - J = 0
Info_Crit_table[1,1] <- WAIC_J0_gn1
Info_Crit_table[2,1] <- p0_W_gn1
#Info_Crit_table[3,1] <- WBIC_J0_gn1
Info_Crit_table[4,1] <- WAIC_J0_gn2
Info_Crit_table[5,1] <- p0_W_gn2
#Info_Crit_table[6,1] <- WBIC_J0_gn2
Info_Crit_table[7,1] <- WAIC_J0_mn1
Info_Crit_table[8,1] <- p0_W_mn1
#Info_Crit_table[9,1] <- WBIC_J0_mn1
Info_Crit_table[10,1] <- WAIC_J0_mn2
Info_Crit_table[11,1] <- p0_W_mn2
#Info_Crit_table[12,1] <- WBIC_J0_mn2

## - J = 1
Info_Crit_table[1,2] <- WAIC_J1_gn1
Info_Crit_table[2,2] <- p1_W_gn1
#Info_Crit_table[3,2] <- WBIC_J1_gn1
Info_Crit_table[4,2] <- WAIC_J1_gn2
Info_Crit_table[5,2] <- p1_W_gn2
#Info_Crit_table[6,2] <- WBIC_J1_gn2
Info_Crit_table[7,2] <- WAIC_J1_mn1
Info_Crit_table[8,2] <- p1_W_mn1
#Info_Crit_table[9,2] <- WBIC_J1_mn1
Info_Crit_table[10,2] <- WAIC_J1_mn2
Info_Crit_table[11,2] <- p1_W_mn2
#Info_Crit_table[12,2] <- WBIC_J1_mn2

## - J = 2
Info_Crit_table[1,3] <- WAIC_J2_gn1
Info_Crit_table[2,3] <- p2_W_gn1
#Info_Crit_table[3,3] <- WBIC_J2_gn1
Info_Crit_table[4,3] <- WAIC_J2_gn2
Info_Crit_table[5,3] <- p2_W_gn2
#Info_Crit_table[6,3] <- WBIC_J2_gn2
Info_Crit_table[7,3] <- WAIC_J2_mn1
Info_Crit_table[8,3] <- p2_W_mn1
#Info_Crit_table[9,3] <- WBIC_J2_mn1
Info_Crit_table[10,3] <- WAIC_J2_mn2
Info_Crit_table[11,3] <- p2_W_mn2
#Info_Crit_table[12,3] <- WBIC_J2_mn2

## - J = 3
Info_Crit_table[1,4] <- WAIC_J3_gn1
Info_Crit_table[2,4] <- p3_W_gn1
#Info_Crit_table[3,4] <- WBIC_J3_gn1
Info_Crit_table[4,4] <- WAIC_J3_gn2
Info_Crit_table[5,4] <- p3_W_gn2
#Info_Crit_table[6,4] <- WBIC_J3_gn2
Info_Crit_table[7,4] <- WAIC_J3_mn1
Info_Crit_table[8,4] <- p3_W_mn1
#Info_Crit_table[9,4] <- WBIC_J3_mn1
Info_Crit_table[10,4] <- WAIC_J3_mn2
Info_Crit_table[11,4] <- p3_W_mn2
#Info_Crit_table[12,4] <- WBIC_J3_mn2

## - J = 4
Info_Crit_table[1,5] <- WAIC_J4_gn1
Info_Crit_table[2,5] <- p4_W_gn1
#Info_Crit_table[3,5] <- WBIC_J4_gn1
Info_Crit_table[4,5] <- WAIC_J4_gn2
Info_Crit_table[5,5] <- p4_W_gn2
#Info_Crit_table[6,5] <- WBIC_J4_gn2
Info_Crit_table[7,5] <- WAIC_J4_mn1
Info_Crit_table[8,5] <- p4_W_mn1
#Info_Crit_table[9,5] <- WBIC_J4_mn1
Info_Crit_table[10,5] <- WAIC_J4_mn2
Info_Crit_table[11,5] <- p4_W_mn2
#Info_Crit_table[12,5] <- WBIC_J4_mn2

xtable(Info_Crit_table, digits=c(0,2,2,2,2,2))

############################# - Deviance comparison plot from Richardson and Green - ############################################

# J=0

phi0 <- par_postJ0_gn1_th[,1]
b <- par_postJ0_gn1_th[,2]

log_lik_J0 <- rep(0, 9000)

for(s in 1:9000){
  
  for (i in 1:nrow(DS_C)){ 
    log_lik_J0[s] <- log_lik_J0[s] + log((exp( -((exp( b[s] * t[i]) - 1) /  b[s] ) * phi0[s] * exp( b[s] * (x[i] - 77.5) ) ) * ((phi0[s] * exp( b[s] * (x[i] - 77.5 + t[i])  ))^d[i]) ))
  }
  
}

# J=1

phi0 <- par_postJ1_gn1_th[,1]
phi1 <- par_postJ1_gn1_th[,3]
b <- par_postJ1_gn1_th[,2]

z0 <- par_postJ1_gn1_th[,4]
z1 <- par_postJ1_gn1_th[,5]
z2 <- par_postJ1_gn1_th[,6]
z3 <- par_postJ1_gn1_th[,7]
z4 <- par_postJ1_gn1_th[,8]
z5 <- par_postJ1_gn1_th[,9]
z6 <- par_postJ1_gn1_th[,10]
z7 <- par_postJ1_gn1_th[,11]
z8 <- par_postJ1_gn1_th[,12]
z9 <- par_postJ1_gn1_th[,13]
z10 <- par_postJ1_gn1_th[,14]
z11 <- par_postJ1_gn1_th[,15]

log_lik_J1_P1 <- rep(0, 9000)
log_lik_J1_P2 <- rep(0, 9000)
log_lik_J1 <- rep(0, 9000)

for(s in 1:9000){
  
  for (i in 1:nrow(DS1)){ 
    log_lik_J1_P1[s] <- log_lik_J1_P1[s] + log((exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s]           * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s]           * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0[s] + z1[s] + z2[s]) * bL[i] + (z6[s] +  z7[s] +  z8[s]) * bH[i]) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3[s] + z4[s] + z5[s]) * bL[i] + (z9[s] + z10[s] + z11[s]) * bH[i])) / (bL[i] * (z0[s] + z1[s] + z2[s] + z3[s] + z4[s] + z5[s]) + bH[i] * (z6[s] + z7[s] + z8[s] + z9[s] + z10[s] + z11[s]))  )
  }
  
}

for(s in 1:9000){
  
  for (i in 1:nrow(DS2)){ 
    log_lik_J1_P2[s] <- log_lik_J1_P2[s] + log((exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s]           * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s]           * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z0[s] + z6[s]) * c0[i] + (z1[s] +  z7[s]) * c1[i] + (z2[s] +  z8[s]) * c2[i]) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * ((z3[s] + z9[s]) * c0[i] + (z4[s] + z10[s]) * c1[i] + (z5[s] + z11[s]) * c2[i])) / (c0[i] * (z0[s] + z3[s] + z6[s] + z9[s]) + c1[i] * (z1[s] + z4[s] + z7[s] + z10[s]) + c2[i] * (z2[s] + z5[s] + z8[s] + z11[s]))  )  
  }
  
}

log_lik_J1 <- log_lik_J1_P1 + log_lik_J1_P2


# J=2

phi0 <- par_postJ2_gn1_th[,1]
phi1 <- par_postJ2_gn1_th[,3]
phi2 <- par_postJ2_gn1_th[,4]
b <- par_postJ2_gn1_th[,2]

z0 <- par_postJ2_gn1_th[,5]
z1 <- par_postJ2_gn1_th[,6]
z2 <- par_postJ2_gn1_th[,7]
z3 <- par_postJ2_gn1_th[,8]
z4 <- par_postJ2_gn1_th[,9]
z5 <- par_postJ2_gn1_th[,10]
z6 <- par_postJ2_gn1_th[,11]
z7 <- par_postJ2_gn1_th[,12]
z8 <- par_postJ2_gn1_th[,13]
z9 <- par_postJ2_gn1_th[,14]
z10 <- par_postJ2_gn1_th[,15]
z11 <- par_postJ2_gn1_th[,16]
z12 <- par_postJ2_gn1_th[,17]
z13 <- par_postJ2_gn1_th[,18]
z14 <- par_postJ2_gn1_th[,19]
z15 <- par_postJ2_gn1_th[,20]
z16 <- par_postJ2_gn1_th[,21]
z17 <- par_postJ2_gn1_th[,22]

log_lik_J2_P1 <- rep(0, 9000)
log_lik_J2_P2 <- rep(0, 9000)
log_lik_J2 <- rep(0, 9000)

for(s in 1:9000){
  
  for (i in 1:nrow(DS1)){ 
    log_lik_J2_P1[s] <- log_lik_J2_P1[s] + log((exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s]                     * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s]                     * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z0[s] + z1[s] + z2[s]) * bL[i] + (z9[s]  + z10[s] + z11[s]) * bH[i]) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s]           * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s]           * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z3[s] + z4[s] + z5[s]) * bL[i] + (z12[s] + z13[s] + z14[s]) * bH[i]) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * ((z6[s] + z7[s] + z8[s]) * bL[i] + (z15[s] + z16[s] + z17[s]) * bH[i]) ) / (bL[i] * (z0[s] + z1[s] + z2[s] + z3[s] + z4[s] + z5[s] + z6[s] + z7[s] + z8[s]) + bH[i] * (z9[s] + z10[s] + z11[s] + z12[s] + z13[s] + z14[s] + z15[s] + z16[s] + z17[s]))  )
  }
  
}

for(s in 1:9000){
  
  for (i in 1:nrow(DS2)){ 
    log_lik_J2_P2[s] <- log_lik_J2_P2[s] + log((exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s]                     * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s]                     * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0[s] +  z9[s]) + c1[i] * (z1[s] + z10[s]) + c2[i] * (z2[s] + z11[s])) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s]           * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s]           * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3[s] + z12[s]) + c1[i] * (z4[s] + z13[s]) + c2[i] * (z5[s] + z14[s])) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6[s] + z15[s]) + c1[i] * (z7[s] + z16[s]) + c2[i] * (z8[s] + z17[s])) ) / (c0[i] * (z0[s] + z3[s] + z6[s] + z9[s] + z12[s] + z15[s]) + c1[i] * (z1[s] + z4[s] + z7[s] + z10[s] + z13[s] + z16[s]) + c2[i] * (z2[s] + z5[s] + z8[s] + z11[s] + z14[s] + z17[s]))  )
  }
  
}

log_lik_J2 <- log_lik_J2_P1 + log_lik_J2_P2


# J=3

phi0 <- par_postJ3_gn1_th[,1]
phi1 <- par_postJ3_gn1_th[,3]
phi2 <- par_postJ3_gn1_th[,4]
phi3 <- par_postJ3_gn1_th[,5]
b <- par_postJ3_gn1_th[,2]

z0 <- par_postJ3_gn1_th[,6]
z1 <- par_postJ3_gn1_th[,7]
z2 <- par_postJ3_gn1_th[,8]
z3 <- par_postJ3_gn1_th[,9]
z4 <- par_postJ3_gn1_th[,10]
z5 <- par_postJ3_gn1_th[,11]
z6 <- par_postJ3_gn1_th[,12]
z7 <- par_postJ3_gn1_th[,13]
z8 <- par_postJ3_gn1_th[,14]
z9 <- par_postJ3_gn1_th[,15]
z10 <- par_postJ3_gn1_th[,16]
z11 <- par_postJ3_gn1_th[,17]
z12 <- par_postJ3_gn1_th[,18]
z13 <- par_postJ3_gn1_th[,19]
z14 <- par_postJ3_gn1_th[,20]
z15 <- par_postJ3_gn1_th[,21]
z16 <- par_postJ3_gn1_th[,22]
z17 <- par_postJ3_gn1_th[,23]
z18 <- par_postJ3_gn1_th[,24]
z19 <- par_postJ3_gn1_th[,25]
z20 <- par_postJ3_gn1_th[,26]
z21 <- par_postJ3_gn1_th[,27]
z22 <- par_postJ3_gn1_th[,28]
z23 <- par_postJ3_gn1_th[,29]

log_lik_J3_P1 <- rep(0, 9000)
log_lik_J3_P2 <- rep(0, 9000)
log_lik_J3 <- rep(0, 9000)

for(s in 1:9000){
  
  for (i in 1:nrow(DS1)){ 
    log_lik_J3_P1[s] <- log_lik_J3_P1[s] + log((exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s]                               * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s]                               * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0[s] + z1[s]  + z2[s] ) + bH[i] * (z12[s] + z13[s] + z14[s])) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s]                     * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s]                     * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3[s] + z4[s]  + z5[s] ) + bH[i] * (z15[s] + z16[s] + z17[s])) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s]           * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s]           * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6[s] + z7[s]  + z8[s] ) + bH[i] * (z18[s] + z19[s] + z20[s])) +
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * phi3[s] * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * phi3[s] * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9[s] + z10[s] + z11[s]) + bH[i] * (z21[s] + z22[s] + z23[s]))  ) / (bL[i] * (z0[s] + z1[s] + z2[s] + z3[s] + z4[s] + z5[s] + z6[s] + z7[s] + z8[s] + z9[s] + z10[s] + z11[s]) + bH[i] * (z12[s] + z13[s] + z14[s] + z15[s] + z16[s] + z17[s] + z18[s] + z19[s] + z20[s] + z21[s] + z22[s] + z23[s]))  ) 
  }
  
}

for(s in 1:9000){
  
  for (i in 1:nrow(DS2)){ 
    log_lik_J3_P2[s] <- log_lik_J3_P2[s] + log((exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s]                               * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s]                               * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0[s] + z12[s]) + c1[i] * ( z1[s] + z13[s]) + c2[i] * ( z2[s] + z14[s])) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s]                     * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s]                     * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3[s] + z15[s]) + c1[i] * ( z4[s] + z16[s]) + c2[i] * ( z5[s] + z17[s])) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s]           * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s]           * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6[s] + z18[s]) + c1[i] * ( z7[s] + z19[s]) + c2[i] * ( z8[s] + z20[s])) +
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * phi3[s] * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * phi3[s] * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9[s] + z21[s]) + c1[i] * (z10[s] + z22[s]) + c2[i] * (z11[s] + z23[s]))  ) / (c0[i] * (z0[s] + z3[s] + z6[s] + z9[s] + z12[s] + z15[s] + z18[s] + z21[s]) + c1[i] * (z1[s] + z4[s] + z7[s] + z10[s] + z13[s] + z16[s] + z19[s] + z22[s]) + c2[i] * (z2[s] + z5[s] + z8[s] + z11[s] + z14[s] + z17[s] + z20[s] + z23[s]))  ) 
  }
  
}

log_lik_J3 <- log_lik_J3_P1 + log_lik_J3_P2


# J=4

phi0 <- par_postJ4_gn1_th[,1]
phi1 <- par_postJ4_gn1_th[,3]
phi2 <- par_postJ4_gn1_th[,4]
phi3 <- par_postJ4_gn1_th[,5]
phi4 <- par_postJ4_gn1_th[,6]
b <- par_postJ4_gn1_th[,2]

z0 <- par_postJ4_gn1_th[,7]
z1 <- par_postJ4_gn1_th[,8]
z2 <- par_postJ4_gn1_th[,9]
z3 <- par_postJ4_gn1_th[,10]
z4 <- par_postJ4_gn1_th[,11]
z5 <- par_postJ4_gn1_th[,12]
z6 <- par_postJ4_gn1_th[,13]
z7 <- par_postJ4_gn1_th[,14]
z8 <- par_postJ4_gn1_th[,15]
z9 <- par_postJ4_gn1_th[,16]

z10 <- par_postJ4_gn1_th[,17]
z11 <- par_postJ4_gn1_th[,18]
z12 <- par_postJ4_gn1_th[,19]
z13 <- par_postJ4_gn1_th[,20]
z14 <- par_postJ4_gn1_th[,21]
z15 <- par_postJ4_gn1_th[,22]
z16 <- par_postJ4_gn1_th[,23]
z17 <- par_postJ4_gn1_th[,24]
z18 <- par_postJ4_gn1_th[,25]
z19 <- par_postJ4_gn1_th[,26]
z20 <- par_postJ4_gn1_th[,27]
z21 <- par_postJ4_gn1_th[,28]
z22 <- par_postJ4_gn1_th[,29]
z23 <- par_postJ4_gn1_th[,30]
z24 <- par_postJ4_gn1_th[,31]
z25 <- par_postJ4_gn1_th[,32]
z26 <- par_postJ4_gn1_th[,33]
z27 <- par_postJ4_gn1_th[,34]
z28 <- par_postJ4_gn1_th[,35]
z29 <- par_postJ4_gn1_th[,36]

log_lik_J4_P1 <- rep(0, 9000)
log_lik_J4_P2 <- rep(0, 9000)
log_lik_J4 <- rep(0, 9000)

for(s in 1:9000){
  
  for (i in 1:nrow(DS1)){  
    log_lik_J4_P1[s] <- log_lik_J4_P1[s] + log((exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s]                                         * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s]                                         * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z0[s]  + z1[s]  + z2[s] ) + bH[i] * (z15[s] + z16[s] + z17[s])) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s]                               * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s]                               * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z3[s]  + z4[s]  + z5[s] ) + bH[i] * (z18[s] + z19[s] + z20[s])) + 
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s]                     * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s]                     * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z6[s]  + z7[s]  + z8[s] ) + bH[i] * (z21[s] + z22[s] + z23[s])) +
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * phi3[s]           * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * phi3[s]           * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z9[s]  + z10[s] + z11[s]) + bH[i] * (z24[s] + z25[s] + z26[s])) +
                                                  exp( -((exp(b[s] * t1[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * phi3[s] * phi4[s] * exp(b[s] * (x1[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * phi3[s] * phi4[s] * exp(b[s] * (x1[i] - 77.5 + t1[i])))^d1[i]) * (bL[i] * (z12[s] + z13[s] + z14[s]) + bH[i] * (z27[s] + z28[s] + z29[s]))) / (bL[i] * (z0[s] + z1[s] +z2[s] + z3[s] + z4[s] + z5[s] + z6[s] + z7[s] + z8[s] + z9[s] + z10[s] + z11[s] + z12[s] + z13[s] + z14[s]) + bH[i] * (z15[s] + z16[s] + z17[s] + z18[s] + z19[s] + z20[s] + z21[s] + z22[s] + z23[s] + z24[s] + z25[s] + z26[s] +z27[s] + z28[s] + z29[s]))  ) 
  }
  
}

for(s in 1:9000){
  
  for (i in 1:nrow(DS2)){ 
    log_lik_J4_P2[s] <- log_lik_J4_P2[s] + log((exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s]                                         * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s]                                         * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z0[s]  + z15[s]) + c1[i] * (z1[s]  + z16[s]) + c2[i] * (z2[s]  + z17[s])) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s]                               * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s]                               * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z3[s]  + z18[s]) + c1[i] * (z4[s]  + z19[s]) + c2[i] * (z5[s]  + z20[s])) + 
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s]                     * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s]                     * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z6[s]  + z21[s]) + c1[i] * (z7[s]  + z22[s]) + c2[i] * (z8[s]  + z23[s])) +
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * phi3[s]           * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * phi3[s]           * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z9[s]  + z24[s]) + c1[i] * (z10[s] + z25[s]) + c2[i] * (z11[s] + z26[s])) +
                                                  exp( -((exp(b[s] * t2[i]) - 1) / b[s]) * phi0[s] * phi1[s] * phi2[s] * phi3[s] * phi4[s] * exp(b[s] * (x2[i] - 77.5)) ) * ((phi0[s] * phi1[s] * phi2[s] * phi3[s] * phi4[s] * exp(b[s] * (x2[i] - 77.5 + t2[i])))^d2[i]) * (c0[i] * (z12[s] + z27[s]) + c1[i] * (z13[s] + z28[s]) + c2[i] * (z14[s] + z29[s]))) / (c0[i] * (z0[s] + z3[s] + z6[s] + z9[s] + z12[s] + z15[s] + z18[s] + z21[s] + z24[s] + z27[s]) + c1[i] * (z1[s] + z4[s] + z7[s] + z10[s] + z13[s] + z16[s] + z19[s] + z22[s] + z25[s] + z28[s]) + c2[i] * (z2[s] + z5[s] + z8[s] + z11[s] + z14[s] + z17[s] + z20[s] + z23[s] + z26[s] + z29[s]))  ) 
  }
  
}

log_lik_J4 <- log_lik_J4_P1 + log_lik_J4_P2


log_lik_J0_dev <- log_lik_J0
log_lik_J0_dens <- density(- 2 *log_lik_J0_dev)
log_lik_J1_dev <- log_lik_J1
log_lik_J1_dens <- density(- 2 *log_lik_J1_dev)
log_lik_J2_dev <- log_lik_J2
log_lik_J2_dens <- density(- 2 *log_lik_J2_dev)
log_lik_J3_dev <- log_lik_J3
log_lik_J3_dens <- density(- 2 *log_lik_J3_dev)
log_lik_J4_dev <- log_lik_J4
log_lik_J4_dens <- density(- 2 *log_lik_J4_dev)

plot(log_lik_J0_dens, type = "l", xlim=c(41075, 41150), ylim=c(0, 0.1))

par(mfrow=c(1,1), mai=c(0.9,0.9, 0.2,0.2))
plot(log_lik_J1_dens, type="l", lty=1, xlim=c(41095, 41145), ylim=c(0, 0.09), main="", xlab = "-2 log L", lwd=2)
lines(log_lik_J2_dens, type="l", lty=2, lwd=2)
lines(log_lik_J3_dens, type="l", lty=3, lwd=2)
lines(log_lik_J4_dens, type="l", lty=6, lwd=2)
legend('topright', c("J = 1","J = 2", "J = 3", "J = 4"), lty=c(1,2,3,6), bty='n', cex=0.9, lwd=2)


####################################################### - Allocation map - #########################################################

Alloc_table1 <- matrix(NA, nrow=9000, ncol=nrow(DS1))

for(i in 1:9000){
  
  phi0_1 <- par_postJ2_gn1_th[i,1]
  beta_1 <- par_postJ2_gn1_th[i,2]
  phi1_1 <- par_postJ2_gn1_th[i,3]
  phi2_1 <- par_postJ2_gn1_th[i,4]
  
  zeta0_1 <- par_postJ2_gn1_th[i,5]
  zeta1_1 <- par_postJ2_gn1_th[i,6]
  zeta2_1 <- par_postJ2_gn1_th[i,7]
  zeta3_1 <- par_postJ2_gn1_th[i,8]
  zeta4_1 <- par_postJ2_gn1_th[i,9]
  zeta5_1 <- par_postJ2_gn1_th[i,10]
  zeta6_1 <- par_postJ2_gn1_th[i,11]
  zeta7_1 <- par_postJ2_gn1_th[i,12]
  zeta8_1 <- par_postJ2_gn1_th[i,13]
  zeta9_1 <- par_postJ2_gn1_th[i,14]
  zeta10_1 <- par_postJ2_gn1_th[i,15]
  zeta11_1 <- par_postJ2_gn1_th[i,16]
  zeta12_1 <- par_postJ2_gn1_th[i,17]
  zeta13_1 <- par_postJ2_gn1_th[i,18]
  zeta14_1 <- par_postJ2_gn1_th[i,19]
  zeta15_1 <- par_postJ2_gn1_th[i,20]
  zeta16_1 <- par_postJ2_gn1_th[i,21]
  zeta17_1 <- par_postJ2_gn1_th[i,22]
  
  # - Dataset augmentation - Step 1
  
  # - Creation fractional weights
  
  DS1$PostG0C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta0_1 * bL + zeta9_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG0C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta1_1 * bL + zeta10_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG0C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * (zeta2_1 * bL + zeta11_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG1C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta3_1 * bL + zeta12_1  * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG1C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta4_1 * bL + zeta13_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG1C2 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * (zeta5_1 * bL + zeta14_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG2C0 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta6_1 * bL + zeta15_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  DS1$PostG2C1 <- exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * (zeta7_1 * bL + zeta16_1 * bH) / 
                ( exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x1 - 77.5))))                            * ((zeta0_1 + zeta1_1 + zeta2_1) * bL + (zeta9_1  + zeta10_1 + zeta11_1) * bH) + 
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x1 - 77.5)))) * ( phi1_1           ^ d1) * ((zeta3_1 + zeta4_1 + zeta5_1) * bL + (zeta12_1 + zeta13_1 + zeta14_1) * bH) +
                  exp(-(((exp(beta_1 * t1) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x1 - 77.5)))) * ((phi1_1 * phi2_1) ^ d1) * ((zeta6_1 + zeta7_1 + zeta8_1) * bL + (zeta15_1 + zeta16_1 + zeta17_1) * bH) )
  
  
  DS1$PostG2C2 <- 1 - DS1$PostG0C0 - DS1$PostG0C1 - DS1$PostG0C2 - DS1$PostG1C0 - DS1$PostG1C1 - DS1$PostG1C2 - DS1$PostG2C0 - DS1$PostG2C1
  
  

  # Draw a value for C,G
  
  ## - DS1
  g_unif <- runif(nrow(DS1), 0, 1)
  
  DS1$G <- ifelse((g_unif <= (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2) ), 0, ifelse(((g_unif > (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2)) & (g_unif <= (DS1$PostG0C0 + DS1$PostG0C1 + DS1$PostG0C2 + DS1$PostG1C0 + DS1$PostG1C1 + DS1$PostG1C2))), 1, 2))
  
  Alloc_table1[i, c(1:nrow(DS1))] <- DS1$G

}

Alloc_table2 <- matrix(NA, nrow=9000, ncol=nrow(DS2))

for(i in 1:9000){
  
  phi0_1 <- par_postJ2_gn1_th[i,1]
  beta_1 <- par_postJ2_gn1_th[i,2]
  phi1_1 <- par_postJ2_gn1_th[i,3]
  phi2_1 <- par_postJ2_gn1_th[i,4]
  
  zeta0_1 <- par_postJ2_gn1_th[i,5]
  zeta1_1 <- par_postJ2_gn1_th[i,6]
  zeta2_1 <- par_postJ2_gn1_th[i,7]
  zeta3_1 <- par_postJ2_gn1_th[i,8]
  zeta4_1 <- par_postJ2_gn1_th[i,9]
  zeta5_1 <- par_postJ2_gn1_th[i,10]
  zeta6_1 <- par_postJ2_gn1_th[i,11]
  zeta7_1 <- par_postJ2_gn1_th[i,12]
  zeta8_1 <- par_postJ2_gn1_th[i,13]
  zeta9_1 <- par_postJ2_gn1_th[i,14]
  zeta10_1 <- par_postJ2_gn1_th[i,15]
  zeta11_1 <- par_postJ2_gn1_th[i,16]
  zeta12_1 <- par_postJ2_gn1_th[i,17]
  zeta13_1 <- par_postJ2_gn1_th[i,18]
  zeta14_1 <- par_postJ2_gn1_th[i,19]
  zeta15_1 <- par_postJ2_gn1_th[i,20]
  zeta16_1 <- par_postJ2_gn1_th[i,21]
  zeta17_1 <- par_postJ2_gn1_th[i,22]
  
  # - Dataset augmentation - Step 1
  
  DS2$PostG0BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta0_1 * c0 + zeta1_1 * c1 + zeta2_1 * c2) / 
                ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1 ) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG0BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * (zeta9_1 * c0 + zeta10_1 * c1 + zeta11_1 * c2) / 
                ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1 ) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta3_1 * c0 + zeta4_1 * c1 + zeta5_1 * c2) / 
                ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1 ) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG1BH <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * (zeta12_1 * c0 + zeta13_1 * c1 + zeta14_1 * c2) / 
                ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1 ) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BL <- exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * (zeta6_1 * c0 + zeta7_1 * c1 + zeta8_1 * c2) / 
                ( exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1                   * exp(beta_1 * (x2 - 77.5))))                            * ((zeta0_1 + zeta9_1 ) * c0 + (zeta1_1 + zeta10_1) * c1 + (zeta2_1 + zeta11_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1          * exp(beta_1 * (x2 - 77.5)))) *  (phi1_1           ^ d2) * ((zeta3_1 + zeta12_1) * c0 + (zeta4_1 + zeta13_1) * c1 + (zeta5_1 + zeta14_1) * c2) +
                  exp(-(((exp(beta_1 * t2) - 1) / beta_1) * phi0_1 * phi1_1 * phi2_1 * exp(beta_1 * (x2 - 77.5)))) * ((phi1_1 * phi2_1) ^ d2) * ((zeta6_1 + zeta15_1) * c0 + (zeta7_1 + zeta16_1) * c1 + (zeta8_1 + zeta17_1) * c2) )
  
  DS2$PostG2BH <- 1 - DS2$PostG0BL - DS2$PostG0BH - DS2$PostG1BL - DS2$PostG1BH - DS2$PostG2BL
  

  # Draw a value for B,G
  
  ## - DS2
  g_unif <- runif(nrow(DS2), 0, 1)
  
  DS2$G <- ifelse((g_unif <= DS2$PostG0BL+ DS2$PostG0BH), 0, ifelse(((g_unif > DS2$PostG0BL+ DS2$PostG0BH) & (g_unif <= DS2$PostG0BL+ DS2$PostG0BH + DS2$PostG1BL+ DS2$PostG1BH)), 1, 2))
  
  Alloc_table2[i, c(1:nrow(DS2))] <- DS2$G

}

#plot()


Alloc_table1 <- as.data.frame(Alloc_table1)
Alloc_table2 <- as.data.frame(Alloc_table2)

alloc_mean <- rep(0, nrow(DS_C))
alloc_mean[c(1:nrow(DS1))] <- colMeans(Alloc_table1)
alloc_mean[c((nrow(DS1)+1):nrow(DS_C))] <- colMeans(Alloc_table2)

#remove(Alloc_table2)

### - Calc Expected allocation

phi0_p <- par_post_meanJ2_C11[1]
beta_p <- par_post_meanJ2_C11[2]
phi1_p <- par_post_meanJ2_C11[3]
phi2_p <- par_post_meanJ2_C11[4]
z0_p <- par_post_meanJ2_C11[5]
z1_p <- par_post_meanJ2_C11[6]
z2_p <- par_post_meanJ2_C11[7]
z3_p <- par_post_meanJ2_C11[8]
z4_p <- par_post_meanJ2_C11[9]
z5_p <- par_post_meanJ2_C11[10]
z6_p <- par_post_meanJ2_C11[11]
z7_p <- par_post_meanJ2_C11[12]
z8_p <- par_post_meanJ2_C11[13]
z9_p <- par_post_meanJ2_C11[14]
z10_p <- par_post_meanJ2_C11[15]
z11_p <- par_post_meanJ2_C11[16]
z12_p <- par_post_meanJ2_C11[17]
z13_p <- par_post_meanJ2_C11[18]
z14_p <- par_post_meanJ2_C11[19]
z15_p <- par_post_meanJ2_C11[20]
z16_p <- par_post_meanJ2_C11[21]
z17_p <- par_post_meanJ2_C11[22]



################################################## - Allocation vs Expectation - ##################################################


alloc_exp <- rep(0, nrow(DS_C))

alloc_exp[c(1:nrow(DS1))] <- ((bL * (z3_p + z4_p + z5_p) + bH * (z12_p + z13_p + z14_p))                     * exp( -((exp(beta_p * t1) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (x1 - 77.5)) ) * ((phi0_p * phi1_p          * exp(beta_p * (x1 - 77.5 + t1)))^d1) + 
                          2 * (bL * (z6_p + z7_p + z8_p) + bH * (z15_p + z16_p + z17_p))                     * exp( -((exp(beta_p * t1) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (x1 - 77.5)) ) * ((phi0_p * phi1_p * phi2_p * exp(beta_p * (x1 - 77.5 + t1)))^d1)) /
                            ( (bL * (z0_p + z1_p + z2_p) + bH * (z9_p  + z10_p + z11_p))                     * exp( -((exp(beta_p * t1) - 1) / beta_p) * phi0_p                   * exp(beta_p * (x1 - 77.5)) ) * ((phi0_p                   * exp(beta_p * (x1 - 77.5 + t1)))^d1) +
                              (bL * (z3_p + z4_p + z5_p) + bH * (z12_p + z13_p + z14_p))                     * exp( -((exp(beta_p * t1) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (x1 - 77.5)) ) * ((phi0_p * phi1_p          * exp(beta_p * (x1 - 77.5 + t1)))^d1) + 
                              (bL * (z6_p + z7_p + z8_p) + bH * (z15_p + z16_p + z17_p))                     * exp( -((exp(beta_p * t1) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (x1 - 77.5)) ) * ((phi0_p * phi1_p * phi2_p * exp(beta_p * (x1 - 77.5 + t1)))^d1) )

alloc_exp[c((nrow(DS1)+1):nrow(DS_C))] <- ((c0 * (z3_p + z12_p) + c1 * (z4_p + z13_p) + c2 * (z5_p + z14_p)) * exp( -((exp(beta_p * t2) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (x2 - 77.5)) ) * ((phi0_p * phi1_p          * exp(beta_p * (x2 - 77.5 + t2)))^d2) + 
                                       2 * (c0 * (z6_p + z15_p) + c1 * (z7_p + z16_p) + c2 * (z8_p + z17_p)) * exp( -((exp(beta_p * t2) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (x2 - 77.5)) ) * ((phi0_p * phi1_p * phi2_p * exp(beta_p * (x2 - 77.5 + t2)))^d2)) /
                                         ( (c0 * (z0_p + z9_p ) + c1 * (z1_p + z10_p) + c2 * (z2_p + z11_p)) * exp( -((exp(beta_p * t2) - 1) / beta_p) * phi0_p                   * exp(beta_p * (x2 - 77.5)) ) * ((phi0_p                   * exp(beta_p * (x2 - 77.5 + t2)))^d2) + 
                                           (c0 * (z3_p + z12_p) + c1 * (z4_p + z13_p) + c2 * (z5_p + z14_p)) * exp( -((exp(beta_p * t2) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (x2 - 77.5)) ) * ((phi0_p * phi1_p          * exp(beta_p * (x2 - 77.5 + t2)))^d2) + 
                                           (c0 * (z6_p + z15_p) + c1 * (z7_p + z16_p) + c2 * (z8_p + z17_p)) * exp( -((exp(beta_p * t2) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (x2 - 77.5)) ) * ((phi0_p * phi1_p * phi2_p * exp(beta_p * (x2 - 77.5 + t2)))^d2)   )

par(mfrow=c(1,1), mai=c(0.9,0.9, 0.2,0.2)) # save as 8.5 x 6.5
plot(alloc_exp[c(1:nrow(DS1))], alloc_mean[c(1:nrow(DS1))], pch=19, cex=0.1, xlim=c(0,2), ylim=c(0,2), xlab="Expected Allocation", ylab="Mean Allocation", main="")
abline(a=0, b=1, lty=2)
points(alloc_exp[c((nrow(DS1)+1):nrow(DS_C))], alloc_mean[c((nrow(DS1)+1):nrow(DS_C))], pch=19, cex=0.1, col="grey")
#points(alloc_exp[c(1:nrow(DS1))], alloc_mean[c(1:nrow(DS1))], pch=19, cex=0.1)
legend('topleft', c("Pension Scheme P1","Pension Scheme P2"), pch=19, cex=1, bty='n', col=c("black", "grey"))


alloc <- matrix(NA, nrow(DS_C), 2)
alloc[,1] <- alloc_exp
alloc[,2] <- alloc_mean

colnames(alloc) <- c("exp", "mean")

alloc <- as.data.frame(alloc)
alloc <- alloc[order(alloc$exp),]

cor_matrix <- cor(data.frame(par_postJ3_gn1_th))






#################################################### - Cumsums - ###############################################

partl_diff <- matrix(NA, 9000, 22)

for(i in 1:9000){
  partl_diff[i,] <- colSums(par_postJ2_gn1_th[c(1:i),]) - i * par_post_meanJ2_C11
}


plot(c(1:9000), partl_diff[,1], type="l", ylim=c(-5,5))
lines(c(1:9000), partl_diff[,2], type="l")
lines(c(1:9000), partl_diff[,3], type="l")
lines(c(1:9000), partl_diff[,4], type="l")
lines(c(1:9000), partl_diff[,5], type="l")
lines(c(1:9000), partl_diff[,6], type="l")
lines(c(1:9000), partl_diff[,7], type="l")
lines(c(1:9000), partl_diff[,8], type="l")
lines(c(1:9000), partl_diff[,9], type="l")
lines(c(1:9000), partl_diff[,10], type="l")
lines(c(1:9000), partl_diff[,11], type="l")
lines(c(1:9000), partl_diff[,12], type="l")
lines(c(1:9000), partl_diff[,13], type="l")
lines(c(1:9000), partl_diff[,14], type="l")
lines(c(1:9000), partl_diff[,15], type="l")
lines(c(1:9000), partl_diff[,16], type="l")
lines(c(1:9000), partl_diff[,17], type="l")
lines(c(1:9000), partl_diff[,18], type="l")
lines(c(1:9000), partl_diff[,19], type="l")
lines(c(1:9000), partl_diff[,20], type="l")
lines(c(1:9000), partl_diff[,21], type="l")
lines(c(1:9000), partl_diff[,22], type="l")

#################################### - Kolmogorov-Smirnov - ##################################################

kolmogorov_smirnov_table <- matrix(NA, 22, 2)
rownames(kolmogorov_smirnov_table) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")
colnames(kolmogorov_smirnov_table) <- c("Statistic", "p-value")

for(i in 1:22){
kolm_smirnov <- ks.test(par_postJ2_gn1_th[c(1:4500),i],  par_postJ2_gn1_th[c(4501:9000),i], exact = NULL, alternative="two.sided")
kolmogorov_smirnov_table[i,1] <- kolm_smirnov$statistic
kolmogorov_smirnov_table[i,2] <- kolm_smirnov$p.value
}

# - Analytical computation
ks_sum <- function(k, par){
  ks_sum_bar <- abs(sum(ifelse((par_postJ2_gn1_th[c(1:4500), par] < k), 1, 0)) - sum(ifelse((par_postJ2_gn1_th[c(4501:9000), par] < k), 1, 0)))
  return(ks_sum_bar)
}

ks_sum_post_par <- rep(NA, 22)

ks_sum_processing <- matrix(NA, nrow = 10000, ncol=2)
ks_sum_processing[,1] <- seq(0.0001, 1, by=0.0001)

for(i in 1:22){
  
  ks_sum_processing <- as.matrix(ks_sum_processing)
  
  for(j in 1:10000){
    
    ks_sum_processing[j,2] <- ks_sum(ks_sum_processing[j,1],i)
    
  }
  
  ks_sum_processing <- as.data.frame(ks_sum_processing)
  ks_sum_post_par[i] <- ks_sum_processing[which.max(ks_sum_processing[,2]),2]/4500
  
}

sqrt(4500) * ks_sum_post_par

############################################### - Gelman-Rubin - #############################################

R_hat <- rep(NA, 22)

for(par in 1:22){
psi_par_1 <- par_postJ2_gn1_th[c(1:4500),par]
psi_par_2 <- par_postJ2_gn1_th[c(4501:9000),par]
psi_par_3 <- par_postJ2_gn2_th[c(1:4500),par]
psi_par_4 <- par_postJ2_gn2_th[c(4501:9000),par]

psi_d_bar <- (sum(psi_par_1) + sum(psi_par_2) + sum(psi_par_3) + sum(psi_par_4)) / 18000

m=4
n=4500

B <- (n / (m-1)) * (((mean(psi_par_1) - psi_d_bar)^2) + ((mean(psi_par_2) - psi_d_bar)^2) + ((mean(psi_par_3) - psi_d_bar)^2) + ((mean(psi_par_4) - psi_d_bar)^2))
W <- (1 / m) * (var(psi_par_1) + var(psi_par_2) + var(psi_par_3) + var(psi_par_4)) * (n / (n-1))

var_post <- ((n - 1) / n) * W + B / n

R_hat[par] <- sqrt(var_post / W)

}

# - Split of chain in 3 parts

R_hat <- rep(NA, 22)

for(par in 1:22){
  psi_par_1 <- par_postJ2_gn1_th[c(1:3000),par]
  psi_par_2 <- par_postJ2_gn1_th[c(3001:6000),par]
  psi_par_3 <- par_postJ2_gn1_th[c(6001:9000),par]
  psi_par_4 <- par_postJ2_gn2_th[c(1:3000),par]
  psi_par_5 <- par_postJ2_gn2_th[c(3001:6000),par]
  psi_par_6 <- par_postJ2_gn2_th[c(6001:9000),par]
  
  psi_d_bar <- (sum(psi_par_1) + sum(psi_par_2) + sum(psi_par_3) + sum(psi_par_4) + sum(psi_par_5) + sum(psi_par_6)) / 18000
  
  m=6
  n=3000
  
  B <- (n / (m-1)) * (((mean(psi_par_1) - psi_d_bar)^2) + ((mean(psi_par_2) - psi_d_bar)^2) + ((mean(psi_par_3) - psi_d_bar)^2) + ((mean(psi_par_4) - psi_d_bar)^2) + ((mean(psi_par_5) - psi_d_bar)^2) + ((mean(psi_par_6) - psi_d_bar)^2))
  W <- (1 / m) * (var(psi_par_1) + var(psi_par_2) + var(psi_par_3) + var(psi_par_4) + var(psi_par_5) + var(psi_par_6)) * (n / (n-1))
  
  var_post <- ((n - 1) / n) * W + B / n
  
  R_hat[par] <- sqrt(var_post / W)
  
}

##################### - Ergodic mean - ############################
par_postJ2_gn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ2_gibbs_norm1))
par_postJ2_gn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ2_gibbs_norm2))
par_postJ2_mn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ2_metr_norm1))
par_postJ2_mn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ2_metr_norm2))

par_postJ2_gn1_erg_m <- as.data.frame(par_postJ2_gn1_erg_m)
par_postJ2_gn2_erg_m <- as.data.frame(par_postJ2_gn2_erg_m)
par_postJ2_mn1_erg_m <- as.data.frame(par_postJ2_mn1_erg_m)
par_postJ2_mn2_erg_m <- as.data.frame(par_postJ2_mn2_erg_m)

#################### - ergodic low post var - ######################
par_postJ2_gn1_erg_lb <- matrix(NA, nrow=100000, ncol=ncol(par_postJ2_gibbs_norm1))
par_postJ2_gn1_erg_lb <- as.data.frame(par_postJ2_gn1_erg_lb)

#################### - ergodic upp post var - ######################
par_postJ2_gn1_erg_ub <- matrix(NA, nrow=100000, ncol=ncol(par_postJ2_gibbs_norm1))
par_postJ2_gn1_erg_ub <- as.data.frame(par_postJ2_gn1_erg_ub)

for (i in 1:100000){
  #par_postJ2_gn1_erg_m[i,] <- colSums(par_postJ2_gibbs_norm1[c(1:i), , drop = FALSE])/i
  par_postJ2_gn1_erg_lb[i,] <- colQuantiles(as.matrix(par_postJ2_gibbs_norm1[c(1:i),]), rows = NULL, cols = NULL, probs = 0.025)
  par_postJ2_gn1_erg_ub[i,] <- colQuantiles(as.matrix(par_postJ2_gibbs_norm1[c(1:i),]), rows = NULL, cols = NULL, probs = 0.975)
  #par_postJ2_gn2_erg_m[i,] <- colSums(par_postJ2_gibbs_norm2[c(1:i), , drop = FALSE])/i
  #par_postJ2_mn1_erg_m[i,] <- colSums(par_postJ2_metr_norm1[c(1:i), , drop = FALSE])/i
  #par_postJ2_mn2_erg_m[i,] <- colSums(par_postJ2_metr_norm2[c(1:i), , drop = FALSE])/i
}

#### - Plots for paper - #####
par(mfrow=c(4,3), mai=c(0.6,0.5, 0.2,0.2))
## - Ergodic plots

plot(1:100000, log(par_postJ2_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha", ylim=c(-3.5, -2.5)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(log(par_post_meanJ2_C11[1]), 100000), type="l", lty=2)

plot(c(1:100000), log(par_postJ2_gn1_erg_m[,1]), type="l", xlab="Iteration", ylab="alpha", ylim=c(-3.5, -2.5))
lines(c(1:100000), log(par_postJ2_gn1_erg_lb[,1]), type="l", lty=3)
lines(c(1:100000), log(par_postJ2_gn1_erg_ub[,1]), type="l", lty=3)
lines(c(1:100000), rep(par_post_meanJ2_C11[1], 100000), type="l", lty=2)

plot(alpha_J2_gn1, xlab = "alpha", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(log(par_postJ2_gn1_th[,1]), probs = 0.025), lty=3)
abline(v=quantile(log(par_postJ2_gn1_th[,1]), probs = 0.975), lty=3)
abline(v=log(mean(par_postJ2_gn1_th[,1])), lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim=c(0.07, 0.15)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(par_post_meanJ2_C11[2], 100000), type="l", lty=2)

plot(c(1:100000), par_postJ2_gn1_erg_m[,2], type="l", xlab="Iteration", ylab="beta", ylim=c(0.07,0.15))
lines(c(1:100000), par_postJ2_gn1_erg_lb[,2], type="l", lty=3)
lines(c(1:100000), par_postJ2_gn1_erg_ub[,2], type="l", lty=3)
lines(c(1:100000), rep(par_post_meanJ2_C11[2], 100000), type="l", lty=2)

plot(beta_J2_gn1, xlab = "beta", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(par_postJ2_gn1_th[,2], probs = 0.025), lty=3)
abline(v=quantile(par_postJ2_gn1_th[,2], probs = 0.975), lty=3)
abline(v=mean(par_postJ2_gn1_th[,2]), lty=2)

plot(1:100000, log(par_postJ2_gibbs_norm1[,3]), type = "l", xlab="Iteration", ylab = "psi1", ylim=c(-1,0)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(log(par_post_meanJ2_C11[3]), 100000), type="l", lty=2)

plot(c(1:100000), log(par_postJ2_gn1_erg_m[,3]), type="l", xlab="Iteration", ylab="psi1", ylim=c(-1,0))
lines(c(1:100000), log(par_postJ2_gn1_erg_lb[,3]), type="l", lty=3)
lines(c(1:100000), log(par_postJ2_gn1_erg_ub[,3]), type="l", lty=3)
lines(c(1:100000), rep(par_post_meanJ2_C11[3], 100000), type="l", lty=2)

plot(gamma1_J2_gn1, xlab = "psi1", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(log(par_postJ2_gn1_th[,3]), probs = 0.025), lty=3)
abline(v=quantile(log(par_postJ2_gn1_th[,3]), probs = 0.975), lty=3)
abline(v=log(mean(par_postJ2_gn1_th[,3])), lty=2)

plot(1:100000, log(par_postJ2_gibbs_norm1[,4]), type = "l", xlab="Iteration", ylab = "psi2", ylim=c(-1,0)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(log(par_post_meanJ2_C11[4]), 100000), type="l", lty=2)

plot(c(1:100000), log(par_postJ2_gn1_erg_m[,4]), type="l", xlab="Iteration", ylab="psi2", ylim=c(-1,0))
lines(c(1:100000), log(par_postJ2_gn1_erg_lb[,4]), type="l", lty=3)
lines(c(1:100000), log(par_postJ2_gn1_erg_ub[,4]), type="l", lty=3)
lines(c(1:100000), rep(par_post_meanJ2_C11[4], 100000), type="l", lty=2)

plot(gamma2_J2_gn1, xlab = "psi2", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(log(par_postJ2_gn1_th[,4]), probs = 0.025), lty=3)
abline(v=quantile(log(par_postJ2_gn1_th[,4]), probs = 0.975), lty=3)
abline(v=log(mean(par_postJ2_gn1_th[,4])), lty=2)

# - Traceplots
plot(1:100000, log(par_postJ2_gibbs_norm1[,1]), type = "l", xlab="Iteration", ylab = "alpha", ylim=c(-3.5, -2.5)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(log(par_post_meanJ2_C11[1]), 100000), type="l", lty=2)

plot(1:100000, par_postJ2_gibbs_norm1[,2], type = "l", xlab="Iteration", ylab = "beta", ylim=c(0.07, 0.15)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(par_post_meanJ2_C11[2], 100000), type="l", lty=2)

plot(1:100000, log(par_postJ2_gibbs_norm1[,3]), type = "l", xlab="Iteration", ylab = "psi1", ylim=c(-1,0)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(log(par_post_meanJ2_C11[3]), 100000), type="l", lty=2)

plot(1:100000, log(par_postJ2_gibbs_norm1[,4]), type = "l", xlab="Iteration", ylab = "psi2", ylim=c(-1,0)) #ylim=c(-3.2, -1), 
lines(c(1:100000), rep(log(par_post_meanJ2_C11[4]), 100000), type="l", lty=2)


# - Densities
plot(alpha_J2_gn1, xlab = "alpha", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(log(par_postJ2_gn1_th[,1]), probs = 0.025), lty=3)
abline(v=quantile(log(par_postJ2_gn1_th[,1]), probs = 0.975), lty=3)
abline(v=log(mean(par_postJ2_gn1_th[,1])), lty=2)

plot(beta_J2_gn1, xlab = "beta", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(par_postJ2_gn1_th[,2], probs = 0.025), lty=3)
abline(v=quantile(par_postJ2_gn1_th[,2], probs = 0.975), lty=3)
abline(v=mean(par_postJ2_gn1_th[,2]), lty=2)

plot(gamma1_J2_gn1, xlab = "psi1", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(log(par_postJ2_gn1_th[,3]), probs = 0.025), lty=3)
abline(v=quantile(log(par_postJ2_gn1_th[,3]), probs = 0.975), lty=3)
abline(v=log(mean(par_postJ2_gn1_th[,3])), lty=2)

plot(gamma2_J2_gn1, xlab = "psi2", main=" ") #, main="Density of the posterior distribution of alpha"
abline(v=quantile(log(par_postJ2_gn1_th[,4]), probs = 0.025), lty=3)
abline(v=quantile(log(par_postJ2_gn1_th[,4]), probs = 0.975), lty=3)
abline(v=log(mean(par_postJ2_gn1_th[,4])), lty=2)


####################

plot(c(1:100000), par_postJ2_gn1_erg_m[,1], type="l") #, ylim=c(0,1)
plot(c(1:100000), par_postJ2_gn1_erg_m[,2], type="l") #, ylim=c(0,1)
plot(c(1:100000), par_postJ2_gn1_erg_m[,3], type="l", ylim=c(0,1)) #, ylim=c(0,1)
lines(c(1:100000), par_postJ2_gn1_erg_m[,4], type="l", col="red") #, ylim=c(0,1)

plot(c(1:100000), par_postJ2_gn1_erg_m[,5], type="l", ylim=c(0,0.25))

# - Egodic mean zeta - Save as 8 x 6
par(mfrow=c(1,1), mai=c(0.8,0.8, 0.2,0.2))
plot(c(1:100000), par_postJ2_gn1_erg_m[,5], type="l", ylim=c(0,0.27), xlab="Iteration", ylab="zeta")

for(ind in 6:22){
  lines(c(1:100000), par_postJ2_gn1_erg_m[,ind], type="l", ylim=c(0,0.27))
}
text(locator(), labels = c("zeta4", "zeta0", "zeta1"))

##################################### - Analysis of the hazard function - ###############################################

##################### - Ergodic means for thesis - ############################

# - J=0

par_postJ0_gn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ0_gibbs_norm1))
par_postJ0_gn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ0_gibbs_norm2))
par_postJ0_mn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ0_metr_norm1))
par_postJ0_mn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ0_metr_norm2))

par_postJ0_gn1_erg_m <- as.data.frame(par_postJ0_gn1_erg_m)
par_postJ0_gn2_erg_m <- as.data.frame(par_postJ0_gn2_erg_m)
par_postJ0_mn1_erg_m <- as.data.frame(par_postJ0_mn1_erg_m)
par_postJ0_mn2_erg_m <- as.data.frame(par_postJ0_mn2_erg_m)

for (i in 1:100000){
  par_postJ0_gn1_erg_m[i,] <- colSums(par_postJ0_gibbs_norm1[c(1:i), , drop = FALSE])/i
  par_postJ0_gn2_erg_m[i,] <- colSums(par_postJ0_gibbs_norm2[c(1:i), , drop = FALSE])/i
  par_postJ0_mn1_erg_m[i,] <- colSums(par_postJ0_metr_norm1[c(1:i), , drop = FALSE])/i
  par_postJ0_mn2_erg_m[i,] <- colSums(par_postJ0_metr_norm2[c(1:i), , drop = FALSE])/i
}

# - J=1

par_postJ1_gn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ1_gibbs_norm1))
par_postJ1_gn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ1_gibbs_norm2))
par_postJ1_mn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ1_metr_norm1))
par_postJ1_mn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ1_metr_norm2))

par_postJ1_gn1_erg_m <- as.data.frame(par_postJ1_gn1_erg_m)
par_postJ1_gn2_erg_m <- as.data.frame(par_postJ1_gn2_erg_m)
par_postJ1_mn1_erg_m <- as.data.frame(par_postJ1_mn1_erg_m)
par_postJ1_mn2_erg_m <- as.data.frame(par_postJ1_mn2_erg_m)

for (i in 1:100000){
  par_postJ1_gn1_erg_m[i,] <- colSums(par_postJ1_gibbs_norm1[c(1:i), , drop = FALSE])/i
  par_postJ1_gn2_erg_m[i,] <- colSums(par_postJ1_gibbs_norm2[c(1:i), , drop = FALSE])/i
  par_postJ1_mn1_erg_m[i,] <- colSums(par_postJ1_metr_norm1[c(1:i), , drop = FALSE])/i
  par_postJ1_mn2_erg_m[i,] <- colSums(par_postJ1_metr_norm2[c(1:i), , drop = FALSE])/i
}

# - J=3

par_postJ3_gn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ3_gibbs_norm1))
par_postJ3_gn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ3_gibbs_norm2))
par_postJ3_mn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ3_metr_norm1))
par_postJ3_mn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ3_metr_norm2))

par_postJ3_gn1_erg_m <- as.data.frame(par_postJ3_gn1_erg_m)
par_postJ3_gn2_erg_m <- as.data.frame(par_postJ3_gn2_erg_m)
par_postJ3_mn1_erg_m <- as.data.frame(par_postJ3_mn1_erg_m)
par_postJ3_mn2_erg_m <- as.data.frame(par_postJ3_mn2_erg_m)

for (i in 1:100000){
  par_postJ3_gn1_erg_m[i,] <- colSums(par_postJ3_gibbs_norm1[c(1:i), , drop = FALSE])/i
  par_postJ3_gn2_erg_m[i,] <- colSums(par_postJ3_gibbs_norm2[c(1:i), , drop = FALSE])/i
  par_postJ3_mn1_erg_m[i,] <- colSums(par_postJ3_metr_norm1[c(1:i), , drop = FALSE])/i
  par_postJ3_mn2_erg_m[i,] <- colSums(par_postJ3_metr_norm2[c(1:i), , drop = FALSE])/i
}

# - J=4

par_postJ4_gn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ4_gibbs_norm1))
par_postJ4_gn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ4_gibbs_norm2))
par_postJ4_mn1_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ4_metr_norm1))
par_postJ4_mn2_erg_m <- matrix(NA, nrow=100000, ncol=ncol(par_postJ4_metr_norm2))

par_postJ4_gn1_erg_m <- as.data.frame(par_postJ4_gn1_erg_m)
par_postJ4_gn2_erg_m <- as.data.frame(par_postJ4_gn2_erg_m)
par_postJ4_mn1_erg_m <- as.data.frame(par_postJ4_mn1_erg_m)
par_postJ4_mn2_erg_m <- as.data.frame(par_postJ4_mn2_erg_m)

for (i in 1:100000){
  par_postJ4_gn1_erg_m[i,] <- colSums(par_postJ4_gibbs_norm1[c(1:i), , drop = FALSE])/i
  par_postJ4_gn2_erg_m[i,] <- colSums(par_postJ4_gibbs_norm2[c(1:i), , drop = FALSE])/i
  par_postJ4_mn1_erg_m[i,] <- colSums(par_postJ4_metr_norm1[c(1:i), , drop = FALSE])/i
  par_postJ4_mn2_erg_m[i,] <- colSums(par_postJ4_metr_norm2[c(1:i), , drop = FALSE])/i
}
