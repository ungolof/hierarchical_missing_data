#############################################

# Datasets loading

#############################################

#Read dataset
wd <- setwd("H:\\.windows_settings\\Desktop\\Survival Analysis\\R workflow")

#wd_home <- setwd("/Users/francescou/Desktop/Survival Analysis/R workflow/")

SurvIn1<-read.csv(file="SurvivalInput_41751.csv", header=TRUE)   #nrow(SurvIn[SurvIn$Status==1,])
SurvIn1$Gender<-ifelse(SurvIn1$Gender=="M", 1, 0)
SurvIn1$Status<-ifelse(SurvIn1$Status=="DEATH", 1, 0)

SurvIn2<-read.csv(file="SurvivalInput_41755.csv", header=TRUE)   #nrow(SurvIn[SurvIn$Status==1,])
SurvIn2$Gender<-ifelse(SurvIn2$Gender=="M", 1, 0)
SurvIn2$Status<-ifelse(SurvIn2$Status=="DEATH", 1, 0)
SurvIn2<-as.data.frame(SurvIn2)
SurvIn2$Reserve <- NULL    #erase column

SurvIn1 <- unique(SurvIn1)
SurvIn2 <- unique(SurvIn2)

library(xtable)
######################################

# Dataset definition

######################################

#Verify deaths occurrence
###SurvIn1$EndObs<-SurvIn1$EntryYear+SurvIn1$TimeObserved
###SurvIn2$EndObs<-SurvIn2$EntryYear+SurvIn2$TimeObserved

SurvIn <- merge(SurvIn1, SurvIn2, by=c("EntryYear", "EntryAge", "TimeObserved", "EntryDuration", "Gender", "Status"), all=TRUE)

SurvIn$Benefit.x[is.na(SurvIn$Benefit.x)] <- -1
SurvIn$Benefit.y[is.na(SurvIn$Benefit.y)] <- -1

#Considering the eventuality that the common individual between these two datasets gets more benefits
SurvIn[,"Benefit"] <- apply(SurvIn[,c(7,8)], 1, max)

#Elimination of duplicated rows
SurvIn <- SurvIn[,c(1,2,3,4,5,6,9,10)]
SurvIn <- unique(SurvIn)
#####################################

# - Initial characteristics of the working dataset

## - Number of males and females
nrow(SurvIn[SurvIn$Gender==0,])

## - Total time observed
sum(SurvIn$TimeObserved)
nrow(SurvIn[(SurvIn$EntryYear+SurvIn$TimeObserved)<2000,])

## - Observed deaths
nrow(SurvIn[SurvIn$Status==1,])
nrow(SurvIn[SurvIn$Status==1 & SurvIn$Gender==0,])

## - Benefit amount histogram
hist(SurvIn$Benefit, breaks=500, main="Histogram of Benefit distribution", xlab="Benefit amount")

# Dataset cleansing

#####################################

#Erasing Exposures under 60 and adjusting the time observed, the entry year and the entry duration, erase exposures after 31/12/2009, and consider
#just male population

##Erase exposures under 60
StartingAge<-60
SurvIn <- SurvIn[!((SurvIn$EntryAge+SurvIn$TimeObserved) < StartingAge),]

##Variable adjustment
SurvIn$EntryYear[((SurvIn$EntryAge+SurvIn$TimeObserved) >= StartingAge & SurvIn$EntryAge < StartingAge)] <- SurvIn$EntryYear[((SurvIn$EntryAge+SurvIn$TimeObserved)>= StartingAge & SurvIn$EntryAge< StartingAge)] + StartingAge - SurvIn$EntryAge[((SurvIn$EntryAge+SurvIn$TimeObserved)>= StartingAge & SurvIn$EntryAge< StartingAge)]
SurvIn$EntryDuration[((SurvIn$EntryAge+SurvIn$TimeObserved) >= StartingAge & SurvIn$EntryAge < StartingAge)] <- SurvIn$EntryDuration[((SurvIn$EntryAge+SurvIn$TimeObserved) >= StartingAge & SurvIn$EntryAge < StartingAge)] + StartingAge - SurvIn$EntryAge[((SurvIn$EntryAge+SurvIn$TimeObserved)>= StartingAge & SurvIn$EntryAge < StartingAge)]

###Creation of two new variables TObs for TimeObserved and EntAge for EntryAge
SurvIn$TObs<-SurvIn$TimeObserved
SurvIn$EntAge<-SurvIn$EntryAge
SurvIn$TObs[((SurvIn$EntryAge+SurvIn$TimeObserved) >= StartingAge & SurvIn$EntryAge < StartingAge)] <- SurvIn$TimeObserved[((SurvIn$EntryAge+SurvIn$TimeObserved)>= StartingAge & SurvIn$EntryAge< StartingAge)] - StartingAge + SurvIn$EntryAge[((SurvIn$EntryAge+SurvIn$TimeObserved)>= StartingAge & SurvIn$EntryAge< StartingAge)]
SurvIn$EntAge[((SurvIn$EntryAge+SurvIn$TimeObserved) >= StartingAge & SurvIn$EntryAge < StartingAge)] <- StartingAge

##Erase exposures after 31/12/2009
SurvIn<-SurvIn[!SurvIn$EntryYear>=2010,]
SurvIn$TObs[(SurvIn$EntryYear+SurvIn$TObs)>=2010]<- 2010 - SurvIn$EntryYear[(SurvIn$EntryYear+SurvIn$TObs)>=2010]

#Indicator check for benefit size
######SurvIn$Ind<-ifelse(SurvIn$Benefit.x>=SurvIn$Benefit.y, 1, 0)

##########################################

# Further elaborations before starting

##########################################

##Erasing female population
SurvInM <- SurvIn[!SurvIn$Gender==0,]

#Creation of the holed datasets (then consider nested holes)
SurvIn$Group<-as.character(SurvIn$Group)
SurvInM$Group<-as.character(SurvInM$Group)

#Dataset with absent geodemographic profile
SurvIn$Group[is.na(SurvIn$Group)]<- -1
SurvInM$Group[is.na(SurvInM$Group)]<- -1

#Treat code 98 and 99 as missing
SurvIn$Group[SurvIn$Group=="98" | SurvIn$Group=="99"] <- -1
SurvInM$Group[SurvInM$Group=="98" | SurvInM$Group=="99"] <- -1


#Working Dataset (only complete observations)
SurvInMW <- SurvInM[SurvInM$Group != "-1",]

#####################################################################################################

# - Descriptive characteristics of the mortality experience

#####################################################################################################

sum(SurvInMW$TObs)
sum(SurvInMW$Status)
min(SurvInMW$EntryYear)
max(SurvInMW$EntryYear + SurvInMW$TObs)

sum(SurvInM$TObs)
sum(SurvInM$Status)
min(SurvInM$EntryYear)
max(SurvInM$EntryYear + SurvInM$TObs)

quantile(SurvInMW$Benefit, 0.7) 
nrow(SurvInMW[SurvInMW$Benefit < 8500,])/nrow(SurvInMW)

############################################################################################################

#Recoding - numbering for each geo-demographic

#Recoding - numbering for each geo-demographic
SurvInMW$Numbering <- ifelse(SurvInMW$Group=="A",1, ifelse(SurvInMW$Group=="B",2,ifelse(SurvInMW$Group=="C",3,ifelse(SurvInMW$Group=="D",4, ifelse(SurvInMW$Group=="E",5,ifelse(SurvInMW$Group=="F",6,
                                                                                                                                                                                ifelse(SurvInMW$Group=="G",7, ifelse(SurvInMW$Group=="H",8,ifelse(SurvInMW$Group=="I",9,ifelse(SurvInMW$Group=="J",10,ifelse(SurvInMW$Group=="K",11,ifelse(SurvInMW$Group=="L",12,ifelse(SurvInMW$Group=="M",13,ifelse(SurvInMW$Group=="N",14,
                                                                                                                                                                                                                                                                                                                                                                                                       ifelse(SurvInMW$Group=="O",15,16)) )) )))))))))))

##########################################################################

# - Complete data comparison

##########################################################################

# - 0 Creation of DS_C, the working dataset

# Complete dataset
DS_C <- SurvInMW[,c(1,6,7,8,9,10)]
DS_C$GeoGroup <- ifelse((DS_C$Group=="A" | DS_C$Group=="B" | DS_C$Group=="C" | DS_C$Group=="D" | DS_C$Group=="E"),2, ifelse((DS_C$Group=="O" | DS_C$Group=="J" | DS_C$Group=="K" | DS_C$Group=="90" | DS_C$Group=="91" | DS_C$Group=="92"),0,1))
DS_C$BL <- ifelse(DS_C$Benefit < 8500,1,0)
DS_C$BH <- 1 - DS_C$BL
DS_C$C0 <- ifelse(DS_C$GeoGroup==0, 1, 0)
DS_C$C1 <- ifelse(DS_C$GeoGroup==1, 1, 0)
DS_C$C2 <- ifelse(DS_C$GeoGroup==2, 1, 0)

# - Socio economic characteristics table
socio_econ <- matrix(NA, 4, 3)
socio_econ[1,1] <- nrow(DS_C[DS_C$C0==1 & DS_C$BL==1,])
socio_econ[1,2] <- nrow(DS_C[DS_C$C1==1 & DS_C$BL==1,])
socio_econ[1,3] <- nrow(DS_C[DS_C$C2==1 & DS_C$BL==1,])
socio_econ[3,1] <- nrow(DS_C[DS_C$C0==1 & DS_C$BL==0,])
socio_econ[3,2] <- nrow(DS_C[DS_C$C1==1 & DS_C$BL==0,])
socio_econ[3,3] <- nrow(DS_C[DS_C$C2==1 & DS_C$BL==0,])

socio_econ[2,1] <- nrow(DS_C[DS_C$C0==1 & DS_C$BL==1,])/nrow(DS_C)
socio_econ[2,2] <- nrow(DS_C[DS_C$C1==1 & DS_C$BL==1,])/nrow(DS_C)
socio_econ[2,3] <- nrow(DS_C[DS_C$C2==1 & DS_C$BL==1,])/nrow(DS_C)
socio_econ[4,1] <- nrow(DS_C[DS_C$C0==1 & DS_C$BL==0,])/nrow(DS_C)
socio_econ[4,2] <- nrow(DS_C[DS_C$C1==1 & DS_C$BL==0,])/nrow(DS_C)
socio_econ[4,3] <- nrow(DS_C[DS_C$C2==1 & DS_C$BL==0,])/nrow(DS_C)

library(truncnorm)
library(TruncatedDistributions)

t <- DS_C$TObs
x <- DS_C$EntAge
bL <- DS_C$BL
bH <- DS_C$BH
c0 <- DS_C$C0
c1 <- DS_C$C1
c2 <- DS_C$C2
d <- DS_C$Status

########## - Dataset split - ###########

split_percentage <- 0.5


DS_C <- DS_C[sample(nrow(DS_C)),]

DS1 <- DS_C[1:round(nrow(DS_C) * split_percentage),]
DS2 <- DS_C[(round(nrow(DS_C) * (split_percentage)) + 1):nrow(DS_C),]

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

##################################################### - J=1 - #################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00003, 0.095, 0.07)
init_zeta1 <- c(0.2, 0.23, 0.01, 0.01, 0.19, 0.12, 0.01, 0.01, 0.01, 0.01, 0.11, 0.09)

init_tau2 <- c(0.0003, 0.08, 0.9)
init_zeta2 <- c(0.10, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.1)

# - Set prior distributions
phi_0_prior_shape <- 0.0005
phi_scale <- 1
phi_1_prior_shape <- 0.05

zeta_prior1 <- c(15, 20, 3, 4, 23, 10, 1, 3, 2, 2, 10, 7)

# - Parameter storage

par_postJ1_gibbs_norm1 <- matrix(NA, 100000, 15)
par_postJ1_gibbs_norm2 <- matrix(NA, 100000, 15)

par_postJ1_metr_norm1 <- matrix(NA, 100000, 15)
par_postJ1_metr_norm2 <- matrix(NA, 100000, 15)

par_postJ1_gibbs_norm1[1,c(1:3)] <- init_tau1
par_postJ1_gibbs_norm1[1,c(4:15)] <- init_zeta1

par_postJ1_gibbs_norm2[1,c(1:3)] <- init_tau2
par_postJ1_gibbs_norm2[1,c(4:15)] <- init_zeta2

par_postJ1_metr_norm1[1,c(1:3)] <- init_tau1
par_postJ1_metr_norm1[1,c(4:15)] <- init_zeta1

par_postJ1_metr_norm2[1,c(1:3)] <- init_tau2
par_postJ1_metr_norm2[1,c(4:15)] <- init_zeta2

colnames(par_postJ1_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")
colnames(par_postJ1_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")
colnames(par_postJ1_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")
colnames(par_postJ1_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11")

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
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  
  par_postJ1_gibbs_norm1[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ g11)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ g21))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), scale=1/(phi_scale + H))
  
  par_postJ1_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum(g11 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0) + sum(g21 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0)
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * g11) + sum(d2 * g21), scale=1/(phi_scale + H1), b=1)
  
  par_postJ1_gibbs_norm1[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_gibbs_norm1[i, 2] <- beta
  
}

for(i in 2:100000){
  
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
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  
  par_postJ1_gibbs_norm2[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ g11)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ g21))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), scale=1/(phi_scale + H))
  
  par_postJ1_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum(g11 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0) + sum(g21 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0)
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * g11) + sum(d2 * g21), scale=1/(phi_scale + H1), b=1)
  
  par_postJ1_gibbs_norm2[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_gibbs_norm2[i, 2] <- beta
  
}

##### - Full metropolis

for(i in 2:100000){
  
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
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  
  par_postJ1_metr_norm1[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0038)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.0038) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.0038))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ1_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.055)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_star) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_star) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_1   ) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.055) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.055))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ1_metr_norm1[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_metr_norm1[i, 2] <- beta
  
}

for(i in 2:100000){
  
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
  gamma_zeta[1] <- rgamma(1, zeta_prior1[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior1[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior1[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior1[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior1[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior1[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  
  gamma_zeta[ 7] <- rgamma(1, zeta_prior1[ 7] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[ 8] <- rgamma(1, zeta_prior1[ 8] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[ 9] <- rgamma(1, zeta_prior1[ 9] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior1[10] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior1[11] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior1[12] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  
  par_postJ1_metr_norm2[i, 4:15] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0038)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + g11 * log(phi1_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + g21 * log(phi1_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.0038) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.0038))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ1_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.055)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_star) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ g11) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_star) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ g21) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1_1   ) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.055) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.055))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ1_metr_norm2[i, 3] <- phi1
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ g11) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ g11) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + g11 * log(phi1) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ g21) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ g21) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + g21 * log(phi1) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ1_metr_norm2[i, 2] <- beta
  
}

################################################# - J=2 - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.05, 0.09, 0.14, 0.05)
init_zeta1 <- c(12, 10, 2, 4, 20, 6, 2, 10, 4, 1, 2, 1, 1, 10, 1, 2, 4, 8)/100

init_tau2 <- c(0.0005, 0.06, 0.5, 0.15)
init_zeta2 <- c(0.10, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1)

# - Set prior distributions
phi_0_prior_shape <- 0.09
phi_1_prior_shape <- 0.5
phi_2_prior_shape <- 0.9

phi_scale <- 1

zeta_prior2 <- c(12, 10, 2, 4, 20, 6, 2, 10, 4, 1, 2, 1, 1, 10, 1, 2, 4, 8)
#zeta_prior2 <- c(9, 8, 7, 5, 6, 4, 1, 2, 3, 3, 2, 1, 4, 6, 5, 7, 8, 9)

# - Parameter storage

par_postJ2_gibbs_norm1 <- matrix(NA, 100000, 22)
par_postJ2_gibbs_norm2 <- matrix(NA, 100000, 22)

par_postJ2_metr_norm1 <- matrix(NA, 100000, 22)
par_postJ2_metr_norm2 <- matrix(NA, 100000, 22)

par_postJ2_gibbs_norm1[1,c(1:4)] <- init_tau1
par_postJ2_gibbs_norm1[1,c(5:22)] <- init_zeta1

par_postJ2_gibbs_norm2[1,c(1:4)] <- init_tau2
par_postJ2_gibbs_norm2[1,c(5:22)] <- init_zeta2

par_postJ2_metr_norm1[1,c(1:4)] <- init_tau1
par_postJ2_metr_norm1[1,c(5:22)] <- init_zeta1

par_postJ2_metr_norm2[1,c(1:4)] <- init_tau2
par_postJ2_metr_norm2[1,c(5:22)] <- init_zeta2

colnames(par_postJ2_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")
colnames(par_postJ2_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")

colnames(par_postJ2_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")
colnames(par_postJ2_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17")


##### - MCMC with Gibbs steps

for(i in 2:100000){
  
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
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  
  par_postJ2_gibbs_norm1[i, 5:22] <- gamma_zeta / sum(gamma_zeta)

  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), phi_scale + H)
  
  par_postJ2_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5) ) * phi0 * (phi2_1 ^ g12)) + sum( (g21 + g22) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5) ) * phi0 * (phi2_1 ^ g22))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * (g11 + g12)) + sum(d2 * (g21 + g22)), scale=1/(phi_scale + H1), b=1)
  
  par_postJ2_gibbs_norm1[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( g12 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 ) + sum( g22 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 )
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + sum(d1 * g12) + sum(d2 * g22), scale=1/(phi_scale + H2), b=1)
  
  par_postJ2_gibbs_norm1[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_gibbs_norm1[i, 2] <- beta
  
}

for(i in 2:100000){
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
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  
  par_postJ2_gibbs_norm2[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), phi_scale + H)
  
  par_postJ2_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5) ) * phi0 * (phi2_1 ^ g12)) + sum( (g21 + g22) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5) ) * phi0 * (phi2_1 ^ g22))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * (g11 + g12)) + sum(d2 * (g21 + g22)), scale=1/(phi_scale + H1), b=1)
  
  par_postJ2_gibbs_norm2[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( g12 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 ) + sum( g22 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 )
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + sum(d1 * g12) + sum(d2 * g22), scale=1/(phi_scale + H2), b=1)
  
  par_postJ2_gibbs_norm2[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_gibbs_norm2[i, 2] <- beta
  
}


for(i in 2:100000){
  
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
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  
  par_postJ2_metr_norm1[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.004)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.004) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.004))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ2_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_star) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_1   ) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_star) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_1   ) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.049))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ2_metr_norm1[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_star ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_star) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_1    ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_star ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_star) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_1    ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_1   ) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.112))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi1_unif > r, phi2_1, phi2_star)
  
  par_postJ2_metr_norm1[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_metr_norm1[i, 2] <- beta
  
}

for(i in 2:100000){
  
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
  gamma_zeta[1] <- rgamma(1, zeta_prior2[1] + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2] <- rgamma(1, zeta_prior2[2] + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3] <- rgamma(1, zeta_prior2[3] + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4] <- rgamma(1, zeta_prior2[4] + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5] <- rgamma(1, zeta_prior2[5] + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6] <- rgamma(1, zeta_prior2[6] + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7] <- rgamma(1, zeta_prior2[7] + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8] <- rgamma(1, zeta_prior2[8] + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9] <- rgamma(1, zeta_prior2[9] + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  
  gamma_zeta[10] <- rgamma(1, zeta_prior2[10] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior2[11] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior2[12] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior2[13] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior2[14] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior2[15] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior2[16] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior2[17] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior2[18] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  
  par_postJ2_metr_norm2[i, 5:22] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.004)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12) * log(phi1_1) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22) * log(phi1_1) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.004) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.004))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ2_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_star) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12)) * (phi2_1 ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1_1   ) + g12 * log(phi2_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_star) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22)) * (phi2_1 ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1_1   ) + g22 * log(phi2_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.049) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.049))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ2_metr_norm2[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_star ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_star) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12)) * (phi2_1    ^ g12) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_star ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_star) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22)) * (phi2_1    ^ g22) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2_1   ) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.112) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.112))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi1_unif > r, phi2_1, phi2_star)
  
  par_postJ2_metr_norm2[i, 4] <- phi2
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12)) * (phi2 ^ g12) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12) * log(phi1) + g12 * log(phi2) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22)) * (phi2 ^ g22) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22) * log(phi1) + g22 * log(phi2) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ2_metr_norm2[i, 2] <- beta
  
}

######################################################## - J=3  - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00004, 0.1, 0.02, 0.03, 0.02)
init_zeta1 <- c(12, 10, 2, 4, 10, 4, 1, 10, 4, 1, 10, 2, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 3, 7)/100

init_tau2 <- c(0.04, 0.06, 0.8, 0.5, 0.6)
init_zeta2 <- c(0.06, rep(0.04, 22), 0.06)

# - Set prior distributions
phi_0_prior_shape <- 0.09
phi_1_prior_shape <- 0.5
phi_2_prior_shape <- 0.7
phi_3_prior_shape <- 0.9

phi_scale <- 1

zeta_prior3 <- c(12, 10, 2, 4, 10, 4, 1, 10, 4, 1, 10, 2, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 3, 7)

# - Parameter storage

par_postJ3_gibbs_norm1 <- matrix(NA, 100000, 29)
par_postJ3_gibbs_norm2 <- matrix(NA, 100000, 29)

par_postJ3_metr_norm1 <- matrix(NA, 100000, 29)
par_postJ3_metr_norm2 <- matrix(NA, 100000, 29)

par_postJ3_gibbs_norm1[1,c(1:5)] <- init_tau1
par_postJ3_gibbs_norm1[1,c(6:29)] <- init_zeta1

par_postJ3_gibbs_norm2[1,c(1:5)] <- init_tau2
par_postJ3_gibbs_norm2[1,c(6:29)] <- init_zeta2

par_postJ3_metr_norm1[1,c(1:5)] <- init_tau1
par_postJ3_metr_norm1[1,c(6:29)] <- init_zeta1

par_postJ3_metr_norm2[1,c(1:5)] <- init_tau2
par_postJ3_metr_norm2[1,c(6:29)] <- init_zeta2

colnames(par_postJ3_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")
colnames(par_postJ3_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")

colnames(par_postJ3_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")
colnames(par_postJ3_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23")


##### - MCMC with Gibbs steps

for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  
  par_postJ3_gibbs_norm1[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), scale=1/(phi_scale + H))
  
  par_postJ3_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum( (g21 + g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * (g11 + g12 + g13)) + sum(d2 * (g21 + g22 + g23)), scale=1/(phi_scale + H1), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ3_gibbs_norm1[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g13)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g23))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + sum(d1 * (g12 + g13)) + sum(d2 * (g22 + g23)), scale=1/(phi_scale + H2), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_2_prior_shape + sum(d * (g2 + g3)), rate = phi_scale + H2)#rgamma(1, phi_2_prior_shape + sum(d * (g2 + g3)), phi_scale + H2)
  
  par_postJ3_gibbs_norm1[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( g13 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2) + sum( g23 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2)

  phi3 <- rtgamma(1, shape=phi_3_prior_shape + sum(d1 * g13) + sum(d2 * g23), scale=1/(phi_scale + H3), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ3_gibbs_norm1[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_gibbs_norm1[i, 2] <- beta
  
}

for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  
  par_postJ3_gibbs_norm2[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), scale=1/(phi_scale + H))
  
  par_postJ3_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13)) + sum( (g21 + g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * (g11 + g12 + g13)) + sum(d2 * (g21 + g22 + g23)), scale=1/(phi_scale + H1), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ3_gibbs_norm2[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g13)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ g23))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + sum(d1 * (g12 + g13)) + sum(d2 * (g22 + g23)), scale=1/(phi_scale + H2), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_2_prior_shape + sum(d * (g2 + g3)), rate = phi_scale + H2)#rgamma(1, phi_2_prior_shape + sum(d * (g2 + g3)), phi_scale + H2)
  
  par_postJ3_gibbs_norm2[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( g13 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2) + sum( g23 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2)
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + sum(d1 * g13) + sum(d2 * g23), scale=1/(phi_scale + H3), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ3_gibbs_norm2[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_gibbs_norm2[i, 2] <- beta
  
}


for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  
  par_postJ3_metr_norm1[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0034)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi0_star, a=0, mean = phi0_1, sd = 0.0034) / dtruncnorm(phi0_1, a=0, mean = phi0_star, sd = 0.0034))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)  
  
  par_postJ3_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.058)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_star) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_1   ) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_star) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_1   ) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi1_star, a=0, mean = phi1_1, sd = 0.058) / dtruncnorm(phi1_1, a=0, mean = phi1_star, sd = 0.058))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)  
  
  par_postJ3_metr_norm1[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.115)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_star ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_star) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_1    ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_1   ) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_star ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_star) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_1    ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_1   ) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi2_star, a=0, mean = phi2_1, sd = 0.115) / dtruncnorm(phi2_1, a=0, mean = phi2_star, sd = 0.115))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)  
  
  par_postJ3_metr_norm1[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.125)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_star ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_star) + beta_1 * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_1    ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_star ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_star) + beta_1 * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_1    ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_1   ) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi3_star, a=0, mean = phi3_1, sd = 0.125) / dtruncnorm(phi3_1, a=0, mean = phi3_star, sd = 0.125))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)  
  
  par_postJ3_metr_norm1[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_metr_norm1[i, 2] <- beta
  
}

for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior3[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior3[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior3[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior3[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior3[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior3[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior3[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior3[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior3[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior3[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior3[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior3[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  
  gamma_zeta[13] <- rgamma(1, zeta_prior3[13] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior3[14] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior3[15] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[16] <- rgamma(1, zeta_prior3[16] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior3[17] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior3[18] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior3[19] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior3[20] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior3[21] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior3[22] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior3[23] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior3[24] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  
  par_postJ3_metr_norm2[i, 6:29] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.0034)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13) * log(phi1_1) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23) * log(phi1_1) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi0_star, a=0, mean = phi0_1, sd = 0.0034) / dtruncnorm(phi0_1, a=0, mean = phi0_star, sd = 0.0034))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)  
  
  par_postJ3_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.058)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_star) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13)) * (phi2_1 ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1_1   ) + (g12 + g13) * log(phi2_1) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_star) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23)) * (phi2_1 ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1_1   ) + (g22 + g23) * log(phi2_1) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi1_star, a=0, mean = phi1_1, sd = 0.058) / dtruncnorm(phi1_1, a=0, mean = phi1_star, sd = 0.058))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)  
  
  par_postJ3_metr_norm2[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.115)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_star ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_star) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2_1    ^ (g12 + g13)) * (phi3_1 ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2_1   ) + g13 * log(phi3_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_star ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_star) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2_1    ^ (g22 + g23)) * (phi3_1 ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2_1   ) + g23 * log(phi3_1) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi2_star, a=0, mean = phi2_1, sd = 0.115) / dtruncnorm(phi2_1, a=0, mean = phi2_star, sd = 0.115))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)  
  
  par_postJ3_metr_norm2[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.125)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_star ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_star) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3_1    ^ g13) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_star ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_star) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3_1    ^ g23) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3_1   ) + beta_1 * (x2 + t2 - 77.5)))  ) / (dtruncnorm(phi3_star, a=0, mean = phi3_1, sd = 0.125) / dtruncnorm(phi3_1, a=0, mean = phi3_star, sd = 0.125))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)  
  
  par_postJ3_metr_norm2[i, 5] <- phi3
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13)) * (phi2 ^ (g12 + g13)) * (phi3 ^ g13) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13) * log(phi1) + (g12 + g13) * log(phi2) + g13 * log(phi3) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23)) * (phi2 ^ (g22 + g23)) * (phi3 ^ g23) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23) * log(phi1) + (g22 + g23) * log(phi2) + g23 * log(phi3) + beta_1    * (x2 + t2 - 77.5)))  ) / (dtruncnorm(beta_star, a=0, mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=0, mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ3_metr_norm2[i, 2] <- beta
  
}

######################################################## - J=4  - ###################################################################

# - Initialization MCMC - Step 0

init_tau1 <- c(0.00004, 0.1, 0.02, 0.03, 0.02, 0.9)
init_zeta1 <- c(10, 5, 1, 4, 10, 2, 4, 10, 4, 0.5, 10, 3, 0.5, 5, 2, 1, 2, 0.5, 1, 2, 0.5, 1, 4, 1.5, .5, 4, 2.5, 0.5, 4, 5)/100

init_tau2 <- c(0.04, 0.06, 0.8, 0.5, 0.6, 0.04)
init_zeta2 <- c(0.08, rep(0.04, 28), 0.08)

# - Set prior distributions
phi_0_prior_shape <- 0.09
phi_1_prior_shape <- 0.5
phi_2_prior_shape <- 0.75
phi_3_prior_shape <- 0.85
phi_4_prior_shape <- 0.9

phi_scale <- 1

zeta_prior4 <- c(10, 5, 1, 4, 10, 2, 4, 10, 4, 1, 10, 3, 1, 5, 2, 1, 2, 1, 1, 2, 1, 1, 4, 2, 1, 4, 3, 1, 4, 5)

# - Parameter storage

par_postJ4_gibbs_norm1 <- matrix(NA, 100000, 36)
par_postJ4_gibbs_norm2 <- matrix(NA, 100000, 36)

par_postJ4_metr_norm1 <- matrix(NA, 100000, 36)
par_postJ4_metr_norm2 <- matrix(NA, 100000, 36)

par_postJ4_gibbs_norm1[1,c(1:6)] <- init_tau1
par_postJ4_gibbs_norm1[1,c(7:36)] <- init_zeta1

par_postJ4_gibbs_norm2[1,c(1:6)] <- init_tau2
par_postJ4_gibbs_norm2[1,c(7:36)] <- init_zeta2

par_postJ4_metr_norm1[1,c(1:6)] <- init_tau1
par_postJ4_metr_norm1[1,c(7:36)] <- init_zeta1

par_postJ4_metr_norm2[1,c(1:6)] <- init_tau2
par_postJ4_metr_norm2[1,c(7:36)] <- init_zeta2

colnames(par_postJ4_gibbs_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")
colnames(par_postJ4_gibbs_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")

colnames(par_postJ4_metr_norm1) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")
colnames(par_postJ4_metr_norm2) <- c("exp_alpha", "beta", "exp_gamma_1", "exp_gamma_2", "exp_gamma_3", "exp_gamma_4", "zeta_0", "zeta_1", "zeta_2", "zeta_3", "zeta_4", "zeta_5", "zeta_6", "zeta_7", "zeta_8", "zeta_9", "zeta_10", "zeta_11", "zeta_12", "zeta_13", "zeta_14", "zeta_15", "zeta_16", "zeta_17", "zeta_18", "zeta_19", "zeta_20", "zeta_21", "zeta_22", "zeta_23", "zeta_24", "zeta_25", "zeta_26", "zeta_27", "zeta_28", "zeta_29")


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
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + sum(bL * c10 * g14) + sum(b2L * c0 * g24), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + sum(bL * c11 * g14) + sum(b2L * c1 * g24), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + sum(bL * c12 * g14) + sum(b2L * c2 * g24), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + sum(bH * c10 * g14) + sum(b2H * c0 * g24), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + sum(bH * c11 * g14) + sum(b2H * c1 * g24), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + sum(bH * c12 * g14) + sum(b2H * c2 * g24), 1)
  
  par_postJ4_gibbs_norm1[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), scale=1/(phi_scale + H))
  
  par_postJ4_gibbs_norm1[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g21 + g22 + g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * (g11 + g12 + g13 + g14)) + sum(d2 * (g21 + g22 + g23 + g24)), scale=1/(phi_scale + H1), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ4_gibbs_norm1[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + sum(d1 * (g12 + g13 + g14)) + sum(d2 * (g22 + g23 + g24)), scale=1/(phi_scale + H2), b=1)
  
  par_postJ4_gibbs_norm1[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( (g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g14)) + sum( (g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g24))
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + sum(d1 * (g13 + g14)) + sum(d2 * (g23 + g24)), scale=1/(phi_scale + H3), b=1)
  
  par_postJ4_gibbs_norm1[i, 5] <- phi3
  
  # - phi_4
  
  H4 <- sum( g14 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * phi3) + sum( g24 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * phi3)
  
  phi4 <- rtgamma(1, shape=phi_4_prior_shape + sum(d1 * g14) + sum(d2 * g24), scale=1/(phi_scale + H4), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ4_gibbs_norm1[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
           sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
           sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
           sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_gibbs_norm1[i, 2] <- beta
  
  
}

for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + sum(bL * c10 * g14) + sum(b2L * c0 * g24), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + sum(bL * c11 * g14) + sum(b2L * c1 * g24), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + sum(bL * c12 * g14) + sum(b2L * c2 * g24), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + sum(bH * c10 * g14) + sum(b2H * c0 * g24), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + sum(bH * c11 * g14) + sum(b2H * c1 * g24), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + sum(bH * c12 * g14) + sum(b2H * c2 * g24), 1)
  
  par_postJ4_gibbs_norm2[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  H <- sum(((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum(((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi0 <- rgamma(1, phi_0_prior_shape + sum(d1) + sum(d2), scale=1/(phi_scale + H))
  
  par_postJ4_gibbs_norm2[i, 1] <- phi0
  
  # - phi_1
  
  H1 <- sum( (g11 + g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g21 + g22 + g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi1 <- rtgamma(1, shape=phi_1_prior_shape + sum(d1 * (g11 + g12 + g13 + g14)) + sum(d2 * (g21 + g22 + g23 + g24)), scale=1/(phi_scale + H1), b=1) #rgamma(1, phi_1_prior_shape + sum(d * (g1 + g2 + g3)), phi_scale + H1) #rtrunc(1, spec="gamma", b = 1, shape = phi_1_prior_shape + sum(d * (g1 + g2 + g3)), scale=1/(phi_scale + H1)) #
  
  par_postJ4_gibbs_norm2[i, 3] <- phi1
  
  # - phi_2
  
  H2 <- sum( (g12 + g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14)) + sum( (g22 + g23) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24))
  
  phi2 <- rtgamma(1, shape=phi_2_prior_shape + sum(d1 * (g12 + g13 + g14)) + sum(d2 * (g22 + g23 + g24)), scale=1/(phi_scale + H2), b=1)
  
  par_postJ4_gibbs_norm2[i, 4] <- phi2
  
  # - phi_3
  
  H3 <- sum( (g13 + g14) * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g14)) + sum( (g23 + g24) * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * (phi4_1 ^ g24))
  
  phi3 <- rtgamma(1, shape=phi_3_prior_shape + sum(d1 * (g13 + g14)) + sum(d2 * (g23 + g24)), scale=1/(phi_scale + H3), b=1)
  
  par_postJ4_gibbs_norm2[i, 5] <- phi3
  
  # - phi_4
  
  H4 <- sum( g14 * ((exp(beta_1 * t1) - 1) / beta_1) * exp(beta_1 * (x1 - 77.5)) * phi0 * phi1 * phi2 * phi3) + sum( g24 * ((exp(beta_1 * t2) - 1) / beta_1) * exp(beta_1 * (x2 - 77.5)) * phi0 * phi1 * phi2 * phi3)
  
  phi4 <- rtgamma(1, shape=phi_4_prior_shape + sum(d1 * g14) + sum(d2 * g24), scale=1/(phi_scale + H4), b=1)#rtrunc(1, spec="gamma", b=1, shape=phi_3_prior_shape + sum(d * g3), scale=phi_scale + H3)#rgamma(1, phi_3_prior_shape + sum(d * g3), phi_scale + H3)
  
  par_postJ4_gibbs_norm2[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  #beta_star <- rtruncnorm(1, a=0.001, mean=beta1[i-1], sd=0.1)
  #beta_star <- min(0.02, runif(1, beta1[i-1] - 0.05, beta1[i-1] + 0.05))
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_gibbs_norm2[i, 2] <- beta
  
  
}


for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + sum(bL * c10 * g14) + sum(b2L * c0 * g24), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + sum(bL * c11 * g14) + sum(b2L * c1 * g24), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + sum(bL * c12 * g14) + sum(b2L * c2 * g24), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + sum(bH * c10 * g14) + sum(b2H * c0 * g24), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + sum(bH * c11 * g14) + sum(b2H * c1 * g24), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + sum(bH * c12 * g14) + sum(b2H * c2 * g24), 1)
  
  par_postJ4_metr_norm1[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.00375)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.00375) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.00375))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ4_metr_norm1[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.06)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_star) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_1   ) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_star) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_1   ) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.06) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.06))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ4_metr_norm1[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.1)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_star ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_star) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_1    ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_1   ) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_star ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_star) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_1    ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_1   ) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.1) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.1))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)
  
  par_postJ4_metr_norm1[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.11)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_star ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_star) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_1    ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_1   ) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_star ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_star) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_1    ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_1   ) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi3_star, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.11) / dtruncnorm(phi3_1, a=(10^(-20)), b=1, mean = phi3_star, sd = 0.11))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)
  
  par_postJ4_metr_norm1[i, 5] <- phi3
  
  # - phi_4
  
  phi4_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.2)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_star ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_star) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_1    ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_star ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_star) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_1    ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_1   ) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi4_star, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.2) / dtruncnorm(phi4_1, a=(10^(-20)), b=1, mean = phi4_star, sd = 0.2))
  
  phi4_unif <- runif(1, 0, 1)
  
  phi4 <- ifelse( phi4_unif > r, phi4_1, phi4_star)
  
  par_postJ4_metr_norm1[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_metr_norm1[i, 2] <- beta

}

for(i in 2:100000){
  
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
  gamma_zeta[1]  <- rgamma(1, zeta_prior4[1]  + sum(bL * c10 * g10) + sum(b2L * c0 * g20), 1)
  gamma_zeta[2]  <- rgamma(1, zeta_prior4[2]  + sum(bL * c11 * g10) + sum(b2L * c1 * g20), 1)
  gamma_zeta[3]  <- rgamma(1, zeta_prior4[3]  + sum(bL * c12 * g10) + sum(b2L * c2 * g20), 1)
  gamma_zeta[4]  <- rgamma(1, zeta_prior4[4]  + sum(bL * c10 * g11) + sum(b2L * c0 * g21), 1)
  gamma_zeta[5]  <- rgamma(1, zeta_prior4[5]  + sum(bL * c11 * g11) + sum(b2L * c1 * g21), 1)
  gamma_zeta[6]  <- rgamma(1, zeta_prior4[6]  + sum(bL * c12 * g11) + sum(b2L * c2 * g21), 1)
  gamma_zeta[7]  <- rgamma(1, zeta_prior4[7]  + sum(bL * c10 * g12) + sum(b2L * c0 * g22), 1)
  gamma_zeta[8]  <- rgamma(1, zeta_prior4[8]  + sum(bL * c11 * g12) + sum(b2L * c1 * g22), 1)
  gamma_zeta[9]  <- rgamma(1, zeta_prior4[9]  + sum(bL * c12 * g12) + sum(b2L * c2 * g22), 1)
  gamma_zeta[10] <- rgamma(1, zeta_prior4[10] + sum(bL * c10 * g13) + sum(b2L * c0 * g23), 1)
  gamma_zeta[11] <- rgamma(1, zeta_prior4[11] + sum(bL * c11 * g13) + sum(b2L * c1 * g23), 1)
  gamma_zeta[12] <- rgamma(1, zeta_prior4[12] + sum(bL * c12 * g13) + sum(b2L * c2 * g23), 1)
  gamma_zeta[13] <- rgamma(1, zeta_prior4[13] + sum(bL * c10 * g14) + sum(b2L * c0 * g24), 1)
  gamma_zeta[14] <- rgamma(1, zeta_prior4[14] + sum(bL * c11 * g14) + sum(b2L * c1 * g24), 1)
  gamma_zeta[15] <- rgamma(1, zeta_prior4[15] + sum(bL * c12 * g14) + sum(b2L * c2 * g24), 1)
  
  gamma_zeta[16] <- rgamma(1, zeta_prior4[16] + sum(bH * c10 * g10) + sum(b2H * c0 * g20), 1)
  gamma_zeta[17] <- rgamma(1, zeta_prior4[17] + sum(bH * c11 * g10) + sum(b2H * c1 * g20), 1)
  gamma_zeta[18] <- rgamma(1, zeta_prior4[18] + sum(bH * c12 * g10) + sum(b2H * c2 * g20), 1)
  gamma_zeta[19] <- rgamma(1, zeta_prior4[19] + sum(bH * c10 * g11) + sum(b2H * c0 * g21), 1)
  gamma_zeta[20] <- rgamma(1, zeta_prior4[20] + sum(bH * c11 * g11) + sum(b2H * c1 * g21), 1)
  gamma_zeta[21] <- rgamma(1, zeta_prior4[21] + sum(bH * c12 * g11) + sum(b2H * c2 * g21), 1)
  gamma_zeta[22] <- rgamma(1, zeta_prior4[22] + sum(bH * c10 * g12) + sum(b2H * c0 * g22), 1)
  gamma_zeta[23] <- rgamma(1, zeta_prior4[23] + sum(bH * c11 * g12) + sum(b2H * c1 * g22), 1)
  gamma_zeta[24] <- rgamma(1, zeta_prior4[24] + sum(bH * c12 * g12) + sum(b2H * c2 * g22), 1)
  gamma_zeta[25] <- rgamma(1, zeta_prior4[25] + sum(bH * c10 * g13) + sum(b2H * c0 * g23), 1)
  gamma_zeta[26] <- rgamma(1, zeta_prior4[26] + sum(bH * c11 * g13) + sum(b2H * c1 * g23), 1)
  gamma_zeta[27] <- rgamma(1, zeta_prior4[27] + sum(bH * c12 * g13) + sum(b2H * c2 * g23), 1)
  gamma_zeta[28] <- rgamma(1, zeta_prior4[28] + sum(bH * c10 * g14) + sum(b2H * c0 * g24), 1)
  gamma_zeta[29] <- rgamma(1, zeta_prior4[29] + sum(bH * c11 * g14) + sum(b2H * c1 * g24), 1)
  gamma_zeta[30] <- rgamma(1, zeta_prior4[30] + sum(bH * c12 * g14) + sum(b2H * c2 * g24), 1)
  
  par_postJ4_metr_norm2[i, 7:36] <- gamma_zeta / sum(gamma_zeta)
  
  # - Step 3 - Generate tau
  
  # - phi0
  
  phi0_star <- rtruncnorm(1, a=(10^(-20)), mean = phi0_1, sd = 0.00375)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_star) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0_1   ) + (g11 + g12 + g13 + g14) * log(phi1_1) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_star * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_star) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0_1    * (phi1_1 ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0_1   ) + (g21 + g22 + g23 + g24) * log(phi1_1) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi0_star, a=(10^(-20)), mean = phi0_1, sd = 0.00375) / dtruncnorm(phi0_1, a=(10^(-20)), mean = phi0_star, sd = 0.00375))
  
  phi0_unif <- runif(1, 0, 1)
  
  phi0 <- ifelse( phi0_unif > r, phi0_1, phi0_star)
  
  par_postJ4_metr_norm2[i, 1] <- phi0
  
  # - phi_1
  
  phi1_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.06)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_star ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_star) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1_1    ^ (g11 + g12 + g13 + g14)) * (phi2_1 ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1_1   ) + (g12 + g13 + g14) * log(phi2_1) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_star ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_star) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1_1    ^ (g21 + g22 + g23 + g24)) * (phi2_1 ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1_1   ) + (g22 + g23 + g24) * log(phi2_1) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi1_star, a=(10^(-20)), b=1, mean = phi1_1, sd = 0.06) / dtruncnorm(phi1_1, a=(10^(-20)), b=1, mean = phi1_star, sd = 0.06))
  
  phi1_unif <- runif(1, 0, 1)
  
  phi1 <- ifelse( phi1_unif > r, phi1_1, phi1_star)
  
  par_postJ4_metr_norm2[i, 3] <- phi1
  
  # - phi_2
  
  phi2_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.1)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_star ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_star) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2_1    ^ (g12 + g13 + g14)) * (phi3_1 ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2_1   ) + (g13 + g14) * log(phi3_1) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_star ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_star) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2_1    ^ (g22 + g23 + g24)) * (phi3_1 ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2_1   ) + (g23 + g24) * log(phi3_1) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi2_star, a=(10^(-20)), b=1, mean = phi2_1, sd = 0.1) / dtruncnorm(phi2_1, a=(10^(-20)), b=1, mean = phi2_star, sd = 0.1))
  
  phi2_unif <- runif(1, 0, 1)
  
  phi2 <- ifelse( phi2_unif > r, phi2_1, phi2_star)
  
  par_postJ4_metr_norm2[i, 4] <- phi2
  
  # - phi_3
  
  phi3_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.11)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_star ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_star) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3_1    ^ (g13 + g14)) * (phi4_1 ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3_1   ) + g14 * log(phi4_1) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_star ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_star) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3_1    ^ (g23 + g24)) * (phi4_1 ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3_1   ) + g24 * log(phi4_1) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi3_star, a=(10^(-20)), b=1, mean = phi3_1, sd = 0.11) / dtruncnorm(phi3_1, a=(10^(-20)), b=1, mean = phi3_star, sd = 0.11))
  
  phi3_unif <- runif(1, 0, 1)
  
  phi3 <- ifelse( phi3_unif > r, phi3_1, phi3_star)
  
  par_postJ4_metr_norm2[i, 5] <- phi3
  
  # - phi_4
  
  phi4_star <- rtruncnorm(1, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.2)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_star ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_star) + beta_1 * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1 * t1) - 1) / beta_1) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4_1    ^ g14) * exp(beta_1 * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4_1   ) + beta_1 * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_star ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_star) + beta_1 * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1 * t2) - 1) / beta_1) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4_1    ^ g24) * exp(beta_1 * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4_1   ) + beta_1 * (x2 + t2 - 77.5)))) / (dtruncnorm(phi4_star, a=(10^(-20)), b=1, mean = phi4_1, sd = 0.2) / dtruncnorm(phi4_1, a=(10^(-20)), b=1, mean = phi4_star, sd = 0.2))
  
  phi4_unif <- runif(1, 0, 1)
  
  phi4 <- ifelse( phi4_unif > r, phi4_1, phi4_star)
  
  par_postJ4_metr_norm2[i, 6] <- phi4
  
  # - beta
  
  # - calculate beta star from the proposal distribution
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta_1, sd = 0.008)  # runif(1, 0.01, 0.2)
  
  # - calculate ratio r from uniform
  
  r <- exp(sum(- ((exp(beta_star * t1) - 1) / beta_star) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_star * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_star * (x1 + t1 - 77.5))) - 
             sum(- ((exp(beta_1    * t1) - 1) / beta_1   ) * phi0 * (phi1 ^ (g11 + g12 + g13 + g14)) * (phi2 ^ (g12 + g13 + g14)) * (phi3 ^ (g13 + g14)) * (phi4 ^ g14) * exp(beta_1    * (x1 - 77.5)) + d1 * (log(phi0) + (g11 + g12 + g13 + g14) * log(phi1) + (g12 + g13 + g14) * log(phi2) + (g13 + g14) * log(phi3) + g14 * log(phi4) + beta_1    * (x1 + t1 - 77.5))) + 
             sum(- ((exp(beta_star * t2) - 1) / beta_star) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_star * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_star * (x2 + t2 - 77.5))) - 
             sum(- ((exp(beta_1    * t2) - 1) / beta_1   ) * phi0 * (phi1 ^ (g21 + g22 + g23 + g24)) * (phi2 ^ (g22 + g23 + g24)) * (phi3 ^ (g23 + g24)) * (phi4 ^ g24) * exp(beta_1    * (x2 - 77.5)) + d2 * (log(phi0) + (g21 + g22 + g23 + g24) * log(phi1) + (g22 + g23 + g24) * log(phi2) + (g23 + g24) * log(phi3) + g24 * log(phi4) + beta_1    * (x2 + t2 - 77.5)))) / (dtruncnorm(beta_star, a=(10^(-5)), mean = beta_1, sd = 0.008) / dtruncnorm(beta_1, a=(10^(-5)), mean = beta_star, sd = 0.008))
  
  beta_unif <- runif(1, 0, 1)
  
  beta <- ifelse( beta_unif > r, beta_1, beta_star)
  
  par_postJ4_metr_norm2[i, 2] <- beta
  
}

############################################### - Posterior summaries - #################################################

par_postJ0_gn1_th <- as.data.frame(par_postJ0_gn1_th)
par_postJ0_gn2_th <- as.data.frame(par_postJ0_gn2_th)
par_postJ0_mn1_th <- as.data.frame(par_postJ0_mn1_th)
par_postJ0_mn2_th <- as.data.frame(par_postJ0_mn2_th)

par_postJ1_gn1_th <- as.data.frame(par_postJ1_gn1_th)
par_postJ1_gn2_th <- as.data.frame(par_postJ1_gn2_th)
par_postJ1_mn1_th <- as.data.frame(par_postJ1_mn1_th)
par_postJ1_mn2_th <- as.data.frame(par_postJ1_mn2_th)

par_postJ2_gn1_th <- as.data.frame(par_postJ2_gn1_th)
par_postJ2_gn2_th <- as.data.frame(par_postJ2_gn2_th)
par_postJ2_mn1_th <- as.data.frame(par_postJ2_mn1_th)
par_postJ2_mn2_th <- as.data.frame(par_postJ2_mn2_th)

par_postJ3_gn1_th <- as.data.frame(par_postJ3_gn1_th)
par_postJ3_gn2_th <- as.data.frame(par_postJ3_gn2_th)
par_postJ3_mn1_th <- as.data.frame(par_postJ3_mn1_th)
par_postJ3_mn2_th <- as.data.frame(par_postJ3_mn2_th)

par_postJ4_gn1_th <- as.data.frame(par_postJ4_gn1_th)
par_postJ4_gn2_th <- as.data.frame(par_postJ4_gn2_th)
par_postJ4_mn1_th <- as.data.frame(par_postJ4_mn1_th)
par_postJ4_mn2_th <- as.data.frame(par_postJ4_mn2_th)

#### - Mean

par_post_meanJ0_C11 <- colMeans(par_postJ0_gn1_th)
par_post_meanJ0_C12 <- colMeans(par_postJ0_gn2_th)
par_post_meanJ0_C21 <- colMeans(par_postJ0_mn1_th)
par_post_meanJ0_C22 <- colMeans(par_postJ0_mn2_th)

par_post_meanJ1_C11 <- colMeans(par_postJ1_gn1_th)
par_post_meanJ1_C12 <- colMeans(par_postJ1_gn2_th)
par_post_meanJ1_C21 <- colMeans(par_postJ1_mn1_th)
par_post_meanJ1_C22 <- colMeans(par_postJ1_mn2_th)

par_post_meanJ2_C11 <- colMeans(par_postJ2_gn1_th)
par_post_meanJ2_C12 <- colMeans(par_postJ2_gn2_th)
par_post_meanJ2_C21 <- colMeans(par_postJ2_mn1_th)
par_post_meanJ2_C22 <- colMeans(par_postJ2_mn2_th)

par_post_meanJ3_C11 <- colMeans(par_postJ3_gn1_th)
par_post_meanJ3_C12 <- colMeans(par_postJ3_gn2_th)
par_post_meanJ3_C21 <- colMeans(par_postJ3_mn1_th)
par_post_meanJ3_C22 <- colMeans(par_postJ3_mn2_th)

par_post_meanJ4_C11 <- colMeans(par_postJ4_gn1_th)
par_post_meanJ4_C12 <- colMeans(par_postJ4_gn2_th)
par_post_meanJ4_C21 <- colMeans(par_postJ4_mn1_th)
par_post_meanJ4_C22 <- colMeans(par_postJ4_mn2_th)

#### - Median

par_post_medianJ0_C11 <- apply(par_postJ0_gn1_th, MARGIN=2, FUN="median")
par_post_medianJ0_C12 <- apply(par_postJ0_gn2_th, MARGIN=2, FUN="median")
par_post_medianJ0_C21 <- apply(par_postJ0_mn1_th, MARGIN=2, FUN="median")
par_post_medianJ0_C22 <- apply(par_postJ0_mn2_th, MARGIN=2, FUN="median")

par_post_medianJ1_C11 <- apply(par_postJ1_gn1_th, MARGIN=2, FUN="median")
par_post_medianJ1_C12 <- apply(par_postJ1_gn2_th, MARGIN=2, FUN="median")
par_post_medianJ1_C21 <- apply(par_postJ1_mn1_th, MARGIN=2, FUN="median")
par_post_medianJ1_C22 <- apply(par_postJ1_mn2_th, MARGIN=2, FUN="median")

par_post_medianJ2_C11 <- apply(par_postJ2_gn1_th, MARGIN=2, FUN="median")
par_post_medianJ2_C12 <- apply(par_postJ2_gn2_th, MARGIN=2, FUN="median")
par_post_medianJ2_C21 <- apply(par_postJ2_mn1_th, MARGIN=2, FUN="median")
par_post_medianJ2_C22 <- apply(par_postJ2_mn2_th, MARGIN=2, FUN="median")

par_post_medianJ3_C11 <- apply(par_postJ3_gn1_th, MARGIN=2, FUN="median")
par_post_medianJ3_C12 <- apply(par_postJ3_gn2_th, MARGIN=2, FUN="median")
par_post_medianJ3_C21 <- apply(par_postJ3_mn1_th, MARGIN=2, FUN="median")
par_post_medianJ3_C22 <- apply(par_postJ3_mn2_th, MARGIN=2, FUN="median")

par_post_medianJ4_C11 <- apply(par_postJ4_gn1_th, MARGIN=2, FUN="median")
par_post_medianJ4_C12 <- apply(par_postJ4_gn2_th, MARGIN=2, FUN="median")
par_post_medianJ4_C21 <- apply(par_postJ4_mn1_th, MARGIN=2, FUN="median")
par_post_medianJ4_C22 <- apply(par_postJ4_mn2_th, MARGIN=2, FUN="median")


#### - Standard deviation

# - J=0
par_post_sdJ0_C11 <- rep(0, 2)
par_post_sdJ0_C12 <- rep(0, 2)
par_post_sdJ0_C21 <- rep(0, 2)
par_post_sdJ0_C22 <- rep(0, 2)

for(i in c(1,2)){
  par_post_sdJ0_C11[i] <- sd(par_postJ0_gn1_th[,i])
  par_post_sdJ0_C12[i] <- sd(par_postJ0_gn2_th[,i])
  par_post_sdJ0_C21[i] <- sd(par_postJ0_mn1_th[,i])
  par_post_sdJ0_C22[i] <- sd(par_postJ0_mn2_th[,i])
}

# - J=1
par_post_sdJ1_C11 <- rep(0, 15)
par_post_sdJ1_C12 <- rep(0, 15)
par_post_sdJ1_C21 <- rep(0, 15)
par_post_sdJ1_C22 <- rep(0, 15)

for(i in c(1:15)){
  par_post_sdJ1_C11[i] <- sd(par_postJ1_gn1_th[,i])
  par_post_sdJ1_C12[i] <- sd(par_postJ1_gn2_th[,i])
  par_post_sdJ1_C21[i] <- sd(par_postJ1_mn1_th[,i])
  par_post_sdJ1_C22[i] <- sd(par_postJ1_mn2_th[,i])
}

# - J=2
par_post_sdJ2_C11 <- rep(0, 22)
par_post_sdJ2_C12 <- rep(0, 22)
par_post_sdJ2_C21 <- rep(0, 22)
par_post_sdJ2_C22 <- rep(0, 22)

for(i in c(1:22)){
  par_post_sdJ2_C11[i] <- sd(par_postJ2_gn1_th[,i])
  par_post_sdJ2_C12[i] <- sd(par_postJ2_gn2_th[,i])
  par_post_sdJ2_C21[i] <- sd(par_postJ2_mn1_th[,i])
  par_post_sdJ2_C22[i] <- sd(par_postJ2_mn2_th[,i])
}

# - J=3
par_post_sdJ3_C11 <- rep(0, 29)
par_post_sdJ3_C12 <- rep(0, 29)
par_post_sdJ3_C21 <- rep(0, 29)
par_post_sdJ3_C22 <- rep(0, 29)

for(i in c(1:29)){
  par_post_sdJ3_C11[i] <- sd(par_postJ3_gn1_th[,i])
  par_post_sdJ3_C12[i] <- sd(par_postJ3_gn2_th[,i])
  par_post_sdJ3_C21[i] <- sd(par_postJ3_mn1_th[,i])
  par_post_sdJ3_C22[i] <- sd(par_postJ3_mn2_th[,i])
}

# - J=4
par_post_sdJ4_C11 <- rep(0, 36)
par_post_sdJ4_C12 <- rep(0, 36)
par_post_sdJ4_C21 <- rep(0, 36)
par_post_sdJ4_C22 <- rep(0, 36)

for(i in c(1:36)){
  par_post_sdJ4_C11[i] <- sd(par_postJ4_gn1_th[,i])
  par_post_sdJ4_C12[i] <- sd(par_postJ4_gn2_th[,i])
  par_post_sdJ4_C21[i] <- sd(par_postJ4_mn1_th[,i])
  par_post_sdJ4_C22[i] <- sd(par_postJ4_mn2_th[,i])
}

#### - Quantiles

# - J=0
par_post_lqJ0_C11 <- rep(0, 2)
par_post_lqJ0_C12 <- rep(0, 2)
par_post_lqJ0_C21 <- rep(0, 2)
par_post_lqJ0_C22 <- rep(0, 2)

for(i in c(1,2)){
  par_post_lqJ0_C11[i] <- quantile(par_postJ0_gn1_th[,i], prob=0.025)
  par_post_lqJ0_C12[i] <- quantile(par_postJ0_gn2_th[,i], prob=0.025)
  par_post_lqJ0_C21[i] <- quantile(par_postJ0_mn1_th[,i], prob=0.025)
  par_post_lqJ0_C22[i] <- quantile(par_postJ0_mn2_th[,i], prob=0.025)
}

par_post_uqJ0_C11 <- rep(0, 2)
par_post_uqJ0_C12 <- rep(0, 2)
par_post_uqJ0_C21 <- rep(0, 2)
par_post_uqJ0_C22 <- rep(0, 2)

for(i in c(1,2)){
  par_post_uqJ0_C11[i] <- quantile(par_postJ0_gn1_th[,i], prob=0.975)
  par_post_uqJ0_C12[i] <- quantile(par_postJ0_gn2_th[,i], prob=0.975)
  par_post_uqJ0_C21[i] <- quantile(par_postJ0_mn1_th[,i], prob=0.975)
  par_post_uqJ0_C22[i] <- quantile(par_postJ0_mn2_th[,i], prob=0.975)
}

# - J=1
par_post_lqJ1_C11 <- rep(0, 15)
par_post_lqJ1_C12 <- rep(0, 15)
par_post_lqJ1_C21 <- rep(0, 15)
par_post_lqJ1_C22 <- rep(0, 15)

for(i in c(1:15)){
  par_post_lqJ1_C11[i] <- quantile(par_postJ1_gn1_th[,i], prob=0.025)
  par_post_lqJ1_C12[i] <- quantile(par_postJ1_gn2_th[,i], prob=0.025)
  par_post_lqJ1_C21[i] <- quantile(par_postJ1_mn1_th[,i], prob=0.025)
  par_post_lqJ1_C22[i] <- quantile(par_postJ1_mn2_th[,i], prob=0.025)
}

par_post_uqJ1_C11 <- rep(0, 15)
par_post_uqJ1_C12 <- rep(0, 15)
par_post_uqJ1_C21 <- rep(0, 15)
par_post_uqJ1_C22 <- rep(0, 15)

for(i in c(1:15)){
  par_post_uqJ1_C11[i] <- quantile(par_postJ1_gn1_th[,i], prob=0.975)
  par_post_uqJ1_C12[i] <- quantile(par_postJ1_gn2_th[,i], prob=0.975)
  par_post_uqJ1_C21[i] <- quantile(par_postJ1_mn1_th[,i], prob=0.975)
  par_post_uqJ1_C22[i] <- quantile(par_postJ1_mn2_th[,i], prob=0.975)
}

# - J=2
par_post_lqJ2_C11 <- rep(0, 22)
par_post_lqJ2_C12 <- rep(0, 22)
par_post_lqJ2_C21 <- rep(0, 22)
par_post_lqJ2_C22 <- rep(0, 22)

for(i in c(1:22)){
  par_post_lqJ2_C11[i] <- quantile(par_postJ2_gn1_th[,i], prob=0.025)
  par_post_lqJ2_C12[i] <- quantile(par_postJ2_gn2_th[,i], prob=0.025)
  par_post_lqJ2_C21[i] <- quantile(par_postJ2_mn1_th[,i], prob=0.025)
  par_post_lqJ2_C22[i] <- quantile(par_postJ2_mn2_th[,i], prob=0.025)
}

par_post_uqJ2_C11 <- rep(0, 22)
par_post_uqJ2_C12 <- rep(0, 22)
par_post_uqJ2_C21 <- rep(0, 22)
par_post_uqJ2_C22 <- rep(0, 22)

for(i in c(1:22)){
  par_post_uqJ2_C11[i] <- quantile(par_postJ2_gn1_th[,i], prob=0.975)
  par_post_uqJ2_C12[i] <- quantile(par_postJ2_gn2_th[,i], prob=0.975)
  par_post_uqJ2_C21[i] <- quantile(par_postJ2_mn1_th[,i], prob=0.975)
  par_post_uqJ2_C22[i] <- quantile(par_postJ2_mn2_th[,i], prob=0.975)
}

# - J=3
par_post_lqJ3_C11 <- rep(0, 29)
par_post_lqJ3_C12 <- rep(0, 29)
par_post_lqJ3_C21 <- rep(0, 29)
par_post_lqJ3_C22 <- rep(0, 29)

for(i in c(1:29)){
  par_post_lqJ3_C11[i] <- quantile(par_postJ3_gn1_th[,i], prob=0.025)
  par_post_lqJ3_C12[i] <- quantile(par_postJ3_gn2_th[,i], prob=0.025)
  par_post_lqJ3_C21[i] <- quantile(par_postJ3_mn1_th[,i], prob=0.025)
  par_post_lqJ3_C22[i] <- quantile(par_postJ3_mn2_th[,i], prob=0.025)
}

par_post_uqJ3_C11 <- rep(0, 29)
par_post_uqJ3_C12 <- rep(0, 29)
par_post_uqJ3_C21 <- rep(0, 29)
par_post_uqJ3_C22 <- rep(0, 29)

for(i in c(1:29)){
  par_post_uqJ3_C11[i] <- quantile(par_postJ3_gn1_th[,i], prob=0.975)
  par_post_uqJ3_C12[i] <- quantile(par_postJ3_gn2_th[,i], prob=0.975)
  par_post_uqJ3_C21[i] <- quantile(par_postJ3_mn1_th[,i], prob=0.975)
  par_post_uqJ3_C22[i] <- quantile(par_postJ3_mn2_th[,i], prob=0.975)
}

# - J=4
par_post_lqJ4_C11 <- rep(0, 36)
par_post_lqJ4_C12 <- rep(0, 36)
par_post_lqJ4_C21 <- rep(0, 36)
par_post_lqJ4_C22 <- rep(0, 36)

for(i in c(1:36)){
  par_post_lqJ4_C11[i] <- quantile(par_postJ4_gn1_th[,i], prob=0.025)
  par_post_lqJ4_C12[i] <- quantile(par_postJ4_gn2_th[,i], prob=0.025)
  par_post_lqJ4_C21[i] <- quantile(par_postJ4_mn1_th[,i], prob=0.025)
  par_post_lqJ4_C22[i] <- quantile(par_postJ4_mn2_th[,i], prob=0.025)
}

par_post_uqJ4_C11 <- rep(0, 36)
par_post_uqJ4_C12 <- rep(0, 36)
par_post_uqJ4_C21 <- rep(0, 36)
par_post_uqJ4_C22 <- rep(0, 36)

for(i in c(1:36)){
  par_post_uqJ4_C11[i] <- quantile(par_postJ4_gn1_th[,i], prob=0.975)
  par_post_uqJ4_C12[i] <- quantile(par_postJ4_gn2_th[,i], prob=0.975)
  par_post_uqJ4_C21[i] <- quantile(par_postJ4_mn1_th[,i], prob=0.975)
  par_post_uqJ4_C22[i] <- quantile(par_postJ4_mn2_th[,i], prob=0.975)
}

#################################### - Tables construction - #################################

library(xtable)

######### - J = 0

# - Transformation
par_post_meanJ0_C11[1] <- log(par_post_meanJ0_C11[1])
par_post_medianJ0_C11[1] <- log(par_post_medianJ0_C11[1])
par_post_lqJ0_C11[1] <- log(par_post_lqJ0_C11[1])
par_post_uqJ0_C11[1] <- log(par_post_uqJ0_C11[1])

J0_Table <- matrix(NA, 2, 5)
J0_Table <- as.data.frame(J0_Table)
J0_Table[,1] <- matrix(unlist(par_post_meanJ0_C11), nrow = 2, byrow = TRUE)
J0_Table[,2] <- matrix(unlist(par_post_medianJ0_C11), nrow = 2, byrow = TRUE)
J0_Table[,3] <- par_post_sdJ0_C11
J0_Table[,4] <- par_post_lqJ0_C11
J0_Table[,5] <- par_post_uqJ0_C11

colnames(J0_Table) <- c("Mean", "Median", "St.Dev", "2.5 Q", "97.5 Q")
xtable(J0_Table, digits = c(0,3,3,3,3,3))

# - Backtransformation
par_post_meanJ0_C11[1] <- exp(par_post_meanJ0_C11[1])
par_post_medianJ0_C11[1] <- exp(par_post_medianJ0_C11[1])
par_post_lqJ0_C11[1] <- exp(par_post_lqJ0_C11[1])
par_post_uqJ0_C11[1] <- exp(par_post_uqJ0_C11[1])

######### - J = 1

# - Transformation
par_post_meanJ1_C11[c(1,3)] <- log(par_post_meanJ1_C11[c(1,3)])
par_post_medianJ1_C11[c(1,3)] <- log(par_post_medianJ1_C11[c(1,3)])
par_post_lqJ1_C11[c(1,3)] <- log(par_post_lqJ1_C11[c(1,3)])
par_post_uqJ1_C11[c(1,3)] <- log(par_post_uqJ1_C11[c(1,3)])

J1_Table <- matrix(NA, 15, 5)

J1_Table[,1] <- matrix(unlist(par_post_meanJ1_C11), ncol = 15, byrow = TRUE)
J1_Table[,2] <- matrix(unlist(par_post_medianJ1_C11), ncol = 15, byrow = TRUE)
J1_Table[,3] <- par_post_sdJ1_C11
J1_Table[,4] <- par_post_lqJ1_C11
J1_Table[,5] <- par_post_uqJ1_C11

colnames(J1_Table) <- c("Mean", "Median", "St.Dev", "2.5 Q", "97.5 Q")
xtable(J1_Table, digits = c(0,3,3,3,3,3))

# - Backtransformation
par_post_meanJ1_C11[c(1,3)] <- exp(par_post_meanJ1_C11[c(1,3)])
par_post_medianJ1_C11[c(1,3)] <- exp(par_post_medianJ1_C11[c(1,3)])
par_post_lqJ1_C11[c(1,3)] <- exp(par_post_lqJ1_C11[c(1,3)])
par_post_uqJ1_C11[c(1,3)] <- exp(par_post_uqJ1_C11[c(1,3)])

######### - J = 2

# - Transformation
par_post_meanJ2_C11[c(1,3,4)] <- log(par_post_meanJ2_C11[c(1,3,4)])
par_post_medianJ2_C11[c(1,3,4)] <- log(par_post_medianJ2_C11[c(1,3,4)])
par_post_lqJ2_C11[c(1,3,4)] <- log(par_post_lqJ2_C11[c(1,3,4)])
par_post_uqJ2_C11[c(1,3,4)] <- log(par_post_uqJ2_C11[c(1,3,4)])

J2_Table <- matrix(NA, 22, 5)

J2_Table[,1] <- matrix(unlist(par_post_meanJ2_C11), nrow = 22, byrow = TRUE)
J2_Table[,2] <- matrix(unlist(par_post_medianJ2_C11), nrow = 22, byrow = TRUE)
J2_Table[,3] <- par_post_sdJ2_C11
J2_Table[,4] <- par_post_lqJ2_C11
J2_Table[,5] <- par_post_uqJ2_C11

colnames(J2_Table) <- c("Mean", "Median", "St.Dev", "2.5 Q", "97.5 Q")
rownames(J2_Table) <- c("alpha", "beta", "gamma_1", "gamma_2", "zeta0", "zeta1", "zeta2", "zeta3", "zeta4", "zeta5", "zeta6", "zeta7", "zeta8", "zeta9", "zeta10", "zeta11", "zeta12", "zeta13", "zeta14", "zeta15", "zeta16", "zeta17")
xtable(J2_Table, digits = c(0,3,3,3,3,3))

# - Transformation
par_post_meanJ2_C11[c(1,3,4)] <- exp(par_post_meanJ2_C11[c(1,3,4)])
par_post_medianJ2_C11[c(1,3,4)] <- exp(par_post_medianJ2_C11[c(1,3,4)])
par_post_lqJ2_C11[c(1,3,4)] <- exp(par_post_lqJ2_C11[c(1,3,4)])
par_post_uqJ2_C11[c(1,3,4)] <- exp(par_post_uqJ2_C11[c(1,3,4)])

######### - J = 3

# - Transformation
par_post_meanJ3_C11[c(1,3,4,5)] <- log(par_post_meanJ3_C11[c(1,3,4,5)])
par_post_medianJ3_C11[c(1,3,4,5)] <- log(par_post_medianJ3_C11[c(1,3,4,5)])
par_post_lqJ3_C11[c(1,3,4,5)] <- log(par_post_lqJ3_C11[c(1,3,4,5)])
par_post_uqJ3_C11[c(1,3,4,5)] <- log(par_post_uqJ3_C11[c(1,3,4,5)])

J3_Table <- matrix(NA, 29, 5)

J3_Table[,1] <- matrix(unlist(par_post_medianJ3_C11), nrow = 29, byrow = TRUE)
J3_Table[,2] <- matrix(unlist(par_post_medianJ3_C11), nrow = 29, byrow = TRUE)
J3_Table[,3] <- par_post_sdJ3_C11
J3_Table[,4] <- par_post_lqJ3_C11
J3_Table[,5] <- par_post_uqJ3_C11

colnames(J3_Table) <- c("Mean", "Median", "St.Dev", "2.5 Q", "97.5 Q")
xtable(J3_Table, digits = c(0,3,3,3,3,3))

# - Transformation
par_post_meanJ3_C11[c(1,3,4,5)] <- exp(par_post_meanJ3_C11[c(1,3,4,5)])
par_post_medianJ3_C11[c(1,3,4,5)] <- exp(par_post_medianJ3_C11[c(1,3,4,5)])
par_post_lqJ3_C11[c(1,3,4,5)] <- exp(par_post_lqJ3_C11[c(1,3,4,5)])
par_post_uqJ3_C11[c(1,3,4,5)] <- exp(par_post_uqJ3_C11[c(1,3,4,5)])

######### - J = 4

# - Transformation
par_post_meanJ4_C11[c(1,3,4,5,6)] <- log(par_post_meanJ4_C11[c(1,3,4,5,6)])
par_post_medianJ4_C11[c(1,3,4,5,6)] <- log(par_post_medianJ4_C11[c(1,3,4,5,6)])
par_post_lqJ4_C11[c(1,3,4,5,6)] <- log(par_post_lqJ4_C11[c(1,3,4,5,6)])
par_post_uqJ4_C11[c(1,3,4,5,6)] <- log(par_post_uqJ4_C11[c(1,3,4,5,6)])

J4_Table <- matrix(NA, 36, 5)

J4_Table[,1] <- matrix(unlist(par_post_meanJ4_C11), ncol = 36, byrow = TRUE)
J4_Table[,2] <- matrix(unlist(par_post_medianJ4_C11), ncol = 36, byrow = TRUE)
J4_Table[,3] <- par_post_sdJ4_C11
J4_Table[,4] <- par_post_lqJ4_C11
J4_Table[,5] <- par_post_uqJ4_C11

colnames(J4_Table) <- c("Mean", "Median", "St.Dev", "2.5 Q", "97.5 Q")
xtable(J4_Table, digits = c(0,3,3,3,3,3))

# - Transformation
par_post_meanJ4_C11[c(1,3,4,5,6)] <- exp(par_post_meanJ4_C11[c(1,3,4,5,6)])
par_post_medianJ4_C11[c(1,3,4,5,6)] <- exp(par_post_medianJ4_C11[c(1,3,4,5,6)])
par_post_lqJ4_C11[c(1,3,4,5,6)] <- exp(par_post_lqJ4_C11[c(1,3,4,5,6)])
par_post_uqJ4_C11[c(1,3,4,5,6)] <- exp(par_post_uqJ4_C11[c(1,3,4,5,6)])

##################################################

# - Transformation
par_post_meanJ4_C21[c(1,3,4,5,6)] <- log(par_post_meanJ4_C21[c(1,3,4,5,6)])
par_post_medianJ4_C21[c(1,3,4,5,6)] <- log(par_post_medianJ4_C21[c(1,3,4,5,6)])
par_post_lqJ4_C21[c(1,3,4,5,6)] <- log(par_post_lqJ4_C21[c(1,3,4,5,6)])
par_post_uqJ4_C21[c(1,3,4,5,6)] <- log(par_post_uqJ4_C21[c(1,3,4,5,6)])

J4_Table <- matrix(NA, 36, 5)

J4_Table[,1] <- matrix(unlist(par_post_meanJ4_C21), ncol = 36, byrow = TRUE)
J4_Table[,2] <- matrix(unlist(par_post_medianJ4_C21), ncol = 36, byrow = TRUE)
J4_Table[,3] <- par_post_sdJ4_C21
J4_Table[,4] <- par_post_lqJ4_C21
J4_Table[,5] <- par_post_uqJ4_C21

colnames(J4_Table) <- c("Mean", "Median", "St.Dev", "2.5 Q", "97.5 Q")
xtable(J4_Table, digits = c(0,3,3,3,3,3))

# - Transformation
par_post_meanJ4_C21[c(1,3,4,5,6)] <- exp(par_post_meanJ4_C21[c(1,3,4,5,6)])
par_post_medianJ4_C21[c(1,3,4,5,6)] <- exp(par_post_medianJ4_C21[c(1,3,4,5,6)])
par_post_lqJ4_C21[c(1,3,4,5,6)] <- exp(par_post_lqJ4_C21[c(1,3,4,5,6)])
par_post_uqJ4_C21[c(1,3,4,5,6)] <- exp(par_post_uqJ4_C21[c(1,3,4,5,6)])


################################## - Log-hazard function plot - ##############################

# - Alternative computation posterior mean

## - Transformation
par_postJ2_gn1_th[,c(1,3,4)] <- log(par_postJ2_gn1_th[,c(1,3,4)])
par_post_meanJ2_C11 <- colMeans(par_postJ2_gn1_th)
par_post_meanJ2_C11[c(1,3,4)] <- exp(par_post_meanJ2_C11[c(1,3,4)])

par_postJ2_gn2_th[,c(1,3,4)] <- log(par_postJ2_gn2_th[,c(1,3,4)])
par_post_meanJ2_C11 <- colMeans(par_postJ2_gn2_th)
par_post_meanJ2_C11[c(1,3,4)] <- exp(par_post_meanJ2_C11[c(1,3,4)])

## - Backtransformation
par_postJ2_gn2_th[,c(1,3,4)] <- exp(par_postJ2_gn2_th[,c(1,3,4)])



t_hf <- seq(0, 45, len=226)

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



phi0_ml <- exp(-11.41)
pgam_ml <- exp(-0.21)
pde1_ml <- exp(-0.22)
pde2_ml <- exp(-0.43)
beta_ml <- 0.108

# - Low benefit, Geo-dem 0
mu_bay_bLc0 <- (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
                z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
                z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) / 
               (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))))

mu_mle_bLc0 <- phi0_ml * exp(beta_ml * (60 + t_hf))

# - High benefit, Geo-dem 0

mu_bay_bHc0 <- (z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) / 
               (z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))))

mu_mle_bHc0 <- phi0_ml * exp(beta_ml * (60 + t_hf)) * pgam_ml

# - Low benefit, Geo-dem 1
mu_bay_bLc1 <- (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
                z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
                z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) / 
               (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))))

mu_mle_bLc1 <- phi0_ml * exp(beta_ml * (60 + t_hf)) * pde1_ml

# - High benefit, Geo-dem 1

mu_bay_bHc1 <- (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
                z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
                z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) / 
               (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))))

mu_mle_bHc1 <- phi0_ml * exp(beta_ml * (60 + t_hf)) * pgam_ml * pde1_ml

# - Low benefit, Geo-dem 2
mu_bay_bLc2 <- (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
                z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
                z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) / 
               (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))))

mu_mle_bLc2 <- phi0_ml * exp(beta_ml * (60 + t_hf)) * pde2_ml

# - High benefit, Geo-dem 2

mu_bay_bHc2 <- (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
                z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
                z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) / 
               (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))))

mu_mle_bHc2 <- phi0_ml * exp(beta_ml * (60 + t_hf)) * pgam_ml * pde2_ml


par(mfrow=c(2,3), mai=c(0.6,0.6, 0.2,0.2)) #, mai=c(0.7,0.7, 0.2,0.2) save as 
plot(60 + t_hf, log(mu_bay_bLc0), type = "l", col="black", lwd="2", xlab = "Age", ylab="log-hazard B=Low and C=0", ylim = c(-6, 2))
lines(60 + t_hf, log(mu_mle_bLc0), type = "l", col="black", lwd="2", lty=2)
lines(c(60:103), log(deaths_bl_geo_0/exposures_bl_geo_0), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bl_geo_0/exposures_bl_geo_0 + 1.96 * sqrt(deaths_bl_geo_0)/exposures_bl_geo_0), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bl_geo_0/exposures_bl_geo_0 - 1.96 * sqrt(deaths_bl_geo_0)/exposures_bl_geo_0), type = "l", col="black", lwd="1", lty=3)

plot(60 + t_hf, log(mu_bay_bLc1), type = "l", col="black", lwd="2", xlab = "Age", ylab="log-hazard B=Low and C=1", ylim = c(-6, 2))
lines(60 + t_hf, log(mu_mle_bLc1), type = "l", col="black", lwd="2", lty=2)
lines(c(60:103), log(deaths_bl_geo_1/exposures_bl_geo_1), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bl_geo_1/exposures_bl_geo_1 + 1.96 * sqrt(deaths_bl_geo_1)/exposures_bl_geo_1), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bl_geo_1/exposures_bl_geo_1 - 1.96 * sqrt(deaths_bl_geo_1)/exposures_bl_geo_1), type = "l", col="black", lwd="1", lty=3)

plot(60 + t_hf, log(mu_bay_bLc2), type = "l", col="black", lwd="2", xlab = "Age", ylab="log-hazard B=Low and C=2", ylim = c(-6, 2))
lines(60 + t_hf, log(mu_mle_bLc2), type = "l", col="black", lwd="2", lty=2)
lines(c(60:103), log(deaths_bl_geo_2/exposures_bl_geo_2), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bl_geo_2/exposures_bl_geo_2 + 1.96 * sqrt(deaths_bl_geo_2)/exposures_bl_geo_2), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bl_geo_2/exposures_bl_geo_2 - 1.96 * sqrt(deaths_bl_geo_2)/exposures_bl_geo_2), type = "l", col="black", lwd="1", lty=3)

plot(60 + t_hf, log(mu_bay_bHc0), type = "l", col="black", lwd="2", xlab = "Age", ylab="log-hazard B=High and C=0", ylim = c(-6, 2))
lines(60 + t_hf, log(mu_mle_bHc0), type = "l", col="black", lwd="2", lty=2)
lines(c(60:103), log(deaths_bh_geo_0/exposures_bh_geo_0), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bh_geo_0/exposures_bh_geo_0 + 1.96 * sqrt(deaths_bh_geo_0)/exposures_bh_geo_0), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bh_geo_0/exposures_bh_geo_0 - 1.96 * sqrt(deaths_bh_geo_0)/exposures_bh_geo_0), type = "l", col="black", lwd="1", lty=3)

plot(60 + t_hf, log(mu_bay_bHc1), type = "l", col="black", lwd="2", xlab = "Age", ylab="log-hazard B=High and C=1", ylim = c(-6, 2))
lines(60 + t_hf, log(mu_mle_bHc1), type = "l", col="black", lwd="2", lty=2)
lines(c(60:103), log(deaths_bh_geo_1/exposures_bh_geo_1), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bh_geo_1/exposures_bh_geo_1 + 1.96 * sqrt(deaths_bh_geo_1)/exposures_bh_geo_1), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bh_geo_1/exposures_bh_geo_1 - 1.96 * sqrt(deaths_bh_geo_1)/exposures_bh_geo_1), type = "l", col="black", lwd="1", lty=3)

plot(60 + t_hf, log(mu_bay_bHc2), type = "l", col="black", lwd="2", xlab = "Age", ylab="log-hazard B=High and C=2", ylim = c(-6, 2))
lines(60 + t_hf, log(mu_mle_bHc2), type = "l", col="black", lwd="2", lty=2)
lines(c(60:103), log(deaths_bh_geo_2/exposures_bh_geo_2), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bh_geo_2/exposures_bh_geo_2 + 1.96 * sqrt(deaths_bh_geo_2)/exposures_bh_geo_2), type = "l", col="black", lwd="1", lty=3)
lines(c(60:103), log(deaths_bh_geo_2/exposures_bh_geo_2 - 1.96 * sqrt(deaths_bh_geo_2)/exposures_bh_geo_2), type = "l", col="black", lwd="1", lty=3)

plot(60 + t_hf, log(mu_bay_bLc0 / mu_mle_bLc0), type = "l", lty=1, col="black", lwd="1", xlab = "Age", ylab="log-hazard bayes-mle diff", ylim=c(-0.22, 0.38)) #, ylim=c(-0.22, 0.38)
lines(60 + t_hf, log(mu_bay_bHc0 / mu_mle_bHc0), type = "l", lty=2, col="black", lwd="1")
lines(60 + t_hf, log(mu_bay_bLc1 / mu_mle_bLc1), type = "l", lty=3, col="black", lwd="1")
lines(60 + t_hf, log(mu_bay_bHc1 / mu_mle_bHc1), type = "l", lty=4, col="black", lwd="1")
lines(60 + t_hf, log(mu_bay_bLc2 / mu_mle_bLc2), type = "l", lty=5, col="black", lwd="1")
lines(60 + t_hf, log(mu_bay_bHc2 / mu_mle_bHc2), type = "l", lty=6, col="black", lwd="1")
legend('topleft', c("Ben = Low, C = 0","Ben = High, C = 0", "Ben = Low, C = 1","Ben = High, C = 1", "Ben = Low, C = 2","Ben = High, C = 2"), lty=c(1,2,3,4,5,6), bty='n', cex=0.9, lwd=1)

plot(60 + t_hf, (mu_bay_bLc0 - mu_mle_bLc0) / mu_mle_bLc0, type = "l", lty=1, col="black", lwd="1", xlab = "Age", ylab="log-hazard bayes-mle diff", ylim=c(-0.22, 0.38)) #, ylim=c(-0.22, 0.38)
lines(60 + t_hf, (mu_bay_bHc0 - mu_mle_bHc0) / mu_mle_bHc0, type = "l", lty=2, col="black", lwd="1")
lines(60 + t_hf, (mu_bay_bLc1 - mu_mle_bLc1) / mu_mle_bLc1, type = "l", lty=3, col="black", lwd="1")
lines(60 + t_hf, (mu_bay_bHc1 - mu_mle_bHc1) / mu_mle_bHc1, type = "l", lty=4, col="black", lwd="1")
lines(60 + t_hf, (mu_bay_bLc2 - mu_mle_bLc2) / mu_mle_bLc2, type = "l", lty=5, col="black", lwd="1")
lines(60 + t_hf, (mu_bay_bHc2 - mu_mle_bHc2) / mu_mle_bHc2, type = "l", lty=6, col="black", lwd="1")
legend('topleft', c("Ben = Low, C = 0","Ben = High, C = 0", "Ben = Low, C = 1","Ben = High, C = 1", "Ben = Low, C = 2","Ben = High, C = 2"), lty=c(1,2,3,4,5,6), bty='n', cex=0.9, lwd=1)



mu_bay_bLc0 - mu_mle_bLc0
mu_bay_bHc0 - mu_mle_bHc0
mu_bay_bLc1 - mu_mle_bLc1
mu_bay_bHc1 - mu_mle_bHc1
mu_bay_bLc2 - mu_mle_bLc2
mu_bay_bHc2 - mu_mle_bHc2


plot(60 + t_hf, mu_bay_bLc0 - mu_mle_bLc0, type = "l", lty=1, col="black", lwd="1", xlab = "Age", ylab="Haz. fun Bayes - Haz. fun. MLE", ylim=c(-0.5, 0.5)) #
lines(60 + t_hf, mu_bay_bHc0 - mu_mle_bHc0, type = "l", lty=2, col="black", lwd="1")
lines(60 + t_hf, mu_bay_bLc1 - mu_mle_bLc1, type = "l", lty=3, col="black", lwd="1")
lines(60 + t_hf, mu_bay_bHc1 - mu_mle_bHc1, type = "l", lty=4, col="black", lwd="1")
lines(60 + t_hf, mu_bay_bLc2 - mu_mle_bLc2, type = "l", lty=5, col="black", lwd="1")
lines(60 + t_hf, mu_bay_bHc2 - mu_mle_bHc2, type = "l", lty=6, col="black", lwd="1")
legend('topleft', c("Ben = Low, C = 0","Ben = High, C = 0", "Ben = Low, C = 1","Ben = High, C = 1", "Ben = Low, C = 2","Ben = High, C = 2"), lty=c(1,2,3,4,5,6), bty='n', cex=0.9, lwd=1)


##################################### - Analysis of the hazard function - ###############################################

pdf_g0_cond_bLc0 <- z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) / 
                   (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g0_cond_bLc1 <- z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) / 
                   (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g0_cond_bLc2 <- z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) / 
                   (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g1_cond_bLc0 <- z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) / 
                   (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g1_cond_bLc1 <- z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) / 
                   (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g1_cond_bLc2 <- z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) / 
                   (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g2_cond_bLc0 <- z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) / 
                   (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g2_cond_bLc1 <- z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) / 
                   (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g2_cond_bLc2 <- z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) / 
                   (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )


pdf_g0_cond_bHc0 <-  z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) / 
                   ( z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g0_cond_bHc1 <- z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) / 
                   (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g0_cond_bHc2 <- z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) / 
                   (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g1_cond_bHc0 <- z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) / 
                   (z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g1_cond_bHc1 <- z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) / 
                   (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g1_cond_bHc2 <- z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) / 
                   (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g2_cond_bHc0 <- z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) / 
                   ( z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g2_cond_bHc1 <- z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) / 
                   (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

pdf_g2_cond_bHc2 <- z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) / 
                   (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                    z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                    z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) )

# Calculation of posterior expectation of G

exp_g_cond_bLc0 <- pdf_g0_cond_bLc0 + phi1_p * pdf_g1_cond_bLc0 + (phi1_p * phi2_p) * pdf_g2_cond_bLc0
exp_g_cond_bLc1 <- pdf_g0_cond_bLc1 + phi1_p * pdf_g1_cond_bLc1 + (phi1_p * phi2_p) * pdf_g2_cond_bLc1
exp_g_cond_bLc2 <- pdf_g0_cond_bLc2 + phi1_p * pdf_g1_cond_bLc2 + (phi1_p * phi2_p) * pdf_g2_cond_bLc2
exp_g_cond_bHc0 <- pdf_g0_cond_bHc0 + phi1_p * pdf_g1_cond_bHc0 + (phi1_p * phi2_p) * pdf_g2_cond_bHc0
exp_g_cond_bHc1 <- pdf_g0_cond_bHc1 + phi1_p * pdf_g1_cond_bHc1 + (phi1_p * phi2_p) * pdf_g2_cond_bHc1
exp_g_cond_bHc2 <- pdf_g0_cond_bHc2 + phi1_p * pdf_g1_cond_bHc2 + (phi1_p * phi2_p) * pdf_g2_cond_bHc2

# Calculation of posterior variance of G

var_g_cond_bLc0 <- pdf_g0_cond_bLc0 + (phi1_p^2) * pdf_g1_cond_bLc0 + ((phi1_p * phi2_p)^2) * pdf_g2_cond_bLc0 - (exp_g_cond_bLc0 ^ 2)
var_g_cond_bLc1 <- pdf_g0_cond_bLc1 + (phi1_p^2) * pdf_g1_cond_bLc1 + ((phi1_p * phi2_p)^2) * pdf_g2_cond_bLc1 - (exp_g_cond_bLc1 ^ 2)
var_g_cond_bLc2 <- pdf_g0_cond_bLc2 + (phi1_p^2) * pdf_g1_cond_bLc2 + ((phi1_p * phi2_p)^2) * pdf_g2_cond_bLc2 - (exp_g_cond_bLc2 ^ 2)
var_g_cond_bHc0 <- pdf_g0_cond_bHc0 + (phi1_p^2) * pdf_g1_cond_bHc0 + ((phi1_p * phi2_p)^2) * pdf_g2_cond_bHc0 - (exp_g_cond_bHc0 ^ 2)
var_g_cond_bHc1 <- pdf_g0_cond_bHc1 + (phi1_p^2) * pdf_g1_cond_bHc1 + ((phi1_p * phi2_p)^2) * pdf_g2_cond_bHc1 - (exp_g_cond_bHc1 ^ 2)
var_g_cond_bHc2 <- pdf_g0_cond_bHc2 + (phi1_p^2) * pdf_g1_cond_bHc2 + ((phi1_p * phi2_p)^2) * pdf_g2_cond_bHc2 - (exp_g_cond_bHc2 ^ 2)

# - deriv log_mu wrt t

der_t_bLc0 <- beta_p - (var_g_cond_bLc0 / exp_g_cond_bLc0) * (exp(beta_p * (t_hf + 60 - 77.5)))
der_t_bLc1 <- beta_p - (var_g_cond_bLc1 / exp_g_cond_bLc1) * (exp(beta_p * (t_hf + 60 - 77.5)))
der_t_bLc2 <- beta_p - (var_g_cond_bLc2 / exp_g_cond_bLc2) * (exp(beta_p * (t_hf + 60 - 77.5)))
der_t_bHc0 <- beta_p - (var_g_cond_bHc0 / exp_g_cond_bHc0) * (exp(beta_p * (t_hf + 60 - 77.5)))
der_t_bHc1 <- beta_p - (var_g_cond_bHc1 / exp_g_cond_bHc1) * (exp(beta_p * (t_hf + 60 - 77.5)))
der_t_bHc2 <- beta_p - (var_g_cond_bHc2 / exp_g_cond_bHc2) * (exp(beta_p * (t_hf + 60 - 77.5)))

# - Plot

plot(60 + t_hf, der_t_bLc0, type = "l", lty=1, col="black", lwd="1", xlab = "Age") #, ylim=c(0.05, 0.115)
abline(h=beta_ml)
lines(60 + t_hf, der_t_bLc1, type = "l", lty=1, lwd="1", xlab = "Age", col="red")
lines(60 + t_hf, der_t_bLc2, type = "l", lty=1, lwd="1", xlab = "Age", col="blue")
lines(60 + t_hf, der_t_bHc0, type = "l", lty=1, lwd="1", xlab = "Age", col="green")
lines(60 + t_hf, der_t_bHc1, type = "l", lty=1, lwd="1", xlab = "Age", col="grey")
lines(60 + t_hf, der_t_bHc2, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
text(locator(), labels = c("BL-C0", "BL-C1", "BL-C2", "BH-C0", "BH-C1", "BH-C2"))


plot(60 + t_hf, der_t_bLc1, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bLc2, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bHc0, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bHc1, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bHc2, type = "l", lty=1, col="black", lwd="1", xlab = "Age")

##################################### - Analysis of the hazard function_v2 - ###############################################

mu_g0 <- exp(beta_p * (60 - 77.5 + t_hf)) * phi0_p
mu_g1 <- exp(beta_p * (60 - 77.5 + t_hf)) * phi0_p * phi1_p
mu_g2 <- exp(beta_p * (60 - 77.5 + t_hf)) * phi0_p * phi1_p * phi2_p


der_t_bLc0 <-( z0_p * mu_g0 * (beta_p - mu_g0) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z3_p * mu_g1 * (beta_p - mu_g1) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z6_p * mu_g2 * (beta_p - mu_g2) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / 
              (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) + mu_bay_bLc0

der_t_bLc1 <-( z1_p * mu_g0 * (beta_p - mu_g0) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z4_p * mu_g1 * (beta_p - mu_g1) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z7_p * mu_g2 * (beta_p - mu_g2) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / 
              (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) + mu_bay_bLc1


der_t_bLc2 <-( z2_p * mu_g0 * (beta_p - mu_g0) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z5_p * mu_g1 * (beta_p - mu_g1) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z8_p * mu_g2 * (beta_p - mu_g2) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / 
              (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) + mu_bay_bLc2

der_t_bHc0 <-( z9_p  * mu_g0 * (beta_p - mu_g0) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z12_p * mu_g1 * (beta_p - mu_g1) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z15_p * mu_g2 * (beta_p - mu_g2) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / 
              (z9_p  * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) + mu_bay_bHc0

der_t_bHc1 <-( z10_p * mu_g0 * (beta_p - mu_g0) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z13_p * mu_g1 * (beta_p - mu_g1) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z16_p * mu_g2 * (beta_p - mu_g2) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / 
              (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) + mu_bay_bHc1


der_t_bHc2 <-( z11_p * mu_g0 * (beta_p - mu_g0) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
               z14_p * mu_g1 * (beta_p - mu_g1) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
               z17_p * mu_g2 * (beta_p - mu_g2) * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / 
              (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) * phi0_p                   * exp(beta_p * (60 + t_hf - 77.5)) + 
               z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p          * exp(beta_p * (60 + t_hf - 77.5)) + 
               z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5))) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 + t_hf - 77.5))) + mu_bay_bHc2


# - Plot

plot(60 + t_hf, der_t_bLc0, type = "l", lty=1, col="black", lwd="1", xlab = "Age") #, ylim=c(0.05, 0.115)
abline(h=beta_ml)
abline(h=beta_p, col="red")
lines(60 + t_hf, der_t_bLc1, type = "l", lty=1, lwd="1", xlab = "Age", col="red")
lines(60 + t_hf, der_t_bLc2, type = "l", lty=1, lwd="1", xlab = "Age", col="blue")
lines(60 + t_hf, der_t_bHc0, type = "l", lty=1, lwd="1", xlab = "Age", col="green")
lines(60 + t_hf, der_t_bHc1, type = "l", lty=1, lwd="1", xlab = "Age", col="grey")
lines(60 + t_hf, der_t_bHc2, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
text(locator(), labels = c("BL-C0", "BL-C1", "BL-C2", "BH-C0", "BH-C1", "BH-C2"))


plot(60 + t_hf, der_t_bLc1, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bLc2, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bHc0, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bHc1, type = "l", lty=1, col="black", lwd="1", xlab = "Age")
plot(60 + t_hf, der_t_bHc2, type = "l", lty=1, col="black", lwd="1", xlab = "Age")

############################################## - Analysis of the Survival function - ###########################################

t_hf <- seq(0, 45, len=226)

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



phi0_ml <- exp(-11.41)
pgam_ml <- exp(-0.21)
pde1_ml <- exp(-0.22)
pde2_ml <- exp(-0.43)
beta_ml <- 0.108

# - Low benefit, Geo-dem 0
surv_bay_bLc0 <- (z0_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                  z3_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                  z6_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / (z0_p + z3_p + z6_p)

surv_mle_bLc0 <- exp(-((exp(beta_ml * t_hf) - 1) / beta_ml) * phi0_ml * exp(beta_ml * 60))

# - High benefit, Geo-dem 0

surv_bay_bHc0 <- (  z9_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                   z12_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                   z15_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / (z9_p + z12_p + z15_p)

surv_mle_bHc0 <- exp(-((exp(beta_ml * t_hf) - 1) / beta_ml) * phi0_ml * exp(beta_ml * 60) * pgam_ml)

# - Low benefit, Geo-dem 1
surv_bay_bLc1 <- (z1_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                  z4_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                  z7_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / (z1_p + z4_p + z7_p)

surv_mle_bLc1 <- exp(-((exp(beta_ml * t_hf) - 1) / beta_ml) * phi0_ml * exp(beta_ml * 60) * pde1_ml)

# - High benefit, Geo-dem 1

surv_bay_bHc1 <- (z10_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                  z13_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                  z16_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / (z10_p + z13_p + z16_p)

surv_mle_bHc1 <- exp(-((exp(beta_ml * t_hf) - 1) / beta_ml) * phi0_ml * exp(beta_ml * 60) * pgam_ml * pde1_ml)

# - Low benefit, Geo-dem 2
surv_bay_bLc2 <- (z2_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                  z5_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                  z8_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / (z2_p + z5_p + z8_p)

surv_mle_bLc2 <- exp(-((exp(beta_ml * t_hf) - 1) / beta_ml) * phi0_ml * exp(beta_ml * 60) * pde2_ml)

# - High benefit, Geo-dem 2

surv_bay_bHc2 <- (z11_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p                   * exp(beta_p * (60 - 77.5))) + 
                  z14_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p          * exp(beta_p * (60 - 77.5))) + 
                  z17_p * exp(-((exp(beta_p * t_hf) - 1) / beta_p) * phi0_p * phi1_p * phi2_p * exp(beta_p * (60 - 77.5)))) / (z11_p + z14_p + z17_p)

surv_mle_bHc2 <- exp(-((exp(beta_ml * t_hf) - 1) / beta_ml) * phi0_ml * exp(beta_ml * 60) * pgam_ml * pde2_ml)


par(mfrow=c(2,3), mai=c(0.6,0.6, 0.2,0.2)) #, mai=c(0.7,0.7, 0.2,0.2) save as 9 x 6
plot(60 + t_hf, surv_bay_bLc0, type = "l", col="black", lwd="2", xlab = "Age", ylab="S(t) B=Low and C=0", ylim = c(0, 1))
lines(60 + t_hf, surv_mle_bLc0, type = "l", col="black", lwd="2", lty=2)
# lines(KM_Table_ben_l_geo_0$ExitAge, KM_Table_ben_l_geo_0$S_t, type="l", lty=3)
lines(KM_Table_ben_l_geo_0$ExitAge, KM_Table_ben_l_geo_0$S_t + 1.96 * KM_Table_ben_l_geo_0$sd_S_t, type="l", lty=3)
lines(KM_Table_ben_l_geo_0$ExitAge, KM_Table_ben_l_geo_0$S_t - 1.96 * KM_Table_ben_l_geo_0$sd_S_t, type="l", lty=3)

plot(60 + t_hf, surv_bay_bLc1, type = "l", col="black", lwd="2", xlab = "Age", ylab="S(t) B=Low and C=1", ylim = c(0, 1))
lines(60 + t_hf, surv_mle_bLc1, type = "l", col="black", lwd="2", lty=2)
#lines(KM_Table_ben_l_geo_1$ExitAge, KM_Table_ben_l_geo_1$S_t, type="l", lty=3)
lines(KM_Table_ben_l_geo_1$ExitAge, KM_Table_ben_l_geo_1$S_t + 1.96 * KM_Table_ben_l_geo_1$sd_S_t, type="l", lty=3)
lines(KM_Table_ben_l_geo_1$ExitAge, KM_Table_ben_l_geo_1$S_t - 1.96 * KM_Table_ben_l_geo_1$sd_S_t, type="l", lty=3)

plot(60 + t_hf, surv_bay_bLc2, type = "l", col="black", lwd="2", xlab = "Age", ylab="S(t) B=Low and C=2", ylim = c(0, 1))
lines(60 + t_hf, surv_mle_bLc2, type = "l", col="black", lwd="2", lty=2)
#lines(KM_Table_ben_l_geo_2$ExitAge, KM_Table_ben_l_geo_2$S_t, type="l", lty=3)
lines(KM_Table_ben_l_geo_2$ExitAge, KM_Table_ben_l_geo_2$S_t + 1.96 * KM_Table_ben_l_geo_2$sd_S_t, type="l", lty=3)
lines(KM_Table_ben_l_geo_2$ExitAge, KM_Table_ben_l_geo_2$S_t - 1.96 * KM_Table_ben_l_geo_2$sd_S_t, type="l", lty=3)

plot(60 + t_hf, surv_bay_bHc0, type = "l", col="black", lwd="2", xlab = "Age", ylab="S(t) B=High and C=0", ylim = c(0, 1))
lines(60 + t_hf, surv_mle_bHc0, type = "l", col="black", lwd="2", lty=2)
#lines(KM_Table_ben_h_geo_0$ExitAge, KM_Table_ben_h_geo_0$S_t, type="l", lty=3)
lines(KM_Table_ben_h_geo_0$ExitAge, KM_Table_ben_h_geo_0$S_t + 1.96 * KM_Table_ben_h_geo_0$sd_S_t, type="l", lty=3)
lines(KM_Table_ben_h_geo_0$ExitAge, KM_Table_ben_h_geo_0$S_t - 1.96 * KM_Table_ben_h_geo_0$sd_S_t, type="l", lty=3)

plot(60 + t_hf, surv_bay_bHc1, type = "l", col="black", lwd="2", xlab = "Age", ylab="S(t) B=High and C=1", ylim = c(0, 1))
lines(60 + t_hf, surv_mle_bHc1, type = "l", col="black", lwd="2", lty=2)
#lines(KM_Table_ben_h_geo_1$ExitAge, KM_Table_ben_h_geo_1$S_t, type="l", lty=3)
lines(KM_Table_ben_h_geo_1$ExitAge, KM_Table_ben_h_geo_1$S_t + 1.96 * KM_Table_ben_h_geo_1$sd_S_t, type="l", lty=3)
lines(KM_Table_ben_h_geo_1$ExitAge, KM_Table_ben_h_geo_1$S_t - 1.96 * KM_Table_ben_h_geo_1$sd_S_t, type="l", lty=3)

plot(60 + t_hf, surv_bay_bHc2, type = "l", col="black", lwd="2", xlab = "Age", ylab="S(t) B=High and C=2", ylim = c(0, 1))
lines(60 + t_hf, surv_mle_bHc2, type = "l", col="black", lwd="2", lty=2)
#lines(KM_Table_ben_h_geo_2$ExitAge, KM_Table_ben_h_geo_2$S_t, type="l", lty=3)
lines(KM_Table_ben_h_geo_2$ExitAge, KM_Table_ben_h_geo_2$S_t + 1.96 * KM_Table_ben_h_geo_2$sd_S_t, type="l", lty=3)
lines(KM_Table_ben_h_geo_2$ExitAge, KM_Table_ben_h_geo_2$S_t - 1.96 * KM_Table_ben_h_geo_2$sd_S_t, type="l", lty=3)


# Max

max(abs(surv_bay_bLc0 - surv_mle_bLc0))
max(abs(surv_bay_bLc1 - surv_mle_bLc1))
max(abs(surv_bay_bLc2 - surv_mle_bLc2))
max(abs(surv_bay_bHc0 - surv_mle_bHc0))
max(abs(surv_bay_bHc1 - surv_mle_bHc1))
max(abs(surv_bay_bHc2 - surv_mle_bHc2))



####################### - Var/cov of zeta - #############################

cor_zeta<-round(cor(par_postJ1_gn1_th[,c(4:15)]),2)

diag(cor_zeta) <- -1

max(cor_zeta)


pr_BL <- (par_postJ2_gn1_th[,5] + par_postJ2_gn1_th[,6] + par_postJ2_gn1_th[,7] + par_postJ2_gn1_th[,8] + par_postJ2_gn1_th[,9] + par_postJ2_gn1_th[,10] + par_postJ2_gn1_th[,11] + par_postJ2_gn1_th[,12] + par_postJ2_gn1_th[,13])

pr_C0 <- (par_postJ2_gn1_th[,5] + par_postJ2_gn1_th[,8] + par_postJ2_gn1_th[,11] + par_postJ2_gn1_th[,14] + par_postJ2_gn1_th[,17] + par_postJ2_gn1_th[,20])

pr_G0 <- (par_postJ2_gn1_th[,5] + par_postJ2_gn1_th[,6] + par_postJ2_gn1_th[,7] + par_postJ2_gn1_th[,14] + par_postJ2_gn1_th[,15] + par_postJ2_gn1_th[,16])

#cor_bl_c0_g0 <- 
  cor(cbind(pr_BL, pr_C0, pr_G0))
  

pr_BL_C0 <- (par_postJ2_gn1_th[,5] + par_postJ2_gn1_th[,8] + par_postJ2_gn1_th[,11])
  
pr_BL_G0 <- (par_postJ2_gn1_th[,5] + par_postJ2_gn1_th[,6] + par_postJ2_gn1_th[,7])

pr_C0_G0 <- (par_postJ2_gn1_th[,5] + par_postJ2_gn1_th[,14])

cor(cbind(pr_BL_C0, pr_BL_G0, pr_C0_G0))



pr_BH <- (par_postJ2_gn1_th[,14] + par_postJ2_gn1_th[,15] + par_postJ2_gn1_th[,16] + par_postJ2_gn1_th[,17] + par_postJ2_gn1_th[,18] + par_postJ2_gn1_th[,19] + par_postJ2_gn1_th[,20] + par_postJ2_gn1_th[,21] + par_postJ2_gn1_th[,22])

pr_C2 <- (par_postJ2_gn1_th[,7] + par_postJ2_gn1_th[,10] + par_postJ2_gn1_th[,13] + par_postJ2_gn1_th[,16] + par_postJ2_gn1_th[,19] + par_postJ2_gn1_th[,22])

pr_G2 <- (par_postJ2_gn1_th[,11] + par_postJ2_gn1_th[,12] + par_postJ2_gn1_th[,13] + par_postJ2_gn1_th[,20] + par_postJ2_gn1_th[,21] + par_postJ2_gn1_th[,22])

#cor_bh_c2_g1 <- 
  cor(cbind(pr_BH, pr_C2, pr_G2))
  
  
pr_BH_C2 <- (par_postJ2_gn1_th[,16] + par_postJ2_gn1_th[,19] + par_postJ2_gn1_th[,22])
  
pr_BH_G2 <- (par_postJ2_gn1_th[,20] + par_postJ2_gn1_th[,21] + par_postJ2_gn1_th[,22])
  
pr_C2_G2 <- (par_postJ2_gn1_th[,13] + par_postJ2_gn1_th[,22])
  
#cor_bh_c2_g1 <- 
cor(cbind(pr_BH_C2, pr_BH_G2, pr_C2_G2))
