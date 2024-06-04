#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation Zika virus                        #
#           Copyright Oswaldo Gressani 2024. All rights reserved.              #
#------------------------------------------------------------------------------#


library("xlsx")
library("EpiLPS")
rm(list = ls())

set.seed(1234)

## Data
# Original data downloaded from:
# https://github.com/HopkinsIDD/ZikaLitReview/blob/master/ZikaLitReviewData.csv
# Last accessed 2024-05-03 at 17h01 CEST.
dataraw <- read.csv(file = "Data/Pathogen6-Zika.csv")
ELraw <- dataraw$EL
ERraw <- dataraw$ER
SLraw <- dataraw$SL
SRraw <- dataraw$SR

databounds <- cbind(ELraw, ERraw, SLraw, SRraw)
# Filter 1 (remove NAs)
keeprows <- c(1,4,10,14,19,22,
              28,33,38,48,55,
              62,67,72,75,80,
              85,90,92,96,109,
              111,117,121,127,135)
databounds <- databounds[keeprows,]
databounds <- as.data.frame(databounds)

EL <- c() # Left bound of exposure
ER <- c() # Right bound of exposure
SL <- c() # Left bound of symptom onset
SR <- c() # Right bound of symptom onset
indexspec <- which(databounds$SLraw<databounds$ERraw)

for(i in 1:nrow(databounds)){
  if(any(indexspec == i)){
    maxx <- 1.5
  } else{
    maxx <- 1
  }
  tEL_temp <- databounds$ELraw[i] + runif(1, min = 0, max = maxx)
  tER_temp <- databounds$ERraw[i] + runif(1, min = 0, max = maxx)
  tSL_temp <- databounds$SLraw[i] + runif(1, min = 0, max = maxx)
  tSR_temp <- databounds$SRraw[i] + runif(1, min = 0, max = maxx)
  while(tEL_temp >= tER_temp | tSL_temp >= tSR_temp | tSL_temp <= tER_temp){
    tEL_temp <- databounds$ELraw[i] + runif(1, min = 0, max = maxx)
    tER_temp <- databounds$ERraw[i] + runif(1, min = 0, max = maxx)
    tSL_temp <- databounds$SLraw[i] + runif(1, min = 0, max = maxx)
    tSR_temp <- databounds$SRraw[i] + runif(1, min = 0, max = maxx)
  }
  EL[i] <- tEL_temp
  ER[i] <- tER_temp
  SL[i] <- tSL_temp
  SR[i] <- tSR_temp
}

tL <- SL-ER
tR <- SR-EL

data <- data.frame(tL = tL, tR = tR)
data <- data[-which((data$tR-data$tL)>20),] # Remove obs having incub window > 20
n <- nrow(data)
write.xlsx(round(data,3), file = "Data_continuous/Pathogen6-Zika.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 25)

# Extract estimates
Pathogen6_Zika_Estimates <- matrix(0, nrow = 3, ncol = 3)
colnames(Pathogen6_Zika_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen6_Zika_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)",
                                        "95% CI of incubation time (days)")
Pathogen6_Zika_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)
tdom <- incubfit$tg
fhat <- incubfit$ftg
dt <- tdom[2] - tdom[1]
Fhat <- cumsum(fhat * dt)
Pathogen6_Zika_Estimates[3,] <- c(NA,
       round(c(incubfit$tg[sum(Fhat<=0.025)],incubfit$tg[sum(Fhat<=0.975)]),1))

# Write estimates in Estimates folder
write.xlsx(Pathogen6_Zika_Estimates, 
           file = "Estimates/Pathogen6_Zika_Estimates.xlsx",
           sheetName = "Zika-Semiparametric")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)













