#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation of Influenza 2009                 #
#           Copyright Oswaldo Gressani 2025. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original data can be found in the coarseDataTools package
library("xlsx")
library("EpiLPS")
library("coarseDataTools")
rm(list = ls())

# Data for 2009 influenza
data("nycH1N1")

set.seed(2009)
n <- nrow(nycH1N1)
EL <- c() # Left bound of exposure
ER <- c() # Right bound of exposure
SL <- c() # Left bound of symptom onset
SR <- c() # Right bound of symptom onset
maxx <- 1.5

for(i in 1:n){
tEL_temp <- nycH1N1$EL[i] + runif(1, min = 0, max = maxx)
tER_temp <- nycH1N1$ER[i] + runif(1, min = 0, max = maxx)
tSL_temp <- nycH1N1$SL[i] + runif(1, min = 0, max = maxx)
tSR_temp <- nycH1N1$SR[i] + runif(1, min = 0, max = maxx)
  while(tEL_temp >= tER_temp | tSL_temp >= tSR_temp | tSL_temp <= tER_temp){
    tEL_temp <- nycH1N1$EL[i] + runif(1, min = 0, max = maxx)
    tER_temp <- nycH1N1$ER[i] + runif(1, min = 0, max = maxx)
    tSL_temp <- nycH1N1$SL[i] + runif(1, min = 0, max = maxx)
    tSR_temp <- nycH1N1$SR[i] + runif(1, min = 0, max = maxx)
  }
EL[i] <- tEL_temp
ER[i] <- tER_temp
SL[i] <- tSL_temp
SR[i] <- tSR_temp
}

tL <- SL-ER
tR <- SR-EL

data <- data.frame(tL = tL, tR = tR)
n <- nrow(data)
write.xlsx(round(data,3), file = "Data_continuous/Pathogen5-Influenza2009.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 15)

# Extract estimates
Pathogen5_Influenza2009_Estimates <- matrix(0, nrow = 2, ncol = 3)
colnames(Pathogen5_Influenza2009_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen5_Influenza2009_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)")
Pathogen5_Influenza2009_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)

# Write estimates in Estimates folder
write.xlsx(Pathogen5_Influenza2009_Estimates, 
           file = "Estimates/Pathogen5_Influenza2009_Estimates.xlsx",
           sheetName = "Influenza2009-Semiparametric")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)
















