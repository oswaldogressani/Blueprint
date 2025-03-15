#------------------------------------------------------------------------------#
#                 EpiLPS incubation estimation MERS                            #
#             Oswaldo Gressani 2025. All rights reserved.                      #
#------------------------------------------------------------------------------#

# Original data can be found in Cauchemez et al. (2014)
# https://www.sciencedirect.com/science/article/pii/S1473309913703049

library("xlsx")
library("EpiLPS")

## Data
tL <- c(6,1,9,3,4,3,4)
tR <- c(9,4,12,4,4,3,4)
n <- length(tL)

set.seed(1234)

for(i in 1:n){
  tL_temp <- tL[i] + runif(1)
  tR_temp <- tR[i] + runif(1)
  while(tL_temp >= tR_temp){
    tL_temp <- tL[i] + runif(1)
    tR_temp <- tR[i] + runif(1)
  }
  tL[i] <- tL_temp
  tR[i] <- tR_temp
}

data <- data.frame(tL = tL, tR = tR)
write.xlsx(round(data,3), file = "Data_continuous/Pathogen1-MERS.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 20)

# Extract estimates
Pathogen1_MERS_Estimates <- matrix(0, nrow = 2, ncol = 3)
colnames(Pathogen1_MERS_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen1_MERS_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)")
Pathogen1_MERS_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)


# Write estimates in Estimates folder
write.xlsx(Pathogen1_MERS_Estimates, 
           file = "Estimates/Pathogen1_MERS_Estimates.xlsx",
           sheetName = "MERS-LogNormal")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)


