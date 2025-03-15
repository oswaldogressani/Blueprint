#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation of SARS 2003                      #
#           Copyright Oswaldo Gressani 2025. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original data can be found directly in article:
# Tsang, Kenneth W., et al. "A cluster of cases of severe acute 
# respiratory syndrome in Hong Kong." New England Journal of 
# Medicine 348.20 (2003): 1977-1985.


library("xlsx")
library("EpiLPS")

rm(list = ls())

set.seed(1234)
# Data for 2003 SARS
n <- 9
tL <- c(2,2,6,2,1,5,5,1,2)
tR <- c(2,2,6,2,6,11,11,5,7)

for(i in 1:4){
  u1 <- runif(1)
  u2 <- runif(1)
  umin <- min(u1,u2)
  umax <- max(u1,u2)
  tL[i] <- tL[i] + umin
  tR[i] <- tR[i] + umax
}

for(i in 5:n){
  tL[i] <- tL[i] + runif(1)
  tR[i] <- tR[i] + runif(1)
}

data <- data.frame(tL = tL, tR = tR)
write.xlsx(round(data,3), file = "Data_continuous/Pathogen7-SARS2003.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 15)

# Extract estimates
Pathogen7_SARS2003_Estimates <- matrix(0, nrow = 2, ncol = 3)
colnames(Pathogen7_SARS2003_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen7_SARS2003_Estimates) <- c("Mean incubation period (days)",
                                                 "SD incubation period (days)")
Pathogen7_SARS2003_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)

# Write estimates in Estimates folder
write.xlsx(Pathogen7_SARS2003_Estimates, 
           file = "Estimates/Pathogen7_SARS2003_Estimates.xlsx",
           sheetName = "SARS2003-Semiparametric")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)












