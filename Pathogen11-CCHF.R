#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation CCHF                              #
#           Copyright Oswaldo Gressani 2025. All rights reserved.              #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("xlsx")
library("ggplot2")
rm(list = ls())

# Import raw data from xlsx file
dataraw <- read.xlsx(file = "Data/Pathogen11-CCHF.xlsx", sheetName = "Sheet1")
tEL <- dataraw$tEL[1:32]
tER <- dataraw$tER[1:32]
tSO <- dataraw$tS[1:32]

data <- data.frame(tEL = as.numeric(tEL), tER = as.numeric(tER), 
                   tSO = as.numeric(tSO))
data <- data[-c(which(tSO<tER),27),]
n <- nrow(data)

all(data$tSO >= data$tER) # Check if symptom onset is larger or equal to exposure end
all(data$tER >= data$tEL) # Check if exposure end is larger or equal to exposure start

# Making data continuous
datacts <- matrix(0, nrow = n, ncol = 6)
colnames(datacts) <- c("tEL", "tER", "Expowin", "tSO", "tIL", "tIR")
tL <- c()
tR <- c()
expowindow<- c()
set.seed(1234)
for(i in 1:n){
  tSO_temp <- data$tSO[i] + runif(1)
  tEL_temp <- data$tEL[i] + runif(1)
  tER_temp <- data$tER[i] + runif(1)
  while(tSO_temp <= tER_temp | tER_temp <= tEL_temp){
    tSO_temp <- data$tSO[i] + runif(1)
    tEL_temp <- data$tEL[i] + runif(1)
    tER_temp <- data$tER[i] + runif(1)
  }
 tL[i] <- tSO_temp - tER_temp
 tR[i] <- tSO_temp - tEL_temp
 expowindow[i] <- tER_temp - tEL_temp
 datacts[i, ] <- c(tEL_temp, tER_temp, expowindow[i], tSO_temp, tL[i], tR[i])
}

dataCCHF <- data.frame(tL = tL, tR = tR)
datacts <- as.data.frame(datacts)

# Check data constraints and export data
nsub <- nrow(datacts)
C1 <- sum(datacts$tEL >= 0) == nsub
C2 <- sum(datacts$tER > datacts$tEL) == nsub
C3 <- sum(datacts$tSO > datacts$tER) == nsub
if(C1 & C2 & C3){
  print("Data ok. All constraints satisfied.")
}
write.xlsx(round(datacts,3), file = "Data_continuous/Pathogen11-CCHF.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = dataCCHF, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 20)

# Extract estimates
Pathogen11_CCHF_Estimates <- matrix(0, nrow = 2, ncol = 3)
colnames(Pathogen11_CCHF_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen11_CCHF_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)")
Pathogen11_CCHF_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)

# Write estimates in Estimates folder
write.xlsx(Pathogen11_CCHF_Estimates, 
           file = "Estimates/Pathogen11_CCHF_Estimates.xlsx",
           sheetName = "CCHF-Gamma")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)





