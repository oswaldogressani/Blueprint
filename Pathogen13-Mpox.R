#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation Mpox                              #
#           Copyright Oswaldo Gressani 2024. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original xls data downloaded from Miura et al. (2022)
# https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.24.2200448

library("xlsx")
library("EpiLPS")

## Load Data
data <- read.csv("Data/Pathogen13-Mpox.csv")
n <- nrow(data)

tEL <- data$Start.date.of.exposure
tER <- data$End.date.of.exposure
tSO <- data$Symptom.onset

# Making data continuous
datacts <- matrix(0, nrow = n, ncol = 6)
colnames(datacts) <- c("tEL", "tER", "Expowin", "tSO", "tIL", "tIR")
tL <- c()
tR <- c()
expowindow<- c()
set.seed(123)
for(i in 1:n){
  tSO_temp <- tSO[i] + runif(1)
  tEL_temp <- tEL[i] + runif(1)
  tER_temp <- tER[i] + runif(1)
  while(tSO_temp <= tER_temp | tER_temp <= tEL_temp){
    tSO_temp <- tSO[i] + runif(1)
    tEL_temp <- tEL[i] + runif(1)
    tER_temp <- tER[i] + runif(1)
  }
  tL[i] <- tSO_temp - tER_temp
  tR[i] <- tSO_temp - tEL_temp
  expowindow[i] <- tER_temp - tEL_temp
  datacts[i, ] <- c(tEL_temp, tER_temp, expowindow[i], tSO_temp, tL[i], tR[i])
}

data <- data.frame(tL = tL, tR = tR)
datacts <- as.data.frame(datacts)

# Check data constraints and export data
nsub <- nrow(datacts)
C1 <- sum(datacts$tEL >= 0) == nsub
C2 <- sum(datacts$tER > datacts$tEL) == nsub
C3 <- sum(datacts$tSO > datacts$tER) == nsub
if(C1 & C2 & C3){
  print("Data ok. All constraints satisfied.")
}
write.xlsx(round(datacts,3), file = "Data_continuous/Pathogen13-Mpox.xls")


# Fit with EpiLPS

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 30)

# Extract estimates
Pathogen13_Mpox_Estimates <- matrix(0, nrow = 3, ncol = 3)
colnames(Pathogen13_Mpox_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen13_Mpox_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)",
                                        "95% CI of incubation time (days)")
Pathogen13_Mpox_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)
Pathogen13_Mpox_Estimates[3,] <- c(NA,round(qlnorm(p=c(0.025,0.975), 
                                                   meanlog = incubfit$meanlog, 
                                                   sdlog = incubfit$sdlog), 1))

# Write estimates in Estimates folder
write.xlsx(Pathogen13_Mpox_Estimates, 
           file = "Estimates/Pathogen13_Mpox_Estimates.xlsx",
           sheetName = "Mpox-LogNormal")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)































