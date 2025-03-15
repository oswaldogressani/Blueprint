#------------------------------------------------------------------------------#
#               EpiLPS incubation estimation Wuhan (Wild-type)                 #
#           Copyright Oswaldo Gressani 2025. All rights reserved.              #
#------------------------------------------------------------------------------#

library("EpiLPS")
library("xlsx")
library(tidyverse)
rm(list = ls())

#-- Code snippet adapted from Backer et al. 2020
dataraw <- read_tsv(file = "Data/Pathogen2-SARSCov2-Wild.tsv")

dataraw <- dataraw %>% 
  mutate(tReport = as.integer((`reporting date` %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tSymptomOnset = as.integer((symptom_onset %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tStartExposure = as.integer((exposure_start %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tEndExposure = as.integer((exposure_end %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))) %>%
  # if no start of exposure (i.e. Wuhan residents) use arbitrarily chosen start exposure in far away in past (here half December 2019)
  mutate(tStartExposure = ifelse(is.na(tStartExposure), min(tSymptomOnset)-8, tStartExposure)) %>%
  # if symptom onset in Wuhan, exposure ends at symptom onset
  mutate(tEndExposure = ifelse(tEndExposure >= tSymptomOnset, tSymptomOnset, tEndExposure))
#-- End of code snippet from Backer et al. 2020

tEL <- dataraw$tStartExposure[-which(is.na(dataraw$exposure_start))]
n <- length(tEL)
mintEL <- abs(min(tEL))

tSO <- dataraw$tSymptomOnset[-which(is.na(dataraw$exposure_start))] + mintEL 
tEL <- tEL + mintEL 
tER <- dataraw$tEndExposure[-which(is.na(dataraw$exposure_start))] + mintEL 

all(tSO >= tER) # Check if symptom onset is larger or equal to exposure end
all(tER >= tEL) # Check if exposure end is larger or equal to exposure start

# Making data continuous
datacts <- matrix(0, nrow = n, ncol = 6)
colnames(datacts) <- c("tEL", "tER", "Expowin", "tSO", "tIL", "tIR")
tL <- c()
tR <- c()
expowindow<- c()
set.seed(1234)
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
data <- data[-c(which(expowindow > 19)),]
datacts <- datacts[-c(which(expowindow > 19)),]
datacts <- as.data.frame(datacts)

# Check data constraints and export data
nsub <- nrow(datacts)
C1 <- sum(datacts$tEL >= 0) == nsub
C2 <- sum(datacts$tER > datacts$tEL) == nsub
C3 <- sum(datacts$tSO > datacts$tER) == nsub
if(C1 & C2 & C3){
  print("Data ok. All constraints satisfied.")
}
write.xlsx(round(datacts,3), file = "Data_continuous/Pathogen2-SARSCov2-Wild.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                       tmax = 20)

# Extract estimates
Pathogen2_Covid_Estimates <- matrix(0, nrow = 2, ncol = 3)
colnames(Pathogen2_Covid_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen2_Covid_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)")
Pathogen2_Covid_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)


# Write estimates in Estimates folder
write.xlsx(Pathogen2_Covid_Estimates, 
           file = "Estimates/Pathogen2_Covid_Estimates.xlsx",
           sheetName = "SARSCov2-Wild-LogNormal")

paste0("Sample size is n=", nrow(data))

# Additional summary statistics
round(incubfit$stats,1)



