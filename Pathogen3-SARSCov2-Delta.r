#------------------------------------------------------------------------------#
#         EpiLPS incubation estimation SARS-CoV-2 Delta variant                #
#             Oswaldo Gressani 2024. All rights reserved.                      #
#------------------------------------------------------------------------------#

library("tidyverse")
library("lubridate")
library("xlsx")
library("EpiLPS")
rm(list = ls())

# Code by Jantien Backer--------------------------------------
data_SGTF_incubation <- read_csv("Data/Pathogen3-SARSCov2-Delta.csv")

data <- data_SGTF_incubation %>% 
  filter(type == "non-SGTF") %>%   # SGTF = Omicron
  mutate(tReport = as.integer(infectee_report - as.Date("2021-12-01")),
         tSymptomOnset = as.integer(infectee_SO - as.Date("2021-12-01")),
         tStartExposure = as.integer(exposure_start - as.Date("2021-12-01")) - 0.5,
         tEndExposure = as.integer(exposure_end - as.Date("2021-12-01")) + 0.5) %>%
  mutate(tStartExposure = ifelse(is.na(tStartExposure), min(tSymptomOnset)-21, tStartExposure)) %>%
  mutate(tEndExposure = ifelse(tEndExposure > tSymptomOnset, tSymptomOnset, tEndExposure)) 

input_data <- list(
  N = nrow(data),
  tStartExposure = data$tStartExposure,
  tEndExposure = data$tEndExposure,
  tSymptomOnset = data$tSymptomOnset)
# Code by Jantien Backer ends here --------------------------------------

set.seed(1234)

dataraw <- data.frame(tEL = input_data$tStartExposure,
                      tER = input_data$tEndExposure, 
                      tSO = input_data$tSymptomOnset)
dataraw <- dataraw[-which(dataraw$tEL>dataraw$tER),]

offset <- abs(min(dataraw$tEL)-1)
tSO <- dataraw$tSO+offset
tEL <- dataraw$tEL+offset
tER <- dataraw$tER+offset

n <- nrow(dataraw)
tL <- c()
tR <- c()
dataset <- matrix(0, nrow = n, ncol = 2)
colnames(dataset) <- c("tL","tR")


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
  dataset[i, ] <- c(tL[i], tR[i])
}

dataset <- as.data.frame(dataset)
dataset <- dataset[-which(dataset$tR>20),]
n <- nrow(dataset)
write.xlsx(round(dataset,3), file = "Data_continuous/Pathogen3-SARSCov2-Delta.xlsx")


# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = dataset, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 20)

# Extract estimates
Pathogen3_Covid_Estimates <- matrix(0, nrow = 3, ncol = 3)
colnames(Pathogen3_Covid_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen3_Covid_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)",
                                        "95% CI of incubation time (days)")
Pathogen3_Covid_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)
Pathogen3_Covid_Estimates[3,] <- c(NA,round(qgamma(p=c(0.025,0.975),
                                                   shape = incubfit$shape,
                                                   rate = incubfit$rate), 1))

# Write estimates in Estimates folder
write.xlsx(Pathogen3_Covid_Estimates, 
           file = "Estimates/Pathogen3_Covid_Estimates.xlsx",
           sheetName = "SARSCov2-Delta_variant-Gamma")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)


