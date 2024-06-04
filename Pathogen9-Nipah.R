#------------------------------------------------------------------------------#
#            EpiLPS incubation estimation of Nipah                             #
#           Copyright Oswaldo Gressani 2024. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original data can be found directly in supp material of article:
# Nikolay, B., Salje, H., Hossain, et al. (2019). Transmission of Nipah 
# virus—14 years of investigations in Bangladesh. New England Journal of 
# Medicine, 380(19), 1804-1814.


library("xlsx")
library("EpiLPS")
rm(list = ls())

set.seed(1234)

tSO <- c(12,13,14,17,13,14,13,13,17,15,13)
tEL<- c(6,6,6,6,4,5,5,5,5,1,1)
tER <- c(6,6,6,6,4,5,5,5,5,1,1)
n <- 11
tL <- c()
tR <- c()
data <- matrix(0, nrow = n, ncol = 2)
colnames(data) <- c("tL","tR")


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
  data[i, ] <- c(tL[i], tR[i])
}

data <- as.data.frame(data)
write.xlsx(round(data,3), file = "Data_continuous/Pathogen9-Nipah.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 18)

# Extract estimates
Pathogen9_Nipah_Estimates <- matrix(0, nrow = 3, ncol = 3)
colnames(Pathogen9_Nipah_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen9_Nipah_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)",
                                        "95% CI of incubation time (days)")
Pathogen9_Nipah_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)
Pathogen9_Nipah_Estimates[3,] <- c(NA,round(qgamma(p=c(0.025,0.975), 
                                                  shape = incubfit$shape,
                                                  rate = incubfit$rate), 1))

# Write estimates in Estimates folder
write.xlsx(Pathogen9_Nipah_Estimates, 
           file = "Estimates/Pathogen9_Nipah_Estimates.xlsx",
           sheetName = "Nipah-Gamma")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)













