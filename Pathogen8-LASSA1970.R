#------------------------------------------------------------------------------#
#            EpiLPS incubation estimation of LASSA 1970 Nigeria                #
#           Copyright Oswaldo Gressani 2024. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original data can be found directly in supp material of article:
# Akhmetzhanov Andrei R., Asai Yusuke and Nishiura Hiroshi 2019. Quantifying 
# the seasonal drivers of transmission for Lassa fever in Nigeria. 
# Phil. Trans. R. Soc. B37420180268
# https://doi.org/10.1098/rstb.2018.0268


library("xlsx")
library("EpiLPS")
rm(list = ls())

set.seed(1234)

tSO <- c(17, 19, 21, 21, 21, 22, 22, 23, 23, 23, 
        23, 23, 25, 26, 26, 26, 27, 31, 39, 41, 
        41, 42, 44)
tEL<- c(5,  9,  5,  9,  5,  5,  5,  11, 5, 5,
        13, 11, 12, 5,  11, 11, 5,  5,  12, 25,
        5,  25, 25)
tER <- c(15, 18, 18, 15, 18, 18, 18, 15, 18, 18,
            14, 15, 18, 18, 15, 15, 18, 18, 18, 31,
            18, 31, 31)
n <- 23
tL <- c()
tR <- c()
data<- matrix(0, nrow = n, ncol = 2)
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
write.xlsx(round(data,3), file = "Data_continuous/Pathogen8-LASSA1970.xlsx")

# Fit incubation density with EpiLPS
incubfit <- EpiLPS::estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 42)

# Extract estimates
Pathogen8_LASSA_Estimates <- matrix(0, nrow = 3, ncol = 3)
colnames(Pathogen8_LASSA_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen8_LASSA_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)",
                                        "95% CI of incubation time (days)")
Pathogen8_LASSA_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)
Pathogen8_LASSA_Estimates[3,] <- c(NA,round(qlnorm(p=c(0.025,0.975), 
                                                  meanlog = incubfit$meanlog, 
                                                  sdlog = incubfit$sdlog), 1))

# Write estimates in Estimates folder
write.xlsx(Pathogen8_LASSA_Estimates, 
           file = "Estimates/Pathogen8_LASSA_Estimates.xlsx",
           sheetName = "LASSA-LogNormal")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)













