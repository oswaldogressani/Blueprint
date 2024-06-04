#------------------------------------------------------------------------------#
#            EpiLPS incubation estimation of Marburg                           #
#           Copyright Oswaldo Gressani 2024. All rights reserved.              #
#------------------------------------------------------------------------------#

# Original data can be found directly in supp material of article:
# Pavlin, B.I. Calculation of incubation period and serial interval 
# from multiple outbreaks of Marburg virus disease. 
# BMC Res Notes 7, 906 (2014). https://doi.org/10.1186/1756-0500-7-906


library("xlsx")
library("EpiLPS")
rm(list = ls())

set.seed(1234)

texact <- c(2,2,3,3,4,5,5,7,7,7,7,7,8,9,9,9,9,10,13,13)
n <- length(texact)

tL <- texact - runif(n)
tR <- texact + runif(n)

data <- matrix(0, nrow = n, ncol = 2)
colnames(data) <- c("tL","tR")
data <- as.data.frame(data)
data$tL <- tL
data$tR <- tR
write.xlsx(round(data,3), file = "Data_continuous/Pathogen10-Marburg.xlsx")

# Fit incubation density with EpiLPS
incubfit <- estimIncub(x = data, K = 20, niter = 20000, verbose = TRUE,
                               tmax = 17)

# Extract estimates
Pathogen10_Marburg_Estimates <- matrix(0, nrow = 3, ncol = 3)
colnames(Pathogen10_Marburg_Estimates) <- c("Point estimate", "CI95L", "CI95R")
rownames(Pathogen10_Marburg_Estimates) <- c("Mean incubation period (days)",
                                        "SD incubation period (days)",
                                        "95% CI of incubation time (days)")
Pathogen10_Marburg_Estimates[c(1:2),] <- round(incubfit$stats[c(1,2),c(1,4,5)],1)
Pathogen10_Marburg_Estimates[3,] <- c(NA,round(qweibull(p=c(0.025,0.975), 
                                                  shape = incubfit$shape,
                                                  scale = incubfit$scale), 1))

# Write estimates in Estimates folder
write.xlsx(Pathogen10_Marburg_Estimates, 
           file = "Estimates/Pathogen10_Marburg_Estimates.xlsx",
           sheetName = "Marburg-Weibull")

paste0("Sample size is n=", n)

# Additional summary statistics
round(incubfit$stats,1)















