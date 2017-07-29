# |-------------------------------------------------------|
# | Project: A Study of Association of LDL Reduction and  | 
# |          Cataracts Using Meta-Analysis                |
# | Script: Mean of Truncated Normal Distribution                 |
# | Author: Davit Sargsyan                                |
# | Created: 11/23/2016                                   |
# |-------------------------------------------------------|
# Source:
# https://assessingpsyche.wordpress.com/2014/06/04/using-the-truncated-normal-distribution/

# Data----
# Alirocumab: mean LDL level (mg/dL) at Baseline and Week 24
mu.ali0 <- 122.8
sem.ali0 <- 42.7
N.ali0 <- 1530

mu.ali24 <- 48.3 
sem.ali24 <- 0.9
N.ali24 <- 1386

# Evolocumab: median LDL level (mg/dL) at Baseline and Week 12
# Based on Figure 1 table
mu.evo0 <- 100*73.4/60.9
mu.evo12 <- mu.evo0 - 73.4

# Based on Supplemental Table S1 (Appendix)
N.evo0 <- 2976
N.evo12 <- 2871
# NOTE: VERY BAD ESTIMATE! 
# Assuming 95% C.I. of achieved LDL level = 4*SD 
std.evo0 <- (150 - 97)/4
std.evo12 <- (71 - 32)/4

#*****************************************************************************************
# Function to calculate truncated distribution mean----
MeanNormalTruncated <- function(mu = 0,
                                sigma = 1,
                                a = -Inf,
                                b = Inf){
  mu + sigma*(dnorm((a - mu)/sigma) - dnorm((b - mu)/sigma))/
    (pnorm((b - mu)/sigma) - pnorm((a - mu)/sigma))
}

#*****************************************************************************************
# Calculate and plot expected values for given threshods----
# drug <- "Alirocumab"
# mu <- mu.ali24
# std <- sem.ali24*sqrt(N.ali24)
# N <- N.ali24

drug <- "Evolocumab"
mu <- mu.evo12
std <- std.evo12
N <- N.evo12

# th <- 15
# th <- 25
th <- 40 

x <- seq(0, 100, 0.1)
y <- dnorm(x = x, 
           mean = mu, 
           sd = std)

mu1 <- MeanNormalTruncated(mu = mu,
                           sigma = std,
                           a = 0,
                           b = th)
mu1

mu2 <- MeanNormalTruncated(mu = mu,
                           sigma = std,
                           a = th,
                           b = 100)
mu2

plot(y ~ x,
     type = "l",
     xlab = "LDL",
     ylab = "Probability",
     main = paste("Simulation of",
                  drug,
                  "Effect on LDL \n Mean =" ,
                  round(mu, 2),
                  ", SD = ",
                  round(std, 2),
                  ", N = ",
                  N,
                  ", Threshold =",
                  th))
polygon(x = c(0, 
              seq(0, th, 0.1),
              th),
        y = c(0, 
              dnorm(x = seq(0, th, 0.1),
                    mean = mu,
                    sd = std),
              0),
        angle = 45,
        density = 10)
polygon(x = c(th, 
              seq(th, 100, 0.1),
              100),
        y = c(0, 
              dnorm(x = seq(th, 100, 0.1),
                    mean = mu,
                    sd = std),
              0),
        angle = -45,
        density = 10)
abline(v = c(mu1, mu2),
       lty = 2)
text(x = c(mu1, mu2),
     y = c(0.008, 0.012),
     labels = round(c(mu1, mu2)))