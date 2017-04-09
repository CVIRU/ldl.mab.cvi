# Code: Meta-Analysis of LDL and CV Events Reduction Through PCSK9 Inhibition by mAb Tretments
# Created: 01/12/2017
# Author: Davit Sargsyan
#*************************************************************
# Header----
require(data.table)
require(metafor)
require(ggplot2)
require(truncnorm)
# Function to calculate truncated distribution mean----
MeanNormalTruncated <- function(mu = 0,
                                sigma = 1,
                                a = -Inf,
                                b = Inf){
  mu + sigma*(dnorm((a - mu)/sigma) - dnorm((b - mu)/sigma))/
    (pnorm((b - mu)/sigma) - pnorm((a - mu)/sigma))
}

# Mean LDL in Boco trials using trunkated mean----
drug <- "Bococizumab"
cutoff <- 150

# From Tables 1 and 2 (baseline * percent change):
# Active:
mu <- 63.6
# # Placebo:
# mu <- 114.6

std <- sqrt(27438)*((47.7 - 46.8)/2)/1.96
th <- 25

x <- seq(0, cutoff, 0.1)
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
                           b = cutoff)
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
              seq(th, cutoff, 0.1),
              cutoff),
        y = c(0, 
              dnorm(x = seq(th, cutoff, 0.1),
                    mean = mu,
                    sd = std),
              0),
        angle = -45,
        density = 10)
abline(v = c(mu1, mu2),
       lty = 2)
text(x = c(mu1, mu2),
     y = c(0.008, 0.008),
     labels = round(c(mu1, mu2)))

#*************************************************************
# Load data----
dt2 <- fread("docs/LDL cataract meta 3April2017.csv",
             sep = ",",
             header = TRUE)
dt2$LDL <- as.numeric(dt2$LDL)

# Model with square term
mm2 <- rma(ai = `Active Cataract`,
           bi = `Active N` - `Active Cataract`,
           ci = `Control Cataract`,
           di = `Control N` - `Control Cataract`, 
           measure = "OR",
           verbose = TRUE,
           method="DL",
           mods = as.matrix(data.table(x = LDL,
                                       x2 = LDL^2)),
           data = dt2)
mm2

# Generic forest plot
forest(mm2,
       slab = dt2$Study)

# Estimates
dt2$lor <- mm2$yi
dt2$lor.var <- mm2$vi
dt2$lor.ll95 <- mm2$yi - 1.96*sqrt(mm2$vi)
dt2$lor.ul95 <- mm2$yi + 1.96*sqrt(mm2$vi)

t1 <- data.table(Study = dt2$Study,
                 LDL = dt2$LDL,
                 `Log OR` = dt2$lor,
                 `Log OR 95% L.L.` = dt2$lor.ll95,
                 `Log OR 95% U.L.` = dt2$lor.ul95,
                 OR = exp(dt2$lor),
                 `OR 95% L.L.` = exp(dt2$lor.ll95),
                 `OR 95% U.L.` = exp(dt2$lor.ul95))
t1
write.csv(t1,
          file = "tmp/t1.csv",
          row.names = FALSE)

# Predicted values
mm2.pred <- predict(mm2)
mm2.pred <- data.table(Study = dt2$Study,
                       LDL = dt2$LDL,
                       `Predicted Log OR` = mm2.pred$pred,
                       `Predicted Log OR 95% L.L.` = mm2.pred$ci.lb,
                       `Predicted Log OR 95% U.L.` = mm2.pred$ci.ub)
mm2.pred
write.csv(mm2.pred,
          file = "tmp/t2.csv",
          row.names = FALSE)

# Polygon vertices
dt.poly <- data.table(id = rep(1:11, each = 4),
                      x = rep(c(0.8, 1, 1.2, 1), 11) + rep(0:10, each = 4),
                      y = c(t(data.table(mm2.pred$`Predicted Log OR`,
                                         mm2.pred$`Predicted Log OR 95% U.L.`,
                                         mm2.pred$`Predicted Log OR`,
                                         mm2.pred$`Predicted Log OR 95% L.L.`))))

# Plot estimates and predicted values
ggplot(dt2, 
       aes(x = Study, 
           y = lor)) + 
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = dt2$lor.ll95, 
                    ymax = dt2$lor.ul95),
                colour = "black", 
                width = .1) +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  geom_polygon(data = dt.poly,
               aes(x = x,
                   y = y,
                   group = id),
               alpha = 0.5) +
  scale_y_continuous("Log Odds Ratio") +
  ggtitle("Cataracts in LDL Lowering Studies") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))

# Fit a curve through the points
xx <- seq(min(dt2$LDL),
          max(dt2$LDL),
          0.01)

# With linear term only
# pdt <- predict(mm2,
#                newmods = xx)

# With quadratic term
pdt <- predict(mm2,
               newmods = as.matrix(data.table(x = xx,
                                              y = xx^2)))
pdt <- data.table(x = xx,
                  y = pdt$pred,
                  y.l = pdt$ci.lb,
                  y.u = pdt$ci.ub)

ggplot(dt2, 
       aes(x = LDL, 
           y = lor)) + 
  geom_point() +
  geom_point(data = mm2.pred,
             aes(x = LDL,
                 y = `Predicted Log OR`),
             size = 2,
             col = "red") +
  geom_line(data = pdt,
            aes(x = x,
                y = y),
            size = 1,
            col = "red") +
  geom_line(data = pdt,
            aes(x = x,
                y = y.l),
            size = 1,
            col = "grey",
            linetype = "dashed") +
  geom_line(data = pdt,
            aes(x = x,
                y = y),
            size = 1,
            col = "red") +
  geom_line(data = pdt,
            aes(x = x,
                y = y.u),
            size = 1,
            col = "grey",
            linetype = "dashed") +
  scale_y_continuous("Log Odds Ratio") +
  scale_x_continuous("Achieved LDL") +
  ggtitle("Linear and Quadratic Terms, with 95% C.I.") 