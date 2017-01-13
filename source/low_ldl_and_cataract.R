# Code: Meta-Analysis of LDL and CV Events Reduction Through PCSK9 Inhibition by mAb Tretments
# Created: 01/12/2017
# Author: Davit Sargsyan
#*****************************************************************************************
# Header----
require(data.table)
require(metafor)
require(ggplot2)

# Linear Regressions----
dt1 <- fread("docs/LDL cataract meta 3Jan2017.csv",
             sep = ",",
             header = TRUE)
dt1

m1 <- lm(dt1$`Odds Ratio` ~ dt1$`Achieved LDL`)
summary(m1)

dt1$ach_ldl_sqr <- (dt1$`Achieved LDL`)^2
m2 <- lm(dt1$`Odds Ratio` ~ dt1$`Achieved LDL` + dt1$ach_ldl_sqr)
summary(m2)

m3 <- lm(dt1$`Odds Ratio` ~ dt1$`Delta Pct`)
summary(m3)

dt1$delta_pct_sqr <- (dt1$`Achieved LDL`)^2
m4 <- lm(dt1$`Odds Ratio` ~ dt1$`Delta Pct` + dt1$delta_pct_sqr)
summary(m4)

# Plots----
par(mfrow = c(2, 2)) 
plot(dt1$`Odds Ratio` ~ dt1$`Achieved LDL`,
     col = as.numeric(factor(dt1$Trial)),
     pch = as.numeric(factor(dt1$Trial)))
abline(m1)
abline(m2, 
       col = "blue")

plot(dt1$`Odds Ratio` ~ dt1$`Delta Delta`,
     col = as.numeric(factor(dt1$Trial)),
     pch = as.numeric(factor(dt1$Trial)))

plot(dt1$`Odds Ratio` ~ dt1$`Delta Pct`,
     col = as.numeric(factor(dt1$Trial)),
     pch = as.numeric(factor(dt1$Trial)))
abline(m3)
abline(m4, 
       col = "blue")

plot.new()
legend("topleft",
       legend = dt1$Trial,
       col = as.numeric(factor(dt1$Trial)),
       pch = as.numeric(factor(dt1$Trial)))
graphics.off()

# Binning by percent change----
dt1$delta_pct.grp <- "<30%"
dt1$delta_pct.grp[dt1$`Delta Pct` >= 0.3 & dt1$`Delta Pct` < 0.5] <- ">=30% & <50%"
dt1$delta_pct.grp[dt1$`Delta Pct` >= 0.5] <- ">=50%"
dt1$delta_pct.grp <- factor(dt1$delta_pct.grp,
                            levels = c("<30%",
                                       ">=30% & <50%",
                                       ">=50%"))

plot(dt1$`Odds Ratio` ~ dt1$delta_pct.grp, type = "p")
mm1 <- lm(dt1$`Odds Ratio` ~ dt1$delta_pct.grp)
anova(mm1)
summary(mm1)

#*****************************************************************************************
# Meta-regressions----
dt2 <- fread("docs/LDL cataract meta 3Jan2017 Counts.csv",
             sep = ",",
             header = TRUE)
dt2 <- dt2[-3, ]
dt2$`BL-FU` <- as.numeric(dt2$`BL-FU`)
dt2$Study <- factor(dt2$Study,
                    levels = dt2$Study)
dt2

# # Log odds ratios of getting a cataract----
# mm2 <- rma(ai = dt2$`Active Cataract`,
#            bi = dt2$`Active N` - dt2$`Active Cataract`,
#            ci = dt2$`Control Cataract`,
#            di = dt2$`Control N` - dt2$`Control Cataract`, 
#            measure = "OR",
#            verbose = TRUE)
# mm2
# summary(mm2)
# forest(mm2,
#        slab = dt2$Study)

# Log odds ratios of getting a cataract, adjusted for percent delta achieved LDL---- 
mm2 <- rma(ai = `Active Cataract`,
           bi = `Active N` - `Active Cataract`,
           ci = `Control Cataract`,
           di = `Control N` - `Control Cataract`, 
           measure = "OR",
           verbose = TRUE,
           method="DL",
           mods = `BL-FU`,
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
                 `Active - Placebo BL - FU % Change` = dt2$`BL-FU`,
                 `Log OR` = dt2$lor,
                 `Log OR 95% L.L.` = dt2$lor.ll95,
                 `Log OR 95% U.L.` = dt2$lor.ul95,
                 OR = exp(dt2$lor),
                 `OR 95% L.L.` = exp(dt2$lor.ll95),
                 `OR 95% U.L.` = exp(dt2$lor.ul95))
write.csv(t1,
          file = "tmp/t1.csv",
          row.names = FALSE)

# Predicted values
mm2.pred <- predict(mm2)
mm2.pred <- data.table(Study = dt2$Study,
                       `Predicted Log OR` = mm2.pred$pred,
                       `Predicted Log OR 95% L.L.` = mm2.pred$ci.lb,
                       `Predicted Log OR 95% U.L.` = mm2.pred$ci.ub)
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