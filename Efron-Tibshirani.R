library(bootstrap)
library(boot)
library(lme4)
library(arm)

set.seed(20130810)

# CHAPTER 5, 6

law.sample <- data.frame(law)
law <- data.frame(law82)

summary(law.sample$GPA)

sd(law.sample$GPA)/sqrt(length(law.sample$GPA))

summary(law$GPA)

# Bootstrap estimation of correlation, replicates table 6.1

cor(law$LSAT, law$GPA) # 0.759
cor(law.sample$LSAT, law.sample$GPA) # 0.776

boot.corr <- function(data, i) {
  cor(data[i, ]$LSAT, data[i, ]$GPA)
}

for (r in c(25, 50, 100, 200, 400, 800, 1600, 3200)) {
  print(law.cor.boot <- boot(data = law.sample, 
                             statistic = boot.corr,
                             R = r, 
                             parallel = "multicore"))
}

# Store the bootstrap replicates of the correlation for 3200 samples; replicates fig 6.2

law.cor.boot <- boot(data = law.sample,
                     statistic = boot.corr,
                     R = 3200,
                     parallel = "multicore")
names(law.cor.boot)
hist(law.cor.boot$t, 
     xlab = "Bootstrap estimate of correlation", 
     ylab = "Count", 
     main ="Histogram of 3200 replicates")

# We can find the mean and variance of the standard error itself to calculate the coefficient of variation

(coef.variation <- sd(law.cor.boot$t)/mean(law.cor.boot$t))

# Problem 6.10

d <- c(1, 2, 3.5, 4, 7, 7, 3, 8.6, 12.4, 13.8, 18.1)
(trimmed.mean <- mean(d, trim = 0.25))

trimmed.mean <- function(data, i) {
  mean(data[i], trim = 0.25)
}

set.seed(20130810)
(boot.trimmed.mean <- boot(data = d, 
                           statistic = trimmed.mean,
                           R = 3200,
                           parallel = "multicore"))
set.seed(1)
(boot.trimmed.mean <- boot(data = d, 
                           statistic = trimmed.mean,
                           R = 3200,
                           parallel = "multicore"))

# CHAPTER 7

# Test scores example

test.scores <- data.frame(scor)
(mean.vector.test <- colMeans(test.scores))
(test.var.cov.matrix <- cov(test.scores))
(test.eig.values <- eigen(test.var.cov.matrix))
(theta.hat <- test.eig.values$values[1]/sum(test.eig.values$values)) # = 0.619, variance explained by the first component

eig.stat <- function(data, i) {
  eigen(cov(data[i, ]))$values[1]/sum(eigen(cov(data[i, ]))$values)
}

(eig.stat.boot <- boot(data = test.scores,
                       statistic = eig.stat,
                       R = 200,
                       parallel = "multicore"))
mean(eig.stat.boot$t) # = 0.624
hist(eig.stat.boot$t, 
     xlab = expression(theta), 
     ylab = "Count", 
     main = "200 bootstrap replicates of eigen value statistic")
abline(v = 0.619, lty = "dotted")

(boot.ci(eig.stat.boot))

# Cholostyramne example

cholost.data <- data.frame(cholost)
summary(cholost.data)
with(cholost.data, plot(z, y, pch = "+", xlab = "compliance", ylab = "improvement")); abline(h = 0) # replicates figure 7.5

(quad.model <- lm(y ~ poly(z, 2), data = cholost.data))
(loess.model <- loess(y ~ z, data = cholost.data, span = 0.3))

# Generate some predicted data to enable plotting the quadratic

quad.data <- data.frame(z = seq(min(cholost.data$z), max(cholost.data$z), by = 0.01))
quad.data$y.hat <- predict(quad.model, quad.data)

with (cholost.data, 
      scatter.smooth(z, y, span = 0.3, 
                     xlab = "compliance", ylab = "improvement"))
lines(x = quad.data$z, y = quad.data$y.hat, type = "l", lty = "dotted") # replicates figure 7.6

# For predicted values at z = 60 and 100, table 7.5, bootstrapped standard errors in the predicted values calculated at 50 replicates

(predict(quad.model, data.frame(z = c(60, 100))))
(predict(loess.model, data.frame(z = c(60, 100))))

pred.z.60.quad <- function(d, i) {
  predict(lm(y ~ poly(z, 2), data = d[i, ]), data.frame(z = 60))
}

pred.z.100.quad <- function(d, i) {
  predict(lm(y ~ poly(z, 2), data = d[i, ]), data.frame(z = 100))
}

pred.z.60.loess <- function(d, i) {
  predict(loess(y ~ z, data = d[i, ], span = 0.3), data.frame(z = 60))
}

pred.z.100.loess <- function(d, i) {
  predict(loess(y ~ z, data = d[i, ], span = 0.3), data.frame(z = 100))
}

(boot.z.60.quad <- boot(data = cholost.data,
                        statistic = pred.z.60.quad,
                        R = 50, 
                        parallel = "multicore"))

(boot.z.100.quad <- boot(data = cholost.data,
                         statistic = pred.z.100.quad,
                         R = 50,
                         parallel = "multicore"))
(boot.z.60.loess <- boot(data = cholost.data,
                         statistic = pred.z.60.loess, 
                         R =50,
                         parallel = "multicore"))

(boot.z.100.loess <- boot(data = cholost.data,
                          statistic = pred.z.100.loess,
                          R = 50,
                          parallel = "multicore"))

# Problem 7.10

# randomly select 30% of the data and fit a quadratic model

(sample.indices <- sample(1:nrow(cholost.data), 0.3 * nrow(cholost.data), replace = FALSE))
(sample.cholost.data <- cholost.data[sample.indices, ])

(sample.quad.model <- lm(y ~ poly(z, 2), data = sample.cholost.data))

(boot.z.60.quad.sample <- boot(data = sample.cholost.data,
                               statistic = pred.z.60.quad,
                               R = 50, 
                               parallel = "multicore"))

# CHAPTER 8

# Two sample problem

mouse.c <- data.frame(mouse.c)
mouse.t <- data.frame(mouse.t)
n.samples <- 1400

mouse.c.boot <- matrix(nrow = nrow(mouse.c), ncol = n.samples)
mouse.t.boot <- matrix(nrow = nrow(mouse.t), ncol = n.samples)

for (i in 1:n.samples) {
  mouse.c.boot[, i] <- sample(x = mouse.c$mouse.c, size = nrow(mouse.c), replace = TRUE)
  mouse.t.boot[, i] <- sample(x = mouse.t$mouse.t, size = nrow(mouse.t), replace = TRUE)
}

(mouse.theta.hat <- colMeans(mouse.t.boot) - colMeans(mouse.c.boot))
hist(mouse.theta.hat, xlab = expression(theta), ylab = "Count", main = "Bootstrap replicates of differences"); 
abline(v = mean(mouse.t$mouse.t) - mean(mouse.c$mouse.c), lty = "dotted") # replicates figure 8.2
(mouse.theta.se <- sd(mouse.theta.hat)) # = 25.92

# Time series data

lutenhorm <- data.frame(lutenhorm)
lutenhorm$t <- 1:48
colnames(lutenhorm) <- c("level", "t")
with(lutenhorm, plot(t, level, type = "l")); abline(h = mean(lutenhorm$level), lty = "dotted") # replicates figure 8.4

lutenhorm$z <- lutenhorm$level - mean(lutenhorm$level) # centers the error terms to a 0 mean
(luten.ar1.fit <- ar(lutenhorm$z, order.max = 1)) # beta.hat = 0.5755

(luten.beta.hat <- luten.ar1.fit$ar)

luten.fit.errors <- rep(NA, nrow(lutenhorm)-1)

for (t in 1:nrow(lutenhorm)-1) {
  luten.fit.errors[t] <- lutenhorm$z[t+1] - luten.beta.hat * lutenhorm$z[t] 
}

hist(luten.fit.errors, xlab = " ", ylab = " ", main = " "); abline(v = mean(luten.fit.errors), lty = "dotted") # replicates figure 8.5

# Generate the bootstrapped time series

n.samples <- 200

luten.series.boot <- matrix(nrow = nrow(lutenhorm), ncol = n.samples)

for (i in 1:n.samples) {
  luten.series.boot[1, i] <- lutenhorm$z[1] 
  
  for (t in 2:nrow(luten.series.boot)) {
    luten.series.boot[t, i] <- luten.series.boot[t-1, i] * luten.beta.hat + sample(x = luten.fit.errors, size = 1, replace = TRUE)
  }
}

# fit the AR(1) model to each of the time series

luten.boot.betas <- apply(luten.series.boot, 2, ar, order.max = 1) # fit the model
luten.boot.betas <- unlist(lapply(luten.boot.betas, FUN = function(x) {x$ar})) # extract the coefficients and unlist them
hist(luten.boot.betas, xlab = " ", ylab = " ", main = " "); abline(v = luten.beta.hat, lty = "dotted") # replicates figure 8.6

# CHAPTER 9

# hormone data

hormone <- data.frame(hormone)

summary(hormone.fit1 <- lm(amount ~ hrs, data = hormone))

with(hormone, plot(hrs, amount, pch = Lot)); abline(hormone.fit1) # reproduces figure 9.1

summary(hormone.fit2 <- lmer(amount ~ (1|Lot) + hrs, data = hormone)) # varying intercept, constant slope model
(coef(hormone.fit2)$Lot)
se.coef(hormone.fit2)

summary(hormone.fit3 <- lm(amount ~ hrs + Lot, data = hormone)) # Table 9.3

# cell survival data

cell.survival.data <- data.frame(cell)

summary(cell.fit1 <- lm(log.surv ~ 0 + dose, data = cell.survival.data))
summary(cell.fit2 <- lm(log.surv ~ 0 + dose + I(dose^2), data = cell.survival.data))  # For all 14 plates
summary(cell.fit3 <- lm(log.surv ~ 0 + dose + I(dose^2), data = cell.survival.data[-13, ])) # Exclude plate 13, replicate Table 9.5

(beta.2.hat <- c(coef(cell.fit2)[2], coef(cell.fit3)[2]))
(beta.2.se <- c(se.coef(cell.fit2)[2], se.coef(cell.fit3)[2]))
(beta.2.ratio <- beta.2.hat/beta.2.se) # Removing one data point changes the conclusion of the study

# Problem 9.8

(cell.fit4 <- MASS::lmsreg(log.surv ~ 0 + dose + I(dose^2), data = cell.survival.data)) # least median squares estimate

# CHAPTER 10

patch.data <- data.frame(patch)

(theta.hat <- mean(patch.data$y)/mean(patch.data$z))

theta.boot <- function(data, i) {
  mean(data[i, ]$y)/mean(data[i, ]$z)
}

boot.theta.patch <- boot(data = patch.data,
                         statistic = theta.boot,
                         R = 400,
                         parallel = "multicore")

hist(boot.theta.patch$t, xlab = "Ratio statistic"); abline(v = theta.hat, lty = "dotted") # Replicates figure 10.1

# CHAPTER 12


mouse.data <- data.frame(c = mouse.c, t = c(mouse.t, rep(NA, 2)))
(theta.hat <- mean(mouse.data$mouse.c))

boot.mean <- function(data, index) {
  mean(data[index, ]$mouse.c)
}

(theta.se <- boot(data = mouse.data,
                  statistic = boot.mean,
                  R = 10000,
                  parallel = "multicore"))

# Problem 12.5
(sample <- rchisq(n = 20, df = 1))
(sample.mean <- mean(sample))
(sample.se <- sd(sample)/sqrt(length(sample)))

# part a
int.a <- qchisq(p = c())





# CHAPTER 13

mouse.t <- data.frame(t = bootstrap::mouse.t)

theta.hat <- mean(mouse.t$t)

boot.mean <- function(data, indices) {
  mean(data[indices, ])
}

(theta.boot <- boot(data = mouse.t, statistic = boot.mean, R = 10000, parallel = "multicore"))

(theta.90ci <- theta.hat + sd(theta.boot$t) * qnorm(p = c(0.05, 0.95))) # 90% CI for the mean of 7 treated mice
hist(theta.boot$t); abline(v = theta.hat); abline(v = theta.90ci, lty = "dotted") # Replicates Figure 13.1

(theta.hat.star.quantiles <- quantile(theta.boot$t, probs = c(0.025, 0.05, 0.1, 0.16, 0.5, 0.84, 0.90, 0.95, 0.975))) # replicates Table 13.1

# Replication of discussion around Figure 2

(x <- rnorm(n = 10))
(theta <- exp(x))

theta.hat <- mean(theta)

boot.mean <- function(data, indices) {
  mean(data[indices])
}

theta.boot <- boot(data = theta, statistic = boot.mean, R = 1000, parallel = "multicore")
hist(theta.boot$t); abline(v = theta.hat + sd(theta.boot$t) * qnorm(p = c(0.05, 0.95)), lty = "dotted")

# Problem 13.4

(x <- rexp(n = 20, rate = 1))
(theta.hat <- mean(x))
theta.true <- 1
n.intervals <- 1000

boot.mean <- function(data, indices) {
  mean(data[indices])
}

theta.ci <- matrix(nrow = n.intervals, ncol = 2)

for (i in 1:n.intervals) {
  theta.boot <- boot(data = x, statistic = boot.mean, R = 1000, parallel = "multicore")
  theta.ci[i, ] <- quantile(x = theta.boot$t, p = c(0.025, 0.975))
}

sum(theta.ci[, 1] > theta.true)
sum(theta.ci[, 2] < theta.true)

# CHAPTER 14

# Spatial Data

spatial.data <- data.frame(spatial)

plot(spatial.data$A, spatial.data$B) # replicates Figure 14.1

(theta.hat <- var(spatial.data$A)) # unbiased estimate, n-1 in the denominator

boot.var <- function(data, indices) {
  var(data[indices, ]$A)
}

theta.boot <- boot(data = spatial.data,
                   statistic = boot.var,
                   R = 2000,
                   parallel = "multicore")
hist(theta.boot$t); abline(v = theta.hat) # replicates left half of Figure 14.3

quantile(x = theta.boot$t, probs = c(0.025, 0.975)) # according to the percentile method from the previous chapter

# the boot object gets these directly

(theta.ci <- boot.ci(boot.out = theta.boot, conf = 0.95, type = c("norm", "perc", "bca")))

bcanon(spatial.data$A, nboot = 2000, theta = var)

# Tooth data

tooth.data <- data.frame(tooth)

summary(tooth.fit1 <- lm(strength ~ D1 + D2, data = tooth.data))
summary(tooth.fit2 <- lm(strength ~ E1 + E2, data = tooth.data))

sd.diff.hat <- sd(tooth.fit2$residuals) - sd(tooth.fit1$residuals)

sd.diff <- function(data, indices) {
  f1 <- lm(strength ~ D1 + D2, data = data[indices, ])
  f2 <- lm(strength ~ E1 + E2, data = data[indices, ])
  
  return(sd(f2$residuals) - sd(f1$residuals))
}

(sd.diff.boot <- boot(data = tooth.data, statistic = sd.diff, R = 2000, parallel = "multicore"))
hist(sd.diff.boot$t); abline(v = sd.diff.hat) # Replicates Figure 14.5

(sd.diff.ci <- boot.ci(boot.out = sd.diff.boot, conf = .95, type = c("perc", "bca"))) # Cant reject the null that the residual variance is different 


