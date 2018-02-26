#' ---
#' title: "Introduction to bootstrap with R "
#' author: Pavan Gurazada
#' output: github_document
#' ---
#' 

library(tidyverse)
library(boot)
library(caret)

set.seed(20130810)

#' Example 2.1
#' Bootstrapping variance of univariate gaussians

n <- 25
x <- rnorm(n) # true value of the variance is 1

#' The statistic we want to explore is the uncorrected variance
varx <- var(x) * (n-1)/n # Uncorrected

nboot <- 5000 # Number of bootstrap samples
bootVarx <- rep(NA, nboot)

for (i in 1:nboot) {
  xStar <- sample(x, size = n, replace = TRUE)
  bootVarx[i] <- var(xStar) * (n - 1)/n
}

cat("Mean of the variance: ", mean(bootVarx))
qplot(bootVarx, geom = "histogram") +
  labs(x = "Variance of the bootstrapped samples")

#' Example of classification

classes <- data.frame(class = c(rep(1, 5), rep(2, 5)),
                      x1 = c(2.052,1.083,0.083,1.278, -1.226, 1.307, -0.548, 
                             2.498,0.832,1.498),
                      x2 = c(0.339, -1.320, -1.524, -0.459, -0.606, 2.268, 
                             1.741, 0.813, 1.409, 2.063))
glimpse(classes)

ldaModel <- train(factor(class) ~ .,
                  data = classes,
                  method = "lda",
                  trControl = trainControl(method = "boot632"))
ldaModel$results

#' Exercise 2.1

x <- c(23, 16, 21, 24, 34, 30, 28, 24, 26, 18, 23, 23, 36,
        37, 49, 50, 51, 56, 46, 41, 54, 30, 40, 31)
n <- length(x)

mean(x); sd(x)

nboot <- 1000
bootMean <- rep(NA, nboot)
bootSD <- rep(NA, nboot)

for (i in 1:nboot) {
  xStar <- sample(x, size = n, replace = TRUE)
  bootMean[i] <- mean(xStar)
  bootSD[i] <- sd(xStar)
}

qplot(bootMean, geom = "histogram") +
  labs(x = "Mean of bootstrappped sample")
qplot(bootSD, geom = "histogram") +
  labs(x = "SD of bootstrapped samples")

#' *Example 3.1* “ Surimi ” is purified fish protein used as a material to make
#' imitation crab and shrimp food products. The strength of surimi gels is a
#' critical factor in production. Each incoming lot of surimi raw material is
#' sampled and a cooked gel is prepared. From these gels, test portions are
#' selected and tested for strength.
#'
#' We suppose our interest here is the mean value of the ultimate stresses and
#' the 95% confidence interval for it.

x <- c(41.28, 45.16, 34.75, 40.76, 43.61, 39.05, 41.20, 41.02, 41.33, 40.61, 
        40.49, 41.77, 42.07, 44.83, 29.12, 45.59, 41.95, 45.78, 42.89, 40.42, 
        49.31, 44.01, 34.87, 38.60, 39.63, 38.52, 38.52, 43.95, 49.08, 50.52, 
        43.85, 40.64, 45.86, 41.25, 50.35, 45.18, 39.67, 43.89, 43.89, 42.16)
summary(x)
var(x); sd(x)
n <- length(x)

qplot(x, geom = "histogram") +
  labs(x = "Data")

nboot <- 1000
bootMean <- rep(NA, nboot)

for (i in 1:nboot) {
  xStar <- sample(x, size = n, replace = TRUE)
  bootMean[i] = mean(xStar)
}

qplot(bootMean, geom = "histogram") +
  labs(x = "Mean of bootstrapped samples")

bootStat <- function(data, i) {
  return(mean(data[i]))
}

bootMean <- boot(data = x, 
                 statistic = bootStat,
                 R = 1000)
boot.ci(bootMean, type = 'bca', conf = 0.95)

#' *Exercise 3.6* Microbiological method comparison: A laboratory wants to
#' determine if two different methods give similar results for quantifying a
#' particular bacterial species in a particular medium. A quick test is
#' performed by a decimal serial dilution assay, using each method in replicate
#' on each dilution. The dilution with the highest usable counts (falling the
#' range 30 – 300) had results (actual data):

results <- data.frame(A = c(176, 125, 152, 180, 159, 168, 160, 151),
                      B = c(164, 121, 137, 169, 144, 145, 156, 139))

#' A and B are the two different methods used
#'
#' The usual procedure in a parametric case to test the hypothesis that the two
#' methods A and B are different is to assume that the differences in the mean
#' dilation form the two samples follows a normal distribution of mean 0 and
#' check the p value. Lets look at a bootstrap method. The bootstrap sample
#' should follow the probability pattern of the data generating process.
#'
#' In this simple case, we resample from the two methods separately and
#' aggregate the difference.
#' 
#' This method scales to AB testing as well.

nboot <- 1000

bootfn <- function(data, i) {
  mA <- mean(data[i, 'A'])
  mB <- mean(data[i, 'B'])
  return(mA-mB)
}

bootAB <- boot(data = results,
               statistic = bootfn,
               R = nboot)
boot.ci(bootAB, conf = 0.95, type = "bca")

#' Since the interval does not include 0, we conclude that the A and B are
#' different
#' 
#' *Alternative implementation* without using the boot package

bootAB <- rep(NA, nboot)

for (i in 1:nboot) {
  mA <- mean(sample(results$A, size = length(results$A), replace = TRUE))
  mB <- mean(sample(results$B, size = length(results$B), replace = TRUE))
  bootAB[i] <- mA - mB
}

qplot(bootAB, geom = "histogram") +
  labs(x = "Difference in bootstrapped samples")


