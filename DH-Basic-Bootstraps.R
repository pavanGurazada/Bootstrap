#' ---
#' title: "Chapter 2: The Basic Bootstraps"
#' author: Davison and Hinkley
#' output: github_document
#' ---

#' last update: Wed Apr 18 10:01:45 2018

library(boot)
library(tidyverse)
library(ggthemes)

set.seed(20130810)

theme_set(theme_few())

#' Here, we are ocncerned about bootstrap estimates from a single, homonegenous
#' set of data
#'
#' The first data set is a set of measurement of times between failures of AC
#' equipment in Boeing 720 aircraft. Possible data generating mechanisms that
#' are proposed are Exponential and Gamma.
#'
#' This data can be analyzed via bootstrap methods: parameteric and
#' non-parametric

data("aircondit")
glimpse(aircondit)

mu <- mean(aircondit$hours)
variance <- var(aircondit$hours)

#' We can start by assuming that these data have been drawn at random from a
#' exponential distribution of the above calculated mean (mu). The parametric
#' assumption is the distribution i.e., exponential.
#'
#' One can generate random samples from the assumed exponential distribution,
#' like so (best guess for mean of the distribution comes from the sample
#' itself):

b_expo <- rexp(length(aircondit), rate = 1/mu)

#' Similar random draws can be made from an assumed gamma distribution. If $X$
#' is a distributed gamma with shape $k$ and scale $\theta$, then $E[X] = k
#' \theta$ and $Var(X) = k \theta^2$. In this case, the best information
#' estimate of the mean and the variance of the distribution comes from the
#' sample itself like before. So, $\theta = Var(X)/E[X] = variance/mu$

theta <- variance/mu
k <- mu/theta

b_gamma <- rgamma(length(aircondit), shape = k, scale = theta)

#' In the above analysis, the fitted model in the first case is the exponential
#' distribution and the second case is the gamma distribution.
#'
#' Now, we move on to generating a large number of bootstrap samples and looking
#' at the the estimates from these bootstrap samples.
#'
#' First, we assume that the fitted model is exponential and see how the Bias,
#' i.e., E[T|F] - t and the Variance, i.e., Var(T) change as the number of
#' bootstrap replicates increase. We wish to see the bias move to 0 as the
#' number of replications increase.

boot_df <- function(data, mu, n_reps) {
  bm <- matrix(nrow = length(data$hours), ncol = n_reps)
  
  for (rep in 1:n_reps) {
    bm[, rep] <- rexp(nrow(bm), rate = 1/mu)
  }
  
  colnames(bm) <- paste0("b", 1:n_reps)
  bdf <- as.data.frame(bm)
  
  return(bdf)
}

bias_boot <- function(boot_df, mu) {
  means <- boot_df %>% 
            summarize_all(mean) %>% 
            as_vector()
  
  return(mean(means - mu))
}

variance_boot <- function(boot_df, mu) {
  means <- boot_df %>% 
            summarize_all(mean) %>% 
            as_vector()
  
  return(var(means))
}

#' We repeat the entire bootstrapping exercise 4 times as shown in the book.

estimate_bias <- function(data, mu) {
  boot_10 <- boot_df(data, mu, 10)
  boot_50 <- boot_df(data, mu, 50)
  boot_100 <- boot_df(data, mu, 100)
  boot_500 <- boot_df(data, mu, 500)
  
  return(c(bias_boot(boot_10, mu), 
           bias_boot(boot_50, mu),
           bias_boot(boot_100, mu),
           bias_boot(boot_500, mu)))
}

estimate_variance <- function(data, mu) {
  boot_10 <- boot_df(data, mu, 10)
  boot_50 <- boot_df(data, mu, 50)
  boot_100 <- boot_df(data, mu, 100)
  boot_500 <- boot_df(data, mu, 500)
  
  return(c(variance_boot(boot_10, mu), 
           variance_boot(boot_50, mu),
           variance_boot(boot_100, mu),
           variance_boot(boot_500, mu)))
}

bias_boot_expo <- data.frame(R = c(10, 50 , 100, 500),
                             rep1 = estimate_bias(aircondit, mu), 
                             rep2 = estimate_bias(aircondit, mu),
                             rep3 = estimate_bias(aircondit, mu),
                             rep4 = estimate_bias(aircondit, mu))

variance_boot_expo <- data.frame(R = c(10, 50 , 100, 500),
                                 rep1 = estimate_variance(aircondit, mu), 
                                 rep2 = estimate_variance(aircondit, mu),
                                 rep3 = estimate_variance(aircondit, mu),
                                 rep4 = estimate_variance(aircondit, mu))

#' Figure 2.1
dev.new()
ggplot(bias_boot_expo %>% gather(replicate, bias, -R)) +
  geom_point(aes(x = R, y = bias, group = replicate, color = replicate)) +
  geom_line(aes(x = R, y = bias, group = replicate, color = replicate)) +
  geom_hline(aes(yintercept = 0)) +
  labs(x  = "Number of bootstrap replicates",
       y = "Bias",
       color = "Replicate",
       title = "Variation in bias with number of bootstrap replicates")

ggplot(variance_boot_expo %>% gather(replicate, bias, -R)) +
  geom_point(aes(x = R, y = bias, group = replicate, color = replicate)) +
  geom_line(aes(x = R, y = bias, group = replicate, color = replicate)) +
  labs(x  = "Number of bootstrap replicates",
       y = "Variance",
       color = "Replicate",
       title = "Variation in variance with number of bootstrap replicates") 

#' In the above plot, as expected the variance converges to the variance of an
#' exponential distribution and the bias converges to 0.

