#' ---
#' title: "Chapter 2: The Basic Bootstraps"
#' author: Davison and Hinkley
#' output: github_document
#' ---

#' last update: Wed Apr 18 10:01:45 2018

library(boot)
library(tidyverse)
library(ggthemes)

library(parallel)
library(doParallel)
registerDoParallel(cores=detectCores(all.tests=TRUE))

set.seed(20130810)

theme_set(theme_few())

#' Here, we are concerned about bootstrap estimates from a single, homonegenous
#' set of data
#'
#' The first data set is a set of measurement of times between failures of AC
#' equipment in Boeing 720 aircraft. Possible data generating mechanisms that
#' are proposed are Exponential and Gamma.

data("aircondit")
glimpse(aircondit)

mu <- mean(aircondit$hours)
variance <- var(aircondit$hours)

#' This data can be analyzed via bootstrap methods: parameteric and
#' non-parametric.
#' 
#' *Parametric bootstrap analysis*
#' 
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
#'
#' A sense of the variation in the bootstrap process can be obtained from
#' histograms of the statistic computed on the bootstrap samples

mean_boot_99 <- boot_df(aircondit, mu, 99) %>% 
                  summarize_all(mean) %>% 
                  as_vector()

mean_boot_999 <- boot_df(aircondit, mu, 999) %>% 
                  summarize_all(mean) %>% 
                  as_vector()

#' Figure 2.3
dev.new()
ggplot(data.frame(t = mean_boot_99)) +
  geom_histogram(aes(x = t), fill = "black", color = "white") +
  labs(x = "Mean of bootstrap sample",
       y = "Count",
       title = "Distribution of bootstrap sample means (R = 99)")

ggplot(data.frame(t = mean_boot_999)) +
  geom_histogram(aes(x = t), fill = "black", color = "white") +
  labs(x = "Mean of bootstrap sample",
       y = "Count",
       title = "Distribution of bootstrap sample means (R = 999)")

#' *Non-parametric boootstrap analysis*
#' 
#' Non-parametric bootstrap works well when we do not have a theoretical reason
#' to expect a well-known distribution as a representative of the data
#' generating process. This often happens when we are interested in an outcome
#' that is a composition of other columns in the data

data("city")
glimpse(city)

data("bigcity")
glimpse(bigcity)

#' This data set containts the populations in 1920 and 1930 for a selection of
#' cities, the smaller data set is the one that is discussed in chapter 2. We
#' are interested in the ratio of the mean populations in the two years, i.e.,
#' (E[X]/E[U]).
#'
#' There is no parametric model, we could use to model this data, so we use the
#' ratio of the sample averages and build an empirical distribution using
#' bootstrap methods.
#'
#' The cool thing is that even though we do not have a 'name' for the
#' distribution followed by the ratio, we can use resampling to generate
#' approximate values.
#'
#' We begin as usual with a boot function that takes the data and resamples it.
#' Here we have a pair of values to be sampled so indices are used (it is most
#' convenient to use the `boot` function!)
#' 

mean_ratio <- function(data, indices) {
  return(mean(data[indices, ]$x)/mean(data[indices, ]$u))
}

b <- boot(city, mean_ratio, R = 999, parallel = "snow")

boot_data <- as.data.frame(b$t)
colnames(boot_data) <- c("mean_ratio")

#' Figure 2.5
dev.new()
ggplot(boot_data) +
  geom_histogram(aes(x = mean_ratio), fill = "black", color = "white") +
  labs(x = expression(bar(X)/bar(U)),
       y = "Count",
       title = "Bootstrap distribution of the mean ratio")

#' The skew is rather apparent from Figure 2.5
#'
#' This brings us naturally to a comparison between parametric and
#' non-parametric methods. We can make a scatter plot of the mean and standard
#' deviation of the bootstrapped samples to see the effect of using both these
#' methods

means <- function(data, indices) {
  return(mean(data[indices, "hours"]))
}

sds <- function(data, indices) {
  return(sd(data[indices, "hours"]))
}

#' non-parametric bootstrap, R = 99
b_means <- boot(aircondit, statistic = means, R = 99, parallel = "snow")
b_sds <- boot(aircondit, statistic = sds, R = 99, parallel = "snow")

boot_mean_nonparam <- as.data.frame(b_means$t)
colnames(boot_mean_nonparam) <- c("mean_time_nonparam")

boot_sd_nonparam <- as.data.frame(b_sds$t)
colnames(boot_sd_nonparam) <- c("sd_time_nonparam")

boot_nonparam <- bind_cols(boot_mean_nonparam, boot_sd_nonparam)

#' parametric bootstrap, R = 99
boot_mean_param <- as.data.frame(mean_boot_99)
colnames(boot_mean_param) <- c("mean_time_param")

boot_sd_param <- boot_df(aircondit, mu, 99) %>% 
                  summarize_all(sd) %>% 
                  as_vector() %>% 
                  as.data.frame() 

colnames(boot_sd_param) <- c("sd_time_param")

boot_param <- bind_cols(boot_mean_param, boot_sd_param)

#' Figure 2.9

dev.new()
ggplot(boot_param) +
  geom_point(aes(x = mean_time_param, y = sd_time_param)) +
  labs(x = "Bootstrap Average",
       y = "Bootstrap SD",
       title = "Variation across the bootstrap samples (parametric, R = 99)")

ggplot(boot_nonparam) +
  geom_point(aes(x = mean_time_nonparam, y = sd_time_nonparam)) +
  labs(x = "Bootstrap Average",
       y = "Bootstrap SD",
       title = "Variation across the bootstrap samples (non-parametric, R = 99)")


#' non-parametric bootstrap, R = 999
b_means <- boot(aircondit, statistic = means, R = 999, parallel = "snow")
b_sds <- boot(aircondit, statistic = sds, R = 999, parallel = "snow")

boot_mean_nonparam <- as.data.frame(b_means$t)
colnames(boot_mean_nonparam) <- c("mean_time_nonparam")

boot_sd_nonparam <- as.data.frame(b_sds$t)
colnames(boot_sd_nonparam) <- c("sd_time_nonparam")

boot_nonparam <- bind_cols(boot_mean_nonparam, boot_sd_nonparam)

#' parametric bootstrap, R = 999
boot_mean_param <- as.data.frame(mean_boot_999)
colnames(boot_mean_param) <- c("mean_time_param")

boot_sd_param <- boot_df(aircondit, mu, 999) %>% 
                  summarize_all(sd) %>% 
                  as_vector() %>% 
                  as.data.frame() 

colnames(boot_sd_param) <- c("sd_time_param")

boot_param <- bind_cols(boot_mean_param, boot_sd_param)

#' Figure 2.9

dev.new()
ggplot(boot_param) +
  geom_point(aes(x = mean_time_param, y = sd_time_param)) +
  labs(x = "Bootstrap Average",
       y = "Bootstrap SD",
       title = "Variation across the bootstrap samples (parametric, R = 999)")

ggplot(boot_nonparam) +
  geom_point(aes(x = mean_time_nonparam, y = sd_time_nonparam)) +
  labs(x = "Bootstrap Average",
       y = "Bootstrap SD",
       title = "Variation across the bootstrap samples (non-parametric, R = 999)")


