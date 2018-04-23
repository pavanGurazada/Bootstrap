#' ---
#' title: "Chapter 4: Tests"
#' author: Daivson and Hinkley
#' output: github_document
#' ---

#' last update: Mon Apr 23 05:23:39 2018

library(boot)
library(tidyverse)
library(ggthemes)

set.seed(20130810)
theme_set(theme_few() + 
            theme(plot.title = element_text(face="bold")))

#' We collected 50 measurements of fir seedlings in small quadrats. Our null
#' hypothesis is that these are from a Poisson distribution of an unknown mean

data("fir")
glimpse(fir)

#' As with Poisson distribution problems, the usual suspect is over-dispersion,
#' i.e., the variance in the data is more than what we expect from a usual
#' Poisson distribution (whose mean and variance are equal).

mu <- mean(fir$count)
variance <- var(fir$count)

if (variance > mu) print("Overdispersed")

#' It appears that the situation is not that dire. To check if the data is
#' really overdispersed we use the dispersion index as the test statistic
#' 
#' If we conduct the usual bootstrap sampling:

dispersion_index <- function(data_vec, indices) {
  
  y_vec <- data_vec[indices]
  
  return(sum((y_vec - mean(y_vec))^2)/mean(y_vec))
}

b <- boot(data = fir$count, 
          statistic = dispersion_index,
          R = 1000,
          parallel = "snow")

dispersion_boot <- as.data.frame(b$t)
colnames(dispersion_boot) <- c("dispersion_stat")

#' Figure 4.1 (not the one in the book)
dev.new()
ggplot(dispersion_boot) +
  geom_histogram(aes(x = dispersion_stat), fill = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = 55.15), linetype = "dashed", size = 1) +
  labs(x = expression(t ^ "*"),
       y = "Count",
       title = "Bootstrap distribution of the dispersion statistic")

#' We now write a multinomial bootstrap function. In this scenario, the sum of
#' the measurements in the resampled data is the same as the original data. This
#' condition was not true in the case of the previous bootstrap calculation.

table(fir$count) 
bm <- rmultinom(n = 1000, size = sum(fir$count), prob = table(fir$count) %>% as.vector()) # incorrect since the sums will not tally to 107

#' *Permutation tests*
#' 
#' The objective in measuring the handedness of individuals is to see
#' if the two kinds of measures - genetic and integer - are independent. 

data("claridge")
glimpse(claridge)

#' Figure 4.6
ggplot(claridge) +
  geom_point(aes(x = dnan, y = hand)) +
  labs(x = "Genetic measure",
       y = "Integer measure",
       title = "Scatter plot of pairs of handedness measurements",
       subtitle = "(n = 37)")

#' A good meeasure to check independence is the correlation coefficient
sample_corr <- with(claridge, cor(dnan, hand))

#' Under the null hypothesis, the two sets of observations are independent. So
#' this should hold under the permutation test where we generate all possible
#' permutations of the `nan` values (keeping the order statistics of the `dnan` 
#' variable constant). For each such data set, we compute the correlation coefficient
#' and check the distribution of these values.

R <- 1000

u <- sort(claridge$dnan)
corr_boot <- numeric(R) # stores the bootstrap estimates from the permutations

for (r in 1:R) {
  x <- sample(x = claridge$hand, 
              size = length(claridge$hand), 
              replace = FALSE) # generates a permutation
  
  corr_boot[r] <- cor(u, x)
}

corr_boot <- data.frame(correlation = corr_boot)

#' Figure 4.7
ggplot(corr_boot) +
  geom_histogram(aes(x = correlation), fill = "black", alpha = 0.5) +
  geom_vline(aes(xintercept = sample_corr), linetype = "dashed", size = 1) +
  labs(x = "Correlation",
       y = "Count",
       title = "Histogram of correlation from permutations")

#' From Figure 4.7 it is clear that the chance of the two being un-correlated is
#' faint, i.e., we reject the null hypothesis. Given that the null is true the
#' probabilty of observing the sample correlation is very less, leading us to
#' the conclusion that our assumption is incorrect.
#' 
#' Now, we return to the case of the non-parametric bootstrap. For the gravity 
#' data we are investigating the null hypothesis that there is no difference in the 
#' measurements from the different groups of measurements of `g`

data(gravity)
glimpse(gravity)

#' Lets pull out two sample set of observations

gravity %>% 
  filter(series == 7) %>% 
  pull(g) -> sample7

gravity %>% 
  filter(series == 8) %>% 
  pull(g) -> sample8

#' Under the null hypothesis, there is no difference between sample7 and sample8
#' So, we can do a pooled resampling, i.e., sampling with replacement within the
#' two subgroups and return the difference in means as the test statistic. We
#' then check how likely is the difference to be what we observe in the sample
#' under the bootstrap distribution
#'
#' We begin with the data subset, resample a total of 26 rows from the data with
#' replacement and check the distribution of the diffeernce in means

gravity_df <- filter(gravity, series == 7 | series == 8)

out_vec <- numeric(1000) # To hold the bootstrapped mean differences 

for (r in 1:1000) {
  gravity_df %>% 
    slice(sample(x = 1:nrow(gravity_df), size = nrow(gravity_df), replace = TRUE)) %>% 
    group_by(series) %>% 
    summarize_at(vars(g), function(col) mean(col)) %>% 
    summarize_at(vars(g), function(col) col[2] - col[1]) %>% 
    as.numeric() -> 
    out_vec[r]
}

qplot(out_vec, geom = "histogram") +
  geom_vline(aes(xintercept = mean(sample8) - mean(sample7)), 
             linetype = "dashed", 
             size = 1) +
  labs(x = "Difference in means",
       y = "Count",
       title = "Distribution of difference in means from the bootstrap sample")

meandiff_boot <- function(s1, s2, R = 1000) {
  
  out_vec <- numeric(R) # allocate a vector to put in the bootstrapped differences
  
  for (r in 1:R) {
    resamp1 <- sample(x = s1, size = length(s1), replace = TRUE)
    resamp2 <- sample(x = s2, size = length(s2), replace = TRUE)
    
    out_vec[r] <- mean(resamp2) - mean(resamp1) 
  }
  
  return(out_vec)
}

b <- meandiff_boot(sample7, sample8)

#' Figure 4.9
dev.new()
qplot(b, geom = "histogram") +
  geom_vline(aes(xintercept = mean(sample8) - mean(sample7)), 
             linetype = "dashed", 
             size = 1) +
  labs(x = "Difference in means",
       y = "Count",
       title = "Distribution of difference in means from the bootstrap sample")

#' Now let us look at the null hypothesis that there is no difference between
#' the means from the eight series of measurements. We summarize the data in the
#' table below.

gravity %>% 
  group_by(series) %>% 
  summarize(n = n(),
            mean_y = mean(g),
            var_y = var(g),
            w_y = n/var_y) ->
  gravity_stats 

gravity_stats %>% 
  summarize(mu0 = sum(w_y * mean_y)/sum(w_y)) %>% 
  as.numeric() ->
  mu0

#' Table 4.3
gravity_stats %>% 
  mutate(sigma_y2 = (n-1)/n * var_y + (mean_y - mu0)^2,
         resid_denom = sqrt(sigma_y2 - 1/sum(w_y))) %>% 
  select(series, n, mean_y, var_y, sigma_y2, w_y, resid_denom) ->
  gravity_stats

#' We bootstrap from the model: y_ij = mu0 + sigma_y * e_ij
#' 
#' Here we sample with replacement from the residuals

bm <- data.frame(series = c(rep(1, gravity_stats$n[1]), rep(2, gravity_stats$n[2]),
                            rep(3, gravity_stats$n[3]), rep(4, gravity_stats$n[4]),
                            rep(5, gravity_stats$n[5]), rep(6, gravity_stats$n[6]),
                            rep(7, gravity_stats$n[7]), rep(8, gravity_stats$n[8])),
                 y_ij = gravity$g,
                 mu0 = rep(mu0, sum(gravity_stats$n)),
                 sigmai0 = c(rep(sqrt(gravity_stats$sigma_y2[1]), gravity_stats$n[1]), 
                             rep(sqrt(gravity_stats$sigma_y2[2]), gravity_stats$n[2]),
                             rep(sqrt(gravity_stats$sigma_y2[3]), gravity_stats$n[3]),
                             rep(sqrt(gravity_stats$sigma_y2[4]), gravity_stats$n[4]),
                             rep(sqrt(gravity_stats$sigma_y2[5]), gravity_stats$n[5]),
                             rep(sqrt(gravity_stats$sigma_y2[6]), gravity_stats$n[6]),
                             rep(sqrt(gravity_stats$sigma_y2[7]), gravity_stats$n[7]),
                             rep(sqrt(gravity_stats$sigma_y2[8]), gravity_stats$n[8])),
                 resid_denom = c(rep(gravity_stats$resid_denom[1], gravity_stats$n[1]), 
                                 rep(gravity_stats$resid_denom[2], gravity_stats$n[2]),
                                 rep(gravity_stats$resid_denom[3], gravity_stats$n[3]),
                                 rep(gravity_stats$resid_denom[4], gravity_stats$n[4]),
                                 rep(gravity_stats$resid_denom[5], gravity_stats$n[5]),
                                 rep(gravity_stats$resid_denom[6], gravity_stats$n[6]),
                                 rep(gravity_stats$resid_denom[7], gravity_stats$n[7]),
                                 rep(gravity_stats$resid_denom[8], gravity_stats$n[8])))

bm <- mutate(bm, e_ij = (y_ij - mu0)/resid_denom)

boot_samples <- matrix(nrow = nrow(bm), ncol = 1000 + 3)
boot_samples[, 1] <- bm$series
boot_samples[, 2] <- bm$mu0
boot_samples[, 3] <- bm$sigmai0

for (r in 1:1000) {
  boot_samples[, r + 3] <- 
}