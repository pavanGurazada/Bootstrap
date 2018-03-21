#' ---
#' title: "Exploring the chisquare distribution"
#' author: Pavan Gurazada
#' output: github_document
#' ---
#' last update: Wed Mar 21 10:45:57 2018
#' 

library(tidyverse)
library(recipes)
library(ggthemes)

#' A chisq distribution with `k` degrees of freedom is the sum of `k` standard
#' normal variables. This provides us a way to quickly generate samples from
#' this distribution.
#' 
#' As an example lets generate data for a chi-square distribution with degree
#' of freedom 1, 2, 3, 4

n_samples <- 1e4

x <- rnorm(n_samples, mean = 0, sd = 1)
y <- rnorm(n_samples, mean = 0, sd = 1)
z <- rnorm(n_samples, mean = 0, sd = 1)
t <- rnorm(n_samples, mean = 0, sd = 1)

chi1df <- x^2
chi2df <- x^2 + y^2
chi3df <- x^2 + y^2 + z^2
chi4df <- x^2 + y^2 + z^2 + t^2

data_chisq <- data.frame(chi1df, chi2df, chi3df, chi4df)

dev.new()
ggplot(data_chisq %>% gather(Variable, Value)) +
  geom_histogram(aes(x = Value), fill = "black", color = "white") +
  facet_wrap(~Variable, ncol = 2) +
  labs(x = "Value",
       y = "Count",
       title = "Distribution of Chi-squared random variables",
       subtitle = "(generated from Normal random variables)") +
  theme_few()

#' There is a close relationship between the Poisson distribution and the normal
#' distribution as well. Lets illustrate with more fake data

pois_means <- seq(10, 50, length = 6) # mean and variance of a poisson distribution are the same
pois_data <- data.frame(n_samples)

for (lambda in pois_means) {
  x <- rpois(n_samples, lambda)
  pois_data <- cbind(pois_data, x)
}

pois_data <- pois_data[, -1]
colnames(pois_data) <- c("pois10", "pois18", "pois26", "pois34", "pois42", "pois50")

ggplot(pois_data %>% gather(Variable, Value)) +
  geom_histogram(aes(x = Value), fill = 'black', color = 'white') +
  facet_wrap(~Variable, ncol = 2) +
  labs(x = "Value",
       y = "Count",
       title = "Distribution of Poisson random variables") +
  theme_few()

#' The above plot looks normal-ish, but let us standardize each of them so that 
#' we can compare with the standard normal distribution

ggplot(scale(pois_data) %>% as.data.frame() %>% gather(Variable, Value)) +
  geom_histogram(aes(x = Value), fill = 'black', color = 'white') +
  facet_wrap(~Variable, ncol = 2) +
  labs(x = "Value",
       y = "Count",
       title = "Distribution of standardized Poisson random variables") +
  theme_few()

#' Given this correlation between Poisson and Normal distribution, the summation 
#' of squares of Poisson random variables should also look like chisquare
#' Once standardized there is no real difference between each of the Poisson's 
#' since they share the same mean and variance
#' 
pois_std <- as.data.frame(scale(pois_data))

pois_to_chisq <- pois_std %>% mutate(chisq1 = pois10^2,
                                     chisq2 = pois10^2 + pois18 ^2,
                                     chisq3 = pois10^2 + pois18^2 + pois26^2,
                                     chisq4 = pois10^2 + pois18^2 + pois26^2 + pois34^2) %>% 
                              select(chisq1, chisq2, chisq3, chisq4)

ggplot(pois_to_chisq %>% gather(Variable, Value)) +
  geom_histogram(aes(x = Value), fill = 'black', color = 'white') +
  facet_wrap(~Variable, ncol = 2) +
  labs(x = "Value",
       y = "Count",
       title = "Distribution of Chi-squared random variables",
       subtitle = "(generated from Poisson random variables)") +
  theme_few()

#' The plots from the data generated fromm a Normal back-end and a Poisson
#' back end are remarkably similar!

#' Now consider the example mentioned [here](https://www.r-bloggers.com/exploring-the-underlying-theory-of-the-chi-square-test-through-simulation-part-1/)
#' 
#' There are three buildings with the number of bottles deposited for recycling
#' daily being 20, 40 and 80. So, on an average the total number of bottles to
#' be collected for recycling from these buildings is 140. Each building is
#' Poisson distributed with means mentioned above and the total is also Poisson
#' with mean 140. We can generate 10,000 days of data by drawing from these 
#' individual predictions

bottles_recycling <- data.frame(bin_1 = rpois(n_samples, 20),
                                bin_2 = rpois(n_samples, 40),
                                bin_3 = rpois(n_samples, 80))

bottles_recycling <- bottles_recycling %>% mutate(bins_total = bin_1 + bin_2 + bin_3) 

ggplot(bottles_recycling %>% gather(Variable, Value)) +
  geom_point(aes(x = Variable, y = Value), position = "jitter") +
  labs(x = "Building",
       y = "Bottles for recycling",
       title = "Bottles collected from each building and in total for 10,000 days") +
  theme_few()

#' Now we wish to subset the data for a total of bottles collected to be between
#' 138 and 142

bottles_collected_subset <- bottles_recycling %>% filter(bins_total >= 138 & bins_total <= 142)
summarize_all(bottles_collected_subset, mean)

#' As we can see, this does not change the means of the individual bins or the
#' total However, the variance is considerably reduced as can be seen from the
#' following plot

ggplot(bottles_recycling %>% gather(Variable, Value)) +
  geom_point(aes(x = Variable, y = Value), position = "jitter", color = "gray") +
  geom_point(data = bottles_collected_subset %>% gather(Variable, Value),
             aes(x = Variable, y = Value), position = "jitter", color = "red") +
  labs(x = "Building",
       y = "Bottles for recycling",
       title = "Bottles collected from each building and in total for 10,000 days",
       subtitle = "(red points show the subset of the data where total bottles collected are between 138 and 142)") +
  theme_few()

#' Thinking along the lines of what is hinted in the article, when we select 
#' a sample of the data (for e.g., when we split a data set into training and 
#' testing), the smaller subset of data necessarily has a smaller variance and 
#' using the reduced variance to standardize this data
#' 
#' So, while the overall distribution remains chi-squared, we can see that 
#' constrained samples from a distribution are a trouble
#' 
#' This is the reason why, one uses the mean and the variance of the training
#' set to standardize the test set!