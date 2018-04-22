#' ---
#' title: "Chapter 3: Further Ideas "
#' author: Davison and Hinkley
#' output: github_document
#' ---
#' last update: Fri Apr 20 06:33:55 2018s
#' 

library(boot)
library(tidyverse)
library(ggthemes)
library(ggfortify)

library(survival)

set.seed(20130810)
theme_set(theme_few() + 
            theme(plot.title = element_text(face="bold")))

data("gravity")
glimpse(gravity)

#' Figure 3.1
dev.new()
ggplot(gravity) +
  geom_boxplot(aes(x = series, y = g), fill = "black", color = "grey")

#' As Figure 3.1 indicates, there seems to be wide variation in the value of g
#' based on the location. This makes it difficullt to estimate the value of `g`
#' based on this data. One good estimator is to use a weighted mean of the
#' values of `g` within each series, with each weight being the reciprocal of
#' the variance. This is a variance stabilizing transformation that is otherwise
#' quite popular.
#' 
#' We use the `boot` library for resampling.

mean_weightd <- function(data, indices) {
  
  wm <- data %>% 
         slice(indices) %>% 
         group_by(series) %>% 
         summarize_at(vars(g), funs(mean, var), na.rm = TRUE) %>% 
         mutate(m = mean/var, wt = 1/var) %>% 
         summarize(wght_mean = sum(m)/sum(wt)) %>% 
         as_vector()
  
  return(wm)
}

b <- boot(data = gravity,
          statistic = mean_weightd,
          R = 1000,
          strata = gravity$series,
          parallel = "snow")

mean_boot <- as.data.frame(b$t)
colnames(mean_boot) <- c("mean_weightd")

ggplot(mean_boot) +
  geom_histogram(aes(x = mean_weightd), 
                 fill = "black", 
                 color = "white", 
                 bins = 20)

#' As can be seen from the plot above, the distribution is really nice now. We
#' use the strata to resample, but for a point estimate we stabilize the
#' variance using a transformation.
#' 
#' *Censored data*

data("aml", package = "boot")
glimpse(aml)

#' Patients are randomly divided into two groups post-development of the
#' condition and the time till symptoms recur are measured.

fit <- survfit(Surv(time, cens) ~ group, data = aml)

#' Figure 3.3
dev.new()
autoplot(fit, facets = TRUE) +
  labs(x = "Time",
       y = "Survival Probability")

summary(fit, time = 20)

#' Our interest is the median remission times of the two groups

st <- function(data) {
  fit <- survfit(Surv(time, cens) ~ group, data)
  out <- NULL
  st <- 1
  
  for (s in 1:length(fit$strata)) {
    inds <- st:(st+fit$strata[s]-1)
    md <- min(fit$time[inds[1-fit$surv[inds] >= 0.5]])
    st <- st + fit$strata[s]
    out <- c(out, md)
  }
  
  return(out)
}

b <- censboot(aml, st, R = 499, strata = aml$group)

