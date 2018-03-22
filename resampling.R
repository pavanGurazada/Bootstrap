#' ---
#' title: "`rsample` package to conduct resampling of data"
#' author: Pavan Gurazada
#' output: github_document
#' ---
#' last update: Thu Mar 22 14:34:48 2018

library(rsample)
library(tidyverse)
library(mlbench)
library(pryr)

set.seed(20130810)

#' The premise of this package is rather interesting and worth exploring
#' When we create bootstrap samples of a data set the bootstraps should not hold
#' a multiple of the memory footprint of the original data se.

data("LetterRecognition")
object_size(LetterRecognition)

boots <- bootstraps(LetterRecognition, times = 1000)

object_size(boots) != 1000 * object_size(LetterRecognition)

class(boots)

#' The main class in this package is the `rset` object which is a collection of 
#' resamples

bt_resamples <- bootstraps(mtcars, times = 3)
glimpse(bt_resamples)

first_resample <- bt_resamples$splits[[1]]
head(as.data.frame(first_resample))
head(as.data.frame(first_resample, data = "assessment"))

data("attrition")
names(attrition)

table(attrition$Attrition)

# Is there a difference in income between the genders?

dev.new()
ggplot(attrition, aes(x = Gender, y = MonthlyIncome)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()

median_diff <- function(splits) {
  x <- analysis(splits)
  median(x$MonthlyIncome[x$Gender == "Female"]) - 
    median(x$MonthlyIncome[x$Gender == "Male"])
}

set.seed(20130810)
bt_resamples <- bootstraps(attrition, times = 2000)
bt_resamples$wage_diff <- map_dbl(bt_resamples$splits, median_diff)

ggplot(bt_resamples, aes(x = wage_diff)) +
  geom_line(stat = "density", adjust = 1.25) +
  xlab("Difference in Median Monthly Income (Female - Male)") +
  theme_bw()
