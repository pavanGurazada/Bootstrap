#' ---
#' title: "Bootstrapping estimates from a power law fit"
#' author: Pavan Gurazada
#' date: "Feb 2018"
#' output: github_document
#' ---
#' 
#' This script comes to existence because the current bootstrap function in the
#' poweRlaw package is too slow for massive data sets. I found that the power 
#' law fit function from the igraph package is much faster and can be incorporated 
#' into the heart of a fatser bootstrap

library(igraph)
library(boot)

set.seed(20130810)

#' * The bootstrap function *
#' 
#' This function follows the prescription input for the boot package. The data
#' and index are passed as arguments, then the power law fit is executed for the
#' resampled data. The power law exponent and the cut-off are returned. 
#'  

plFit <- function(data, i) {
  resampledData <- data[i]
  plfit <- igraph::fit_power_law(resampledData)
  return (c(plfit$alpha, plfit$xmin))
}

sampleData <- degree(sample_pa(n = 1e5, m = 3, directed = FALSE))

samplePLFit <- fit_power_law(sampleData)
cat("The power law coefficient of the fit is : ", samplePLFit$alpha)

nreps <- 200 # Number of bootstrap samples

#' using the bootstrap function from the boot package offers the advantage of 
#' parallel processing to make the process much faster
#' 

samplePLBoot <- boot(data = sampleData,
                     statistic = plFit, 
                     R = nreps, 
                     parallel = "multicore",
                     ncpus = 3)
