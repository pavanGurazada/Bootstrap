library(lsr)
library(boot)

set.seed(20130810)

# This is a fist attempt to compute bootstrap estimates of effect sizes in a  
# traditional ANOVA setting. 

# The first function computes the bootstrapped effect sizes and the second prints
# a traditional R-style ANOVA table.Most of the code is for pretty printing
# Will add in more nuances as the code evolves.
# 
# Sample data needs to be added
# 

esBoot <- function(formula, data, i) {
  
  resampledData <- data[i,] # Resample rows
  
  model <- aov(as.formula(formula), data = resampledData)
  
  es <- lsr::etaSquared(model, type = 1)[, "eta.sq"] # Extract effect sizes
  
  return (as.vector(es))
}

bootES <- function(data, formula, conf.int = 0.95, reps = 1000) {
  
  # Fit non-bootstrapped ANOVA to obtain effect sizes and variable names
  aovFit <- aov(as.formula(formula), data = data)
  # Get DV name
  dvName <- colnames(model.frame(aovFit))[1]                    
  # Get IV names
  ivNames <- attr(aovFit$terms , "term.labels")
  
  # Fit ANOVA model with Type I Sum of Squares
  anovaFit <- anova(aovFit)
  
  
  # Call esboot function
  esBootRun <- boot(statistic = esBoot, 
                    formula = formula, 
                    data = data, 
                    R = reps)
  
  esValues <- summary(esBootRun)$original # original effect sizes
  esSE <- summary(esBootRun)$bootSE       # Bootstrap SE for effect sizes
  
  # Get Bootstrapped confidence intervals
  ciLower <- NULL
  ciUpper <- NULL
  
  for(i in 1:length(esValues)){
    bootBCA <- boot.ci(esBootRun, index = i, type = "bca", conf = conf.int)$bca
    ciLower[i] <- bootBCA[4]
    ciUpper[i] <- bootBCA[5]     
  }
  
  # Column names
  cnames <- c("Effect Size", "SE", "LB", "UB")
  
  # Row names
  rnames <- c(ivNames)
  
  # Turn warnings off as the variables are not equal length when creating output.df
  options(warn = -1)
  
  # Create data frame to use as output
  output <- as.data.frame(cbind(esValues, 
                                esSE, 
                                ciLower,  
                                ciUpper))
  
  
  # Remove redundant values for the Residuals row
  #output[nrow(output), 3:10] <- NA
  
  # Set column and row names
  colnames(output) <- cnames
  rownames(output) <- rnames
  
  # Print results
  cat("\nBootstrap Effect size Table", 
      "\n",
      "\n", "Response: ", dvName,
      "\n",
      sep ="")
  
  print(as.matrix(output), justify = "right", na.print = "" , quote = FALSE )
  
  # Turn warnings on
  options(warn = 0)
}

