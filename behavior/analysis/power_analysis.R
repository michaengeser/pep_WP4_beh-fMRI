# This script performs a power analysis to estimate the required sample size
# to find a partial correlation of a required power. First the required number 
# of pairs for an inter-subject RSA approach will be calculated and than the minimum
# number of subject to reach this number of pairs

# Preparatiom
install.packages("pwr")
library(pwr)

# Values 
partial_r <- 0.164  # partial correlation coefficient
power <- 0.99      # desired power level
alpha <- 0.05     # significance level

# Compute the required pairs of subjects 
result <- pwr.r.test(r = partial_r, power = power, sig.level = alpha)

# Get sample size (n here is number of possible subject pairs)

# Coefficients for the quadratic equation n^2 - n - 2k = 0
k <- result$n
c <- -2 * k

# Calculate the discriminant
discriminant <- (-4) * c

# Calculate the positive root using the quadratic formula
sample_size <- (1 + sqrt(discriminant)) / 2
ceiling(sample_size)
  
  
