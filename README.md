# MultSurvTests
This R package contains multivariate two-sample survival permutation tests, based on the logrank and Gehan statistics. The implemented tests are described in [Persson et al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.1002/pst.1938).

Installation:

```
library(devtools)
install_github("lukketotte/MultSurvTests")
```

Example usage, comparing the bivariate survival times of the two treatment groups in the `diabetes` data (included in the package):

```
library(MultSurvTests)
# Diabetes data:
?diabetes

# Survival times for the two groups:
x <- as.matrix(subset(diabetes, LASER==1)[c(6,8)])
y <- as.matrix(subset(diabetes, LASER==2)[c(6,8)])

# Censoring status for the two groups:
delta.x <- as.matrix(subset(diabetes, LASER==1)[c(7,9)])
delta.y <- as.matrix(subset(diabetes, LASER==2)[c(7,9)])

# Create the input for the test:
z <- rbind(x, y)
delta.z <- rbind(delta.x, delta.y)

# Run the tests with 99 permutations:
perm_gehan(B = 99, z, delta.z, n1 = nrow(x))
perm_mvlogrank(B = 99, z, delta.z, n1 = nrow(x))
```
