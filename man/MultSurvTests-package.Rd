\name{MultSurvTests-package}
\alias{MultSurvTests-package}
\alias{MultSurvTests}
\docType{package}
\title{
\packageTitle{MultSurvTests}
}
\description{
\packageDescription{MultSurvTests}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MultSurvTests}
\packageIndices{MultSurvTests}
}

\author{
\packageAuthor{MultSurvTests}

Maintainer: \packageMaintainer{MultSurvTests}
}
\references{
L.J. Wei & J.M. Lachin, Two-Sample Asymptotically Distribution-Free Tests for Incomplete Multivariate Observations, 1984, Journal of the American Statistical Association, Vol. 79, No.387

I. Persson, L. Arnroth & M. Thulin, Multivariate two‐sample permutation tests for trials with multiple time‐to‐event outcomes, 2019,  Pharmaceutical statistics, Vol. 18, No.4, pp. 476-485
}

\keyword{ package }
\seealso{
 \code{https://github.com/lukketotte/MultSurvTests}
}
\examples{
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
}
