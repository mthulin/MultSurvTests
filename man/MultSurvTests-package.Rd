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
This needs to be fixed
}

\usage{
gehan(x, y, delta.x, delta.y, n1, n2, p)
mvlogrank(x, y, delta.x, delta.y, n1, n2, p)
perm_gehan(B, z, delta.z, n1, n2, g.test)
perm_mvlogrank(B, z, delta.z, n1, n2, mv.test)
}

\arguments{
  \item{x}{a numeric matrix}
  \item{y}{a numeric matrix}
  \item{delta.x}{a numeric matrix}
  \item{delta.y}{a numeric matrix}
  \item{n1}{an integer, nrow(x)}
  \item{n2}{an integer, nrow(y)}
  \item{p}{an integer}
  \item{B}{an integer, number of permutations}
  \item{z}{a numeric matrix (x and y combined)}
  \item{z.delta}{a numeric matrix (delta.x and delta.y             combined)}
  \item{g.test}{a numeric, result from gehan}
  \item{mv.test}{a numeric, result from mvlogrank}
}

\value{
\code{gehan()} Returns a 1x1 matrix containing a numeric value

\code{mvlogrank()} Returns a 1x1 matrix containing a numeric value

\code{perm_gehan()} Returns a numeric as a result from the permutation

\code{perm_mvlogrank()} Returns a numeric as a result from the permutation
}

\author{
\packageAuthor{MultSurvTests}

Maintainer: \packageMaintainer{MultSurvTests}
}
\references{
L.J. Wei & J.M. Lachin, Two-Sample Asymptotically Distribution-Free Tests for Incomplete Multivariate Observations, 1984, Journal of the American Statistical Association, Vol. 79, No.387
}

\keyword{ package }
\seealso{
 \code{Perhaps link to github repository}
}
\examples{
data <- as.matrix(wltestdata)
x <- data[1:23, c(2, 4, 6, 8)]
y <- data[24:47, c(2, 4, 6, 8)]
delta.x <- data[1:23,c(3, 5, 7, 9)]
delta.y <- data[24:47, c(3, 5, 7,9)]
z <- rbind(x, y)
delta.z <- rbind(delta.x, delta.y)

n1<-dim(x)[1]
n2<-dim(y)[1]
p<-dim(x)[2]

g <- gehan(x, y, delta.x, delta.y, n1, n2, p)
mv <- mvlogrank(x, y, delta.x, delta.y, n1, n2, p)

perm_gehan(B = 100, z, delta.z, n1, n2, p, g)
perm_mvlogrank(B = 100, z, delta.z, n1, n2, p, g)

}
