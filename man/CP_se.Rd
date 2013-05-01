\name{CP.se}
\alias{jCP}
\alias{CP.se}

\title{
Compute a conditional probability of observing a set of counts as extreme as the new observations of a subject given the previous observations from the same subject based on the negative binomial mixed effect independent model.
}


\description{
%
Given the parameter estimates of \code{c( alpha,theta,beta0,beta1,...)} of the negative binomial mixed effect independent model,
these functions compute the following conditional probability:

\emph{Pr(q(Y_new) >= q(y_new) | Y_pre = y_pre)}

where \emph{y_new} and \emph{y_pre} are vectors of previous and new observations from a same subject and \code{q()} is a function which provides a scalar summary of the new observations.
%

These functions are subroutines of \code{\link{index.batch}}.
%
When input the maximum likelihood estimates of the paramters, \code{CP.se} returns the estimate of the conditional probability of single subject and the asymptotic standard error of the logit of the estimate of the conditional probability based on the independent model.
The computation for of the probability is done by its subroutine \code{jCP}.
}


\usage{
CP.se(tpar, Y1, Y2, sn1, sn2, XM = NULL, dist = "G", V, pty = "sum")
jCP(tpar, Y1, Y2, sn1, sn2, XM = NULL, dist = "G", LG = FALSE, oth = NULL, type = "sum") 
}



\arguments{
  \item{tpar}{

A vector of length 3 + # covariates, containing the estimates of the model in the order of \code{log(a),log(theta),beta0,beta1...}
If random effects are assumed to be from the log-normal distribution, then \code{theta} is a variance estimate of the random effect.
If random effects are assumed to be from the gamma distribution, then \code{theta} is a \code{shape} parameter estimate.
If the semi-parametric approach is taken, then \code{theta} is a place holder and can be any number.


}
  \item{Y1}{
A scalar containing the sum of the previously observed response counts of a subject.
}
  \item{Y2}{
A scalar containing the summary statistics of the newly observed response counts of a subject \emph{ q(y_new)}
}
  \item{sn1}{
The number of previous observations.
}
  \item{sn2}{
The number of new observations.
}
  \item{XM}{
A \emph{n_i} by # covariates matrix containing the covariate values of subject \emph{i}, where whre \emph{n_i} is the total number of previous and new observations.
If there is no covariate, i.e., the model only has an intercept term, then set \code{XM=NULL}.
}
  \item{dist}{
If \code{dist=="G"}, then the conditional probability is computed by assuming the random effect is from the gamma distribution.
If \code{dist=="N"}, then the conditional probability is computed by assuming the random effect is from the log-normal distribution.

Note that the semiparametric model is NOT accepted in \code{CP.se}. For \code{jCP}, if \code{dist="NoN"}, then the conditional probability is computed by assuming the random effect is from a distribution represented by the argument \code{oth}..
}
  \item{V}{
The variance covariance matrix of the parameter estimates.
}
  \item{pty}{
The summary statistics for the new scans. \emph{q()}. The current options are \code{"sum"} and \code{"max"}.
}

\item{oth}{
The argument only for \code{jCP}.

If \code{dist= "G"} or \code{"N"}, \code{oth=NULL}.
If \code{dist="NoN"}, \code{othr} must be the frequency table of the random effects.
i.e.,\code{dist=obj$gtb} where \code{obj} is the output of \code{\link{mle.a3.fun}}.
}

  \item{LG}{
If \code{LG=TRUE} then the logit of the conditional probability is returned.
}

\item{type}{
Same as \code{pty}.
}
}


\references{
 Zhao, Y., Li, D.K.B., Petkau, J.A., Riddehough, A. & Traboulsee, A.  Detection of unusually large increases in 
MRI lesion counts in multiple sclerosis patients.
}
\author{
 Zhao, Y. and Kondo, Y.
}

\seealso{

The functions to fit the relevant models:

\code{\link{mle.fun}},
\code{\link{mle.ar1.fun}},
\code{\link{mle.a3.fun}},
\code{\link{mle.ar1.non3}},

The other subroutines of \code{\link{index.batch}} to compute the conditional probability index:
%\code{\link{jCP.ar1}},
%\code{\link{CP1.ar1}},
\code{\link{MCCP.ar1}},
\code{\link{CP.ar1.se}},
\code{\link{CP.se}},
\code{\link{jCP}},

The functions to generate simulated datasets:
\code{\link{rNBME.R}}.
}

\examples{

\dontrun{


## tpar contains: log(a),log(theta),beta0
tpar <-  c(0.5, -0.5, 1.3)
## A scalar containing the sum of the response counts in pre scans
Y1 <- 0
## A scalar containing the summary statistics of the response counts in new scans q(y_new)
Y2 <- 1
## The number of scans in the pre scans.
sn1 <- 2
## The number of scans in the new scans.
sn2 <- 3
## the covariate matrix
XM <- NULL
dist <- "G"
## the variance covariance matrix:
V <- matrix(
  c(0.0490673003, -0.0004481864, 0.013279476,
   -0.0004481864,  0.0245814022, 0.001231522,
    0.0132794760,  0.0012315221, 0.023888065),nrow=3)

## the estimate of the conditional probability based on the sum summary statistics and its SE
CP.se(tpar = tpar, Y1=Y1 ,Y2= Y2, sn1 = sn1, sn2 = sn2, XM = XM, dist = dist, V = V, pty = "sum")

## the estimate of the conditional probability based on the max summary statistics and its SE
CP.se(tpar = tpar, Y1=Y1 ,Y2= Y2, sn1 = sn1, sn2 = sn2, XM = XM, dist = dist, V = V, pty = "max")

## jCP calls for CP.se to compute the estimate of the conditional probability
jCP(tpar = tpar, Y1 = Y1, Y2 = Y2, sn1 = sn1, sn2 = sn2, XM = XM, dist = dist, LG = FALSE, oth = NULL, type = "sum") 

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
