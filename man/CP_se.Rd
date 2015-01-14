\name{CP.se}
\alias{jCP}
\alias{CP.se}

\title{
Compute a conditional probability of observing a set of counts as extreme as the new observations of a subject given the previous observations from the same subject based on the negative binomial mixed effect independent model.
}


\description{
Given the parameter estimates of 
\eqn{\alpha,\theta, \beta_0, \beta_1,\cdots}{\alpha,\theta, \beta[0], \beta[1],...} 
of the negative binomial mixed effect AR(1) model,
these functions compute the following conditional probability: 


\eqn{ 
Pr(q(\boldsymbol{Y}_{i,new}) \ge q(\boldsymbol{y}_{i,new})| \boldsymbol{Y}_{i,pre}=\boldsymbol{y}_{i,pre})
}{
Pr(q(Y[i,new])>=q(y[i,new])|Y[i,pre]=y[i,pre])
}


where \eqn{\boldsymbol{y}_{i,new}}{y[i,new]} and \eqn{y_{i,pre}}{y[i,pre]} are vectors of previous and new observations from subject \eqn{i}{i} and \eqn{q()}{q()} is a function which provides a scalar summary of the new observations. These functions are subroutines of \code{\link{index.batch}}. \code{CP.se} returns the estimate of the conditional probability of single subject and the asymptotic standard error of the logit of the estimate of the conditional probability based on the independent model. The computation for of the probability is done by its subroutine \code{jCP}.
}


\usage{
CP.se(tpar, Y1, Y2, sn1, sn2, XM = NULL,
      RE = "G", V, qfun ="sum",i.tol)
jCP(tpar, Y1, Y2, sn1, sn2, XM = NULL, 
    RE = "G", LG =FALSE, oth = NULL, qfun = "sum",i.tol)
}



\arguments{
  \item{tpar}{
A vector of length 3 + # covariates, containing the estimates of the model in the order of 
\eqn{\log(\alpha),\log(\theta),\beta_0,\beta_1,\cdots}{\log(\alpha),\log(\theta),\beta[0],\beta[1],...}.
Note that \eqn{\alpha}{\alpha} is the dispersion parameter and the \eqn{\theta}{\theta} is a variance estimate of the random effect.
If the semi-parametric approach is taken, then \eqn{\log(\theta)}{\log(\theta)} is a place holder and can be any number.
}
  \item{Y1}{
A scalar containing the sum of the previously observed response counts of a subject. 
}
  \item{Y2}{
A scalar containing the summary statistics of the newly observed response counts of a subject \eqn{ q(y_{i,new})}{ q(y[i,new])}
}
  \item{sn1}{
The number of previous observations.
}
  \item{sn2}{
The number of new observations.
}
  \item{XM}{
A \eqn{n_i}{n[i]} by # covariates matrix containing the covariate values of subject \eqn{i}{i}, where whre \eqn{n_i}{n[i]} is the total number of previous and new observations. If there is no covariate, i.e., the model only has an intercept term, then set \code{XM=NULL}.
}
  \item{RE}{
See \code{\link{lmeNB}}.
Note that the semiparametric model is NOT accepted in \code{CP.se}. For \code{jCP}, if \code{dist="NoN"}, then the conditional probability is computed by assuming the random effect is from a distribution represented by the argument \code{oth}.
}
  \item{V}{
The variance covariance matrix of the parameters \code{tpar}.
}
  \item{qfun}{
See \code{\link{index.batch}}.
}

\item{oth}{
If \code{RE= "G"} or \code{"N"}, then \code{oth} should be {NULL}.
If \code{RE="NoN"}, \code{oth} must be the frequency table of the random effects.
i.e.,\code{RE=obj$gtb} where \code{obj} is the output of \code{\link{fitSemiIND}}.
}

  \item{LG}{
If \code{LG=TRUE}, then the logit of the conditional probability is returned.
}
\item{i.tol}{See \code{\link{lmeNB}}}

}


\references{
Detection of unusual increases in MRI lesion counts in individual multiple sclerosis patients. (2013) Zhao, Y., Li, D.K.B., Petkau, A.J., Riddehough, A., Traboulsee, A., Journal of the American Statistical Association.
}
\author{
 Zhao, Y. and Kondo, Y.
}

\seealso{
The main function to fit the model is:
\code{\link{lmeNB}},

The internal functions of \code{\link{lmeNB}} for fitting relevant models:

\code{\link{fitParaIND}},
\code{\link{fitParaAR1}},
\code{\link{fitSemiIND}},
\code{\link{fitSemiAR1}},

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
RE <- "G"
## the variance covariance matrix:
V <- matrix(
  c(0.0490673003, -0.0004481864, 0.013279476,
   -0.0004481864,  0.0245814022, 0.001231522,
    0.0132794760,  0.0012315221, 0.023888065),nrow=3)

## the estimate of the conditional probability based on the sum summary statistics and its SE
CP.se(tpar = tpar, Y1=Y1 ,Y2= Y2, sn1 = sn1, sn2 = sn2, XM = XM, RE = RE, V = V, qfun = "sum")

## the estimate of the conditional probability based on the max summary statistics and its SE
CP.se(tpar = tpar, Y1=Y1 ,Y2= Y2, sn1 = sn1, sn2 = sn2, XM = XM, RE = RE, V = V, qfun = "max")

## jCP calls for CP.se to compute the estimate of the conditional probability
jCP(tpar = tpar, Y1 = Y1, Y2 = Y2, sn1 = sn1, sn2 = sn2, 
    XM = XM, RE = RE, LG = FALSE, oth = NULL, qfun = "sum")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
