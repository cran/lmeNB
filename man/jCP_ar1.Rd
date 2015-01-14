\name{CP.ar1.se}
\alias{CP.ar1.se}
\alias{jCP.ar1}
\alias{CP1.ar1}
\alias{MCCP.ar1}

\title{
Compute a conditional probability of observing a set of counts as extreme as the new observations of a subjectvisit given the previous observations of the same subject based on the negative binomial mixed-effect AR(1) model.
}

\description{
%
Given the parameter estimates of 
\eqn{\alpha,\theta, \delta, \beta_0, \beta_1,\cdots}{\alpha,\theta, \delta, \beta[0], \beta[1],...} 
of the negative binomial mixed effect AR(1) model,
these functions compute the following conditional probability: 


\eqn{ 
Pr(q(\boldsymbol{Y}_{i,new}) \ge q(\boldsymbol{y}_{i,new})| \boldsymbol{Y}_{i,pre}=\boldsymbol{y}_{i,pre})
}{
Pr(q(Y[i,new])>=q(y[i,new])|Y[i,pre]=y[i,pre])
},

where \eqn{\boldsymbol{y}_{i,new}}{y[i,new]} and \eqn{\boldsymbol{y}_{i,pre}}{y[i,pre]} are vectors of previous and new observations from subject \eqn{i}{i} and \eqn{q()}{q()} is a function which provides a scalar summary of the new observations.

%

These functions are subroutines of \code{\link{index.batch}}.
%
\code{CP.ar1.se} returns the estimate of the conditional probability and its asymptotic standard error of a subject based on AR(1) model.
The evaluations for the probability is done by its subroutine \code{jCP.ar1}, which, in turn, has two subroutines \code{CP1.ar1} and \code{MCCP.ar1}. \code{CP1.ar1} computes the probability via the adaptive quadrature while \code{MCCP.ar1} computes the probability via the Monte Carlo integration. 
}

\usage{

CP.ar1.se(tpar, ypre, ynew, y2m = NULL, XM, stp, 
	  RE = "G", V, MC = FALSE, qfun = "sum",i.tol=1E-75)

jCP.ar1(tpar, ypre, ynew, y2m = NULL, XM, stp, RE = "G", LG =FALSE, 
	      MC = FALSE, N.MC = 40000, qfun = "sum", oth =NULL,i.tol=1E-75)

CP1.ar1(ypre, ynew, y2m=NULL, stp, u, th, a, dt, RE = "G", oth,qfun,i.tol=1E-75)

MCCP.ar1(ypre, ynew, stp, u, th, a, dt, RE = "G", N.MC = 1000, oth, qfun = "sum") 

}

%- maybe also 'usage' for other objects documented here.


\arguments{

\item{tpar}{
A vector of length 4 + # covariates, containing the estimates of the model in the order that 
\eqn{\log(\alpha),\log(\theta),logit(\delta),\beta_0,\beta_1,\cdots}{\log(\alpha),\log(\theta),logit(\delta),\beta[0],\beta[1],...}.
%If random effects are assumed to be from the log-normal distribution, then \code{theta} is a variance estimate of the random effect.
%If random effects are assumed to be from the gamma distribution, then \code{theta} is a \code{shape} parameter estimate.
If the semi-parametric approach is taken, then \code{theta} is a place holder and can be any value.
}

\item{ypre}{
A vector of the length the number of previous observations, containing counts on pre-scans. 
}

\item{y2m}{
Internal use only. Set as \code{y2m=NULL}.
}

\item{ynew}{
A vector of length the number of new observations, containing counts on new scans.
}

\item{XM}{
See \code{\link{CP.se}}.
}

\item{stp}{
A vector of length \eqn{n_i}{n[i]}, containing index to indicates missing scans.
The first entry must be zero.
For example, if there is no missing scans and there are five repeated measures, then \code{stp=c(0,1,1,1,1)}.
If the third scan is missing and there are four repeated measures, then \code{stp=c(0,1,2,1)}.
}

\item{RE}{
See \code{\link{lmeNB}}. Note that this option is NOT accepted in \code{CP.ar1.se}.
}

\item{LG}{
See \code{\link{CP.se}}.
}

\item{MC}{
If \code{TRUE} then the function \code{MCCP.ar1} is called and the Monte carlo integration is performed to integrate out the random effect. Fast but could be unreliable; not recommended for computing the confidence intervals.
If \code{FALSE} then the function \code{CP1.ar1} is called and the adaptive quadrature is performed. Slow but reliable.
}

\item{N.MC}{
The number of the Monte Carlo integration. Necessary if \code{MC=TRUE}.  
}

\item{qfun}{
 See \code{\link{index.batch}}.
}

\item{oth}{
See \code{\link{CP.se}}.
If \code{RE="NoN"}, \code{othr} must be the frequency table of the random effects, which can be obtained based on \code{dist=obj$gi} where \code{obj} is the output of \code{\link{fitSemiAR1}}.
}


\item{V}{
See \code{\link{CP.se}}.
} 

\item{th}{The estimated \emph{theta}.}

\item{a}{The estimated \emph{alpha}.} 

\item{dt}{The estimated \emph{delta}.}



\item{u}{
A vector of length the number of repeated measures, containing the estimated mean counts 
\eqn{( \mu_{i1},\cdots,\mu_{i n_i} )}{( \mu[i1],...,\mu[i n[i] ])}.
If the mean of \emph{Y_ij} is modeled linearly on \code{beta} with the log-link function, 
then \code{u=exp( beta0+ XM[,1]*beta1+XM[,2]*beta2+...))}. 
}

\item{i.tol}{See \code{\link{lmeNB}}}
}

\references{
Detection of unusual increases in MRI lesion counts in individual multiple sclerosis patients. (2013) Zhao, Y., Li, D.K.B., Petkau, A.J., Riddehough, A., Traboulsee, A., Journal of the American Statistical Association.
}
\author{
 Zhao, Y. and Kondo, Y.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The main function to fit the Negative Binomial mixed-effect model:
\code{\link{lmeNB}},

The internal functions of \code{\link{lmeNB}} for fitting relevant models:
\code{\link{fitParaIND}},
\code{\link{fitParaAR1}},
\code{\link{fitSemiIND}},
\code{\link{fitSemiAR1}},

The other subroutines of \code{\link{index.batch}} to compute the conditional probability index:
%\code{\link{jCP.ar1}},
%\code{\link{CP1.ar1}},
%\code{\link{MCCP.ar1}},
%\code{\link{CP.ar1.se}},
\code{\link{CP.se}},
\code{\link{jCP}},

The functions to generate simulated datasets:
\code{\link{rNBME.R}}.
}

\examples{

\dontrun{
ilgt <- function (x) 
{
    tem = exp(x)
    res = tem/(1 + tem)
    return(res)
}
lgt <- function (p) 
{
    log(p/(1 - p))
}
## the vector of a parameter estimates if log(a),log(theta),logit(delta),beta0.
tpar  <- c(log(2),log(0.5),lgt(0.5),2)
ypre <- c(0, 1)
ynew <- c(1, 0, 0)
## No covariate
XM <- NULL
## no missing visit
stp <- c(0,1,1,1,1)
RE <- "G"
## The estimate of the variance covariance matrix
V <-
matrix(
c( 0.17720309, -0.240418504,  0.093562548,  0.009141980,
  -0.24041850,  0.605132808, -0.160454773, -0.003978118,
   0.09356255, -0.160454773,  0.095101658,  0.005661923,
   0.00914198, -0.003978118,  0.005661923,  0.007574769),
nrow=4)

## the estimate of the conditional probability based on the sum summary statistics and its SE
CP.ar1.se(tpar = tpar, ypre = ypre, ynew = ynew, 
	  XM =XM, stp = stp, 
	  RE = RE, V = V, MC = FALSE, qfun = "sum")

## the estimate of the conditional probability based on the max summary statistics and its SE
CP.ar1.se(tpar = tpar, ypre = ypre, ynew = ynew, 
	  XM =XM, stp = stp, 
	  RE = RE, V = V, MC = FALSE, qfun = "max")


## CP.ar1.se calls for jCP.ar1 to compute the estimate of the conditional probability
## the estimate of the conditional probability based on the sum summary statistics
jCP.ar1(tpar = tpar, ypre = ypre, ynew = ynew,
	y2m=NULL,  XM =XM, stp = stp,
        RE = RE, LG = FALSE, MC = FALSE, N.MC = 40000, qfun = "sum", oth = NULL)

## jCP.ar1 calls for CP.ar1 to compute the estimate of the conditional probability 
## via the adaptive quadrature (MC=F)
## the estimate of the conditional probability

u <- rep(exp(tpar[4]),length(ypre)+length(ynew))

CP1.ar1(ypre = ypre, ynew =ynew, 
	stp =stp, u = u, th = exp(tpar[2]), a = exp(tpar[1]), 
	dt= ilgt(tpar[3]), RE = RE, qfun = "sum")


## jCP.ar1 calls for CP.ar1 to compute the estimate of the conditional probability 
## via the Monte Carlo method (MC=T)
## the estimate of the conditional probability
MCCP.ar1(ypre = ypre, ynew =ynew, stp = stp, 
	 u = u, th = exp(tpar[2]), a = exp(tpar[1]),  dt = ilgt(tpar[3]), 
	 RE = RE, N.MC = 1000, qfun = "sum")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
