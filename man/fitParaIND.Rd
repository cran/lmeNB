\name{fitParaIND}
\alias{fitParaIND}
\alias{formulaToDat}
%- Also NEED an '\alias' for EACH other topic documented here.
%
\title{
%
Performs the maximum likelihood estimation for the negative binomial mixed-effect independent model 
%
}

\description{
This function fits the parametric negative binomial mixed-effect independent model in the formulation described Zhao et al (2013).
The conditional distribution of response count given random effect is modelled by Negative Binomial as described in description of \code{\link{lmeNB}}. 
The conditional dependence among the response counts of a subject is assumed independent. 
The random effects are modelled with either gamma or log-normal distributions. 
See descriptions of \code{\link{lmeNB}}.
}


\usage{
fitParaIND(formula, data, ID, p.ini = NULL, IPRT = FALSE, RE = "G", 
           i.tol = 1e-75, o.tol = 0.001, COV = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
See \code{\link{lmeNB}}.
}
\item{data}{
See \code{\link{lmeNB}}.
}
%

\item{ID}{
See \code{\link{lmeNB}}.
}
%
  \item{p.ini}{
The initial values of the parameters 
\eqn{\log(\alpha), \log(\theta),\beta_0, \beta_1}{\log(\alpha), \log(\theta),\beta[0], \beta[1]}
\code{NULL} is accepted.
}
  \item{IPRT}{
See \code{\link{lmeNB}}.            
}
  \item{RE}{
The distribution of random effects \eqn{G_i}{G[i]}.
If \code{RE="G"} then the random effects are assumed to be from the gamma distribution.
If \code{RE="N"} then they are assumed to be form the log-normal distribution.
}

\item{i.tol}{
See \code{\link{lmeNB}}.   
%A real number to determine the absolute tolerance for \code{\link{integrate}}.
}
\item{o.tol}{
See \code{\link{lmeNB}}.   
%A real number to determine the relative tolerance for \code{\link{optim}}.
}

\item{COV}{
Internal use only.
}
}
\details{
%
\code{fitParaIND} calls \code{optim} to minimize the negative log-likelihood of 
the negative binomial model with respect to the model parameters: 
\eqn{\log(\alpha), \log(\theta),\beta_0, \beta_1,\cdots}{\log(\alpha), \log(\theta),\beta[0], \beta[1],...}. 
The Nelder-Mead algorithm is employed.
The log-likelihood is obtained by marginalizing out the random effects. 
The numerical integration is carried out using adaptive quadrature.

The missing count responses, if assumed to be missing at random, can be igored.Other types of missing data are currently not accepted.
%




All the computations are done in \code{R}.
}

\value{

\item{call}{ See \code{\link{lmeNB}}. }
%%\item{p.ini}{}   
\item{opt}{ See \code{\link{lmeNB}}.} 
%%\item{formula}{} 
\item{nlk}{ See \code{\link{lmeNB}}.}
\item{V}{ See \code{\link{lmeNB}}.}
\item{est}{ See \code{\link{lmeNB}}.}
%% \item{RE}{}
\item{AR}{\code{FALSE}}

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
  The wrapper function for all the negative binomial mixed effect
  regression:
  \code{\link{lmeNB}}.
  
  The functions to fit the other negative binomial mixed effect models:
  %\code{\link{fitParaIND}},
  \code{\link{fitParaAR1}},
  \code{\link{fitSemiIND}},
  \code{\link{fitSemiAR1}},
  
  The subroutines of \code{\link{index.batch}} to compute the conditional probability index:
  \code{\link{jCP.ar1}},
  \code{\link{CP1.ar1}},
  \code{\link{MCCP.ar1}},
  \code{\link{CP.ar1.se}},
  \code{\link{CP.se}},
  \code{\link{jCP}},
  
  The functions to generate simulated datasets:
  \code{\link{rNBME.R}}.
}

\examples{


\dontrun{

## ================================================================================ ##
## generate a simulated dataset from the negative binomial mixed-effect independent model:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i ~ Gamma(scale=th,shape=1/th)
set.seed(1)
sn <- 5 ## the number of repeated measures of each patient
n <- 80 ## the number of patients
loga <- - 0.5
a <- exp(loga) ## the parameter for the failure probability of the negative binomial distribution
logtheta <- 1.3
th <- exp(logtheta) ## the parameter for the gamma distribution of g_i

## No difference between the means of groups 
## The model only has an intercept term beta0 = 0.5
b0 <- 0.5
u1 <- rep(exp(b0),sn) ## the mean vector for group 1 at time point 1,...,sn
u2 <- rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn


DT0 <- rNBME.R(gdist="G", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn)
ID <- DT0$id
Vcode <- rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
new <- Vcode > 0 
dt1 <- data.frame(CEL=DT0$y)
logitd <- -0.5

## ================================================================================ ##


## [1]: Fit the negative binomial independent model 
## where the random effects are from the gamma distribution. This is the true model.

re.gamma.ind <- fitParaIND(formula=CEL~1,data=dt1,ID=ID,RE="G", 
	                  p.ini=c(loga,logtheta,b0),IPRT=TRUE)
## compute the estimates of the conditional probabilities 
## with sum of the new repeated measure as a summary statistics 
Psum <- index.batch(olmeNB=re.gamma.ind, ID=ID,data=dt1,
                    labelnp=new,qfun="sum", IPRT=TRUE) 

## compute the estimates of the conditional probabilities 
## with max of the new repeated measure as a summary statistics 
Pmax <- index.batch(olmeNB=re.gamma.ind, ID=ID,data=dt1,
                   labelnp=new,qfun="max", IPRT=TRUE)

## Which patient's estimated probabilities 
## based on the sum and max statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax$condProbSummary[,1]-Psum$condProbSummary[,1]))==80) )
## Show the patient's CEL counts  
dt1$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum$condProbSummary[IDBigDif,]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax$condProbSummary[IDBigDif,]



## [2]: Fit the negative binomial independent model 
## where the random effects are from the lognormal distribution. 
re.logn.ind <- fitParaIND(formula=CEL~1,data=dt1,ID=ID, 
                          RE="N", 			   	
		          p.ini=c(loga,logtheta,b0), 		
		          IPRT=TRUE)

Psum <- index.batch(olmeNB=re.logn.ind, ID=ID,data=dt1,C=TRUE,
                   labelnp=new,qfun="sum", IPRT=TRUE) 

## [3]: Fit the semi-parametric negative binomial independent model 


re.semi.ind <- fitSemiIND(formula=CEL~1,data=dt1,ID=ID)

Psum <- index.batch(olmeNB=re.semi.ind,ID=ID,data=dt1, i.se = FALSE,
                   labelnp=new, qfun="sum", IPRT=TRUE) 



## [4]: Fit the negative binomial mixed-effect AR(1) model 
## where random effects are from the gamma distribution


re.gamma.ar1 <- fitParaAR1(formula=CEL~1,data=dt1,ID=ID,
		           p.ini=c(loga,logtheta,logitd,b0), 
		           RE="G", Vcode=Vcode,
		           IPRT=TRUE)

Psum <- index.batch(olmeNB=re.gamma.ar1, ID=ID,data=dt1, labelnp=new,Vcode=Vcode,
                    qfun="sum", IPRT=TRUE,i.se=FALSE) ## i.se=TRUE requires more time...
	 




## ======================================================================== ##
## == Which model performed the best in terms of the estimation of beta0 == ##
## ======================================================================== ##

getpoints <- function(y,estb0,sdb0=NULL,crit=qnorm(0.975))
{	  
points(estb0,y,col="blue",pch=16)
if (!is.null(sdb0))
{
points(c(estb0-crit*sdb0,estb0+crit*sdb0),rep(y,2),col="red",type="l")
}
}

ordermethod <- c("gamma.ind","logn.ind","semi.ind","gamma.ar1")

estb0s <- c(
re.gamma.ind$est[3,1], 
re.logn.ind$est[3,1],
re.semi.ind$est[3],
re.gamma.ar1$est[4,1]
)

## The true beta0 is:
b0
c <- 1.1
plot(0,0,type="n",xlim=c(min(estb0s)-0.5,max(estb0s)*c),
	ylim=c(0,5),yaxt="n",
main="Simulated from the independent model \n with random effect ~ gamma")

legend("topright",
	col="red",
	legend=ordermethod)
abline(v=b0,lty=3)

## [1] gamma.ind
sdb0 <- re.gamma.ind$est[3,2]
getpoints(4,estb0s[1],sdb0)

## [2] logn.ind
sdb0 <- re.logn.ind$est[3,2]
getpoints(3,estb0s[2],sdb0)

## [3] semi.ind
getpoints(2,estb0s[3])

## [4] gamma.ar1
sdb0 <- re.gamma.ar1$est[4,2]
getpoints(1,estb0s[4],sdb0)

}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
