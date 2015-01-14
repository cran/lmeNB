\name{fitSemiIND}
\alias{fitSemiIND}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit the semi-parametric negative binomial mixed-effect independent model.
}

\description{
This function fits the semi-parametric negative binomial mixed-effect independent model to repeated count responses (Zhao et al. 2013). 
The conditional distribution of response count given random effect is modelled by Negative Binomial as described in description of \code{\link{lmeNB}}. 
The conditional dependence among the response counts of a subject is assumed independent. 
The semiparametric procedure is employed for random effects.
See descriptions of \code{\link{lmeNB}}.
}
\usage{
fitSemiIND(formula,data,ID, p.ini = NULL, IPRT = TRUE, deps = 1e-04, 
	   maxit = 100, u.low = 0)
}
\arguments{
\item{formula}{
See \code{\link{lmeNB}}.
}
\item{data}{
See \code{\link{lmeNB}}.
}
\item{ID}{
See \code{\link{lmeNB}}.
}
\item{p.ini}{
A vector of length 3 + # covariates, containing the initial values for
the parameters \eqn{(\log(\alpha), \log(Var(G_i)), \beta_0, \beta_1, ...)}{(\log(alpha), \log(Var(G[i])), beta[0], beta[1],...)}. 
\code{NULL} is accepted.
}
\item{IPRT}{
See \code{\link{lmeNB}}.
}
\item{deps}{
See \code{\link{lmeNB}}.
}
\item{maxit}{
See \code{\link{lmeNB}}.
}	
\item{u.low}{
See \code{\link{lmeNB}}.
}

}
\details{
The algorithm repeats the following four steps until a stoping criterion is satisfied:

Step 1) Given \eqn{\alpha}{\alpha}, Estimate the coefficients of covariates by the method of generalized Least Squares.

That is, this step solves for:
\eqn{
argmin_{\boldsymbol{\beta}} \sum_{i=1}^{N} 
(\boldsymbol{Y}_{i}-E(\boldsymbol{Y}_{i};\boldsymbol{\beta}))^T \boldsymbol{W}_{i} (\boldsymbol{Y}_{i}-E(\boldsymbol{Y}_{i};\boldsymbol{\beta}))
}{
argmin_{\beta}   \sum[i=1]^N (Y[i]-E(Y[i];\beta))^T W[i] (Y[i]-E(Y[i];\beta))
}
where the weight matrix for each patient \eqn{\boldsymbol{W}_{i}}{W[i]} is selected to \eqn{Var(\boldsymbol{Y}_{i})^{-1} }{ Var(Y[i])^{-1}  } 
(which is a function of \eqn{\alpha}{\alpha})
if it exists, else it is set to be an identity matrix.


Step 2) Approximate the distribution of the random effect \eqn{G_i}{G[i]} by \eqn{\gamma}{\gamma}.

Step 3) Estimate \eqn{\alpha}{\alpha} by minimizing the negative psudo-profile likelihood. 
The numerical minimization is carried out using \code{optimize} and the numerical integration is carried out using adaptive quadrature.

Step 4) Estimate \eqn{Var(G_i)}{Var(G[i])} by the medhod of moment and update the weights.

All the computations are done in \code{R}.
}

\value{


\item{opt}{
See \code{\link{lmeNB}}.
}
\item{diffPara}{
The largest absolute difference of parameter vectors between the current and previous iterations.
}
\item{V}{
\code{NULL}
}
\item{est}{
See \code{\link{lmeNB}}.
}
\item{gtb}{
The relative frequency table of \eqn{G_i}{G[i]}, (\eqn{i=1,\cdots,N}{i=1,...,N}).
\code{gh1} (the second column) contains the unique values in ascending order and \code{ghw} (the first column) contains the corresponding relative frequencies.
}

\item{counter}{The number of iterations before the algorithm was terminated}

\item{gi}{
A vector of length \eqn{N}{N}, containing the approximated random effect \eqn{G_i,i=1,\cdots,N}{G[i],i=1,...,N}.
}
\item{RE}{
\code{"NoN"}, denoting that the fitted model is a semi-parametric mixed-effect model.
}
\item{AR}{
\code{FALSE}
}
\item{paraAll}{
Record estimated parameters at every iteration.
}

%\item{weightAll}{
%Record the weighted matrix for each patient at every iteration.
%}




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

The functions to fit the other models:
\code{\link{fitParaIND}},
\code{\link{fitParaAR1}},
%\code{\link{fitSemiIND}},
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
## generate a simulated dataset from the negative binomial 
## mixed-effect independent model:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i is from unknown distribution
## For the simulation purpose, G_i's are from the mixture of 
## the gamma and the log-normal distributions.


sn <- 5 ## the number of repeated measures of each subject
n <- 80 ## the number of subjects
logtheta <- 1.3
th <- exp(logtheta) ## the parameter for the gamma distribution of g_i
loga <- -0.5
## the parameter for the failure probability of the negative binomial distribution
a <- exp(loga) 
b0 <- 0.5
u1 <- rep(exp(b0),sn) ## the mean vector for group 1 at time point 1,...,sn
u2 <- rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn

DT3 <- rNBME.R(gdist="GN", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn,
	       othrp=list(p.mx=0.1,u.n=3,s.n=1,sh.mx = NA) ## 0  < p.mx < 1
	       )

ID <- DT3$id
dt3 <- data.frame(CEL=DT3$y)

Vcode <- rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
new <-Vcode>0         # new scans: 1,2,3


## 1) Fit the negative binomial mixed-effect AR(1) model 
## where random effects is from the gamma distribution

logitd <- -0.2
re.gamma.ar1 <- fitParaAR1(formula=CEL~1,data=dt3,ID=ID,
		          Vcode=Vcode, 
		          p.ini=c(loga,logtheta,logitd,b0), 
		          ## log(a), log(theta), logit(d), b0
		          RE="G", 
		          IPRT=TRUE)


## compute the estimates of the conditional probabilities 
## with sum of the new repeated measure as a summary statistics 
## i.se=FALSE,C=TRUE options for speed up!
Psum <- index.batch(olmeNB=re.gamma.ar1,data=dt3,ID=ID, Vcode=Vcode, 
                    labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE,C=TRUE,i.tol=1E-3)  
		 


## 2) Fit the negative binomial mixed-effect AR(1) model 
## where random effects is from the log-normal distribution


re.logn.ar1 <- fitParaAR1(formula=CEL~1,data=dt3,ID=ID,
		          Vcode=Vcode, RE="N", IPRT=TRUE)

## REQUIRES SOME TIME..
Psum <- index.batch(olmeNB=re.logn.ar1, data=dt3,ID=ID,Vcode=Vcode,
	            labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE,C=TRUE,i.tol=1E-3) 



## 3) Fit the negative binomial independent model 
## where random effects is from the gamma distribution
re.gamma.ind <- fitParaIND(formula=CEL~1,data=dt3,ID=ID, 
                           RE="G",IPRT=TRUE)

Psum <- index.batch(olmeNB=re.gamma.ind, data=dt3,ID=ID,Vcode=Vcode,
                    labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE) 



## 4) Fit the negative binomial independent model 
## where random effects is from the lognormal distribution
re.logn.ind <- fitParaIND(formula=CEL~1,data=dt3,ID=ID, 
                          RE="N", 			   	
		          p.ini=c(loga,logtheta,b0), 		
		          IPRT=TRUE)

Psum <- index.batch(olmeNB=re.logn.ind,data=dt3,ID=ID,labelnp=new,qfun="sum", IPRT=TRUE) 

## 5) Fit the semi-parametric negative binomial AR(1) model 

logvarG <- -0.4

re.semi.ar1 <- fitSemiAR1(formula=CEL~1,data=dt3,ID=ID,Vcode=Vcode)
Psum <- index.batch(olmeNB=re.semi.ar1,data=dt3,ID=ID,Vcode=Vcode, 
                    labelnp=new,qfun="sum", IPRT=TRUE,MC=TRUE,i.se=FALSE) 


## 6) Fit the semi-parametric negative binomial independent model 
## This is closest to the true model
re.semi.ind <- fitSemiIND(formula=CEL~1,data=dt3,ID=ID, p.ini=c(loga, logvarG, b0))

## compute the estimates of the conditional probabilities 
## with sum of the new repeated measure as a summary statistics 
Psum <- index.batch(olmeNB=re.semi.ind,data=dt3,ID=ID, labelnp=new, 
                    qfun="sum", IPRT=TRUE,i.se=FALSE) 
## compute the estimates of the conditional probabilities 
## with max of the new repeated measure as a summary statistics 
Pmax <- index.batch(olmeNB=re.semi.ind, data=dt3,ID=ID,labelnp=new, qfun="max", 
                    IPRT=TRUE,i.se=FALSE) 


## Which patient's estimated probabilities based on the sum and max 
## statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax$condProbSummary[,1]-Psum$condProbSummary[,1]))==80) )
## Show the patient's CEL counts  
dt3$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum$condProbSummary[IDBigDif,]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax$condProbSummary[IDBigDif,]


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
ordermethod <- c("gamma.ar1","logn.ar1","gamma.ind","logn.ind","semi.ar1","semi.ind")

estb0s <- c(
re.gamma.ar1$est[4,1],
re.logn.ar1$est[4,1],
re.gamma.ind$est[3,1],
re.logn.ind$est[3,1],
re.semi.ar1$est[4],
re.semi.ind$est[3]
)

## The true beta0 is:
b0
c <- 1.1
plot(0,0,type="n",xlim=c(min(estb0s)-0.5,max(estb0s)*c),ylim=c(0,7),yaxt="n",
main="Simulated from the independent model \n with random effect ~ mixture of normal and gamma")

legend("topright",
	legend=ordermethod)
abline(v=b0,lty=3)

## 1) gamma.ar1
sdb0 <- re.gamma.ar1$est[4,2]
getpoints(6,estb0s[1],sdb0)

## 2)logn.ar1
sdb0 <- re.logn.ar1$est[4,2]
getpoints(5,estb0s[2],sdb0)

## 3) gamma.ind
sdb0 <- re.gamma.ind$est[3,2]
getpoints(4,estb0s[3],sdb0)

## 4) logn.ind
sdb0 <- re.logn.ind$est[3,2]
getpoints(3,estb0s[4],sdb0)

## 5) semi.ar1
getpoints(2,estb0s[5])

## 6) semi.ind
getpoints(1,estb0s[6])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
