\name{mle.fun}
\alias{mle.fun}
\alias{formulaToDat}
%- Also NEED an '\alias' for EACH other topic documented here.
%
\title{
%
Performs the maximum likelihood estimation for the negative binomial mixed-effect independent model 
%
}

\description{
%
This function fits a negative binomial mixed-effect independent model to repeated count measures (Zhao et al).
%
The model assumes that given the random effect \emph{ G_i=g_i}, the count responses \emph{ Y_ij}s of subject \emph{i}, (\emph{i = 1, ..., N}), at time points \emph{j =1,...,n_i)}
 are independent and follow  the negative binomial distribution:
%

                   \emph{ Y_ij | G_i =g_i ~ NB(r_ij,p_i),} 

%
where \emph{ p_i }, the failure probability  of subject \emph{i} at every time point \emph{j} is parametrized as: 

                    \emph{ p_i= 1/(g_i*a+1),}  

and \emph{ a > 0 }.
% 
The model assumes \emph{ E(G_i) = 1 } 
%
so that \emph{ E(Y_ij|G_i=g_i)=r_ij*g_i*a } and \emph{ E(Y_ij)= r_ij*g_i }.
%
This assumption allows the interpretation of the latent random variable \emph{ G_i } as the subject i's activity rate relative to the overall cohort.
%
The random effect terms \emph{G_i}s are assumed to be either from the log-normal distribution with \emph{E(G_i)=1} and \emph{var(G_i)=theta} or the gamma distribution with \code{scale} = \emph{theta} and \code{shape}=1/\emph{theta}.
%
 
%
The marginal mean \emph{ mu_ij=E(Y_ij) } is modeled with fixed effect coefficients, \emph{beta}: \emph{ mu_ij= exp(X_ij^T beta)}.
%
%


}


\usage{
mle.fun(formula, data, ID, p.ini=NULL, IPRT = FALSE, model = "G",i.tol=1e-75,o.tol=1.e-3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%
An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
%
The formula must contain an intercept term.
%
}
\item{data}{
A data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model.
Each row must contain the data corresponding to the repeated measure j of a subject  and the rows \emph{(i,j)}s must be ordered as 
\emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}. Missing values are not accepted.
}
%

\item{ID}{
A vector of length \emph{ sum_i^N n_i }, containing the patient IDs. i.e., \code{c(rep(ID_1,n_1),rep(ID_2,n_2),...,rep(ID_N,n_N))}.
}
%
  \item{p.ini}{
The initial values of the parameters (log(\emph{a}), log(\emph{theta}),
\emph{beta0}, \emph{beta1}, ...).
\code{NULL} is accepted.

}
  \item{IPRT}{
A logical, passed to Iprt of function \code{optim}. 
If \code{TRUE} then print iterations. 
                
}
  \item{model}{
The distribution of random effects \emph{G_i}.
 If \code{model="G"} then the random effects are assumed to be from the gamma distribution.
If \emph{model="N"} then they are assumed to be form the log-normal distribution.
}

\item{i.tol}{
A real number to determine the tolerance for \code{\link{integrate}}.
}
\item{o.tol}{
A real number to determine the tolerance for \code{\link{optim}}.
}

}
\details{
%
\code{mle.fun} calls \code{optim} to minimize the negative log-likelihood of the negative binomial model with respect to the model parameters: c(log(\emph{a}), log(\emph{theta}), \emph{beta0}, \emph{beta1}, ...). 
%
The Nelder-Mead algorithm is employed.
% 
The log-likelihood is obtained by marginalizing out the random effects. 
%
The numerical integration is carried out using adaptive quadrature.
%

The missing count responses, if assumed to be missing at random, can be igored.Other types of missing data are currently not accepted.
}

\value{
\item{opt}{
The values returned by \code{optim}.
}
\item{nlk}{
The value of the negative log-likelihood corresponding to \code{opt$par} 
}
\item{V}{ 
The approximated asymptotic covariance matrix of the maximum likelihood estimators. \code{V}=\code{solve(opt$hessian)}
}
\item{est}{ 
 A (3 + # covariates) by 2 matrix. The first column contains the estimates of the model parameters, 
log(\emph{a}), log(\emph{theta}),\emph{beta0},\emph{beta1},...)
%
The second column contains the approximated standard deviations of the estimators, i.e., \code{sqrt(diag(V))}
}

\item{mod}{
If \code{model="G"} then \code{mod="G"}.
If \code{model="N"} then \code{mod="N"}.
}

% \item{idat}{
% A dataframe, containing ID, CEL, x.1, x.2, ...  
% %Vcode contains the scan number where -1,0,...,ni correspond to screening, baseline, and follow-up scans 1,...,ni.
% The column labeled as CEL contains the response counts.
% }


\item{cor}{\code{"ind"}, indicating that the model assumes independent structure of the count responses given the random effects.}
}

\references{
 Zhao, Y., Li, D.K.B., Petkau, J.A., Riddehough, A. & Traboulsee, A.  Detection of unusually large increases in 
MRI lesion counts in multiple sclerosis patients.
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
  %\code{\link{mle.fun}},
  \code{\link{mle.ar1.fun}},
  \code{\link{mle.a3.fun}},
  \code{\link{mle.ar1.non3}},
  
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
sn = 5 ## the number of repeated measures of each patient
n = 80 ## the number of patients
loga = - 0.5
a = exp(loga) ## the parameter for the failure probability of the negative binomial distribution
logtheta <- 1.3
th = exp(logtheta) ## the parameter for the gamma distribution of g_i

## No difference between the means of groups 
## The model only has an intercept term beta0 = 0.5
b0 = 0.5
u1 = rep(exp(b0),sn) ## the mean vector for group 1 at time point 1,...,sn
u2 = rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn


DT0=rNBME.R(gdist="G", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn)
ID = DT0$id
Vcode=rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
new = Vcode > 0 
dt1=data.frame(CEL=DT0$y)
logitd = -0.5

## ================================================================================ ##


## [1]: Fit the negative binomial independent model where the random effects are from the gamma distribution. This is the true model.
## This is the true model

re.gamma.ind=mle.fun(formula=CEL~1,data=dt1,ID=ID,model="G", 
	             p.ini=c(loga,logtheta,b0),IPRT=TRUE)
## compute the estimates of the conditional probabilities with sum of the new repeated measure as a summary statistics 
Psum=index.batch(olmeNB=re.gamma.ind, ID=ID,data=dt1,
                 labelnp=new,qfun="sum", IPRT=TRUE) 

## compute the estimates of the conditional probabilities with max of the new repeated measure as a summary statistics 
Pmax=index.batch(olmeNB=re.gamma.ind, ID=ID,data=dt1,
                 labelnp=new,qfun="max", IPRT=TRUE)

## Which patient's estimated probabilities based on the sum and max statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax$condProbSummary[,1]-Psum$condProbSummary[,1]))==80) )
## Show the patient's CEL counts  
dt1$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum$condProbSummary[IDBigDif,]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax$condProbSummary[IDBigDif,]



## [2]: Fit the negative binomial independent model where the random effects are from the lognormal distribution. 
re.logn.ind=mle.fun(formula=CEL~1,data=dt1,ID=ID, 
                   model="N", 			   	
		   p.ini=c(loga,logtheta,b0), 		
		   IPRT=TRUE)

Psum=index.batch(olmeNB=re.logn.ind, ID=ID,data=dt1,
                 labelnp=new,qfun="sum", IPRT=TRUE) 



## [3]: Fit the semi-parametric negative binomial independent model 


re.semi.ind=mle.a3.fun(formula=CEL~1,data=dt1,ID=ID)

Psum=index.batch(olmeNB=re.semi.ind,ID=ID,data=dt1, 
                 labelnp=new, qfun="sum", IPRT=TRUE) 



## [4]: Fit the negative binomial mixed-effect AR(1) model where random effects are from the gamma distribution


re.gamma.ar1=mle.ar1.fun(formula=CEL~1,data=dt1,ID=ID,
		      p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
		       model="G", Vcode=Vcode,
		      IPRT=TRUE)

Psum=index.batch(olmeNB=re.gamma.ar1, ID=ID,data=dt1, labelnp=new,Vcode=Vcode,
                 qfun="sum", IPRT=TRUE,i.se=FALSE) ## i.se=T requires more time...
	 




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
