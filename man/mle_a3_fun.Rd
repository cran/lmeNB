\name{mle.a3.fun}
\alias{mle.a3.fun}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit the semi-parametric negative binomial mixed-effect independent model.
}

\description{

%
This function fits the semi-parametric negative binomial mixed-effect independent model to repeated count responses (Zhao et al).
%
The negative binomial mixed-effect independent model assumes that given the random effect \emph{ G_i=g_i}, the count responses \emph{Y_ij}s of subject \emph{i}, \emph{(i=1, ..., N)}, at time points \emph{j=1,...,n_i} are independent andfollow the negative binomial distribution:
%

                   \emph{ Y_ij | G_i =g_i ~ NB(r_ij,p_i),} 

%
where  \emph{ p_i }, the failure probability of subject \emph{i} at each time point \emph{j} is parametrized as: 

                    \emph{ p_i= 1/(g_i*a+1),}  

and \emph{ a > 0 }.
% 
The model assumes \emph{ E(G_i) = 1 } so that 
%
 \emph{ E(Y_ij|G_i=g_i)=r_ij*g_i*a } and \emph{ E(Y_ij)= r_ij*g_i }.
%
This assumption allows the interpretation of the latent random variable \emph{ G_i } as the subject i's activity rate relative to the overall cohort.
%

%
The marginal mean \emph{ mu_ij = E(Y_ij)} is modeled with fixed effect coefficients, \emph{beta}: \emph{mu_ij= exp(X_ij^T beta)}.

The distribution of random effect terms are not assumed and approximated by the estimated values of the quantity:

\emph{
gamma_i = w_i (y_i+/ mu_i+) + (1-w_i) , i=1,...,N
},

where 
\emph{
y_i+ = sum_j=1^n_i y_ij 
}, 
\emph{
mu_i+ =  sum_j=1^n_i mu_ij 
} and, 
\emph{
w_i = sqrt( Var(G_i)/Var(Y_i+/mu_i+) )
}.
See Zhao et al. for more details.
%

The missing count responses, if assumed to be missing at random, can be igored. Other types of missing data are currently not accepted. 
}


\usage{
mle.a3.fun(formula,data,ID, p.ini = NULL, IPRT = TRUE, deps = 1e-04,maxit = 100)
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
Each row must contain the data corresponding to the repeated measure j of a subject and the rows \emph{(i,j)}s must be ordered as 
\emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}. Missing values are not accepted.
}
%

\item{ID}{
A vector of length \emph{ sum_i^N n_i }, containing the subject IDs. i.e., \code{c(rep(ID_1,n_1),rep(ID_2,n_2),...,rep(ID_N,n_N))}.
}

  \item{p.ini}{
A vector of length 3 + # covariates, containing the initial values for
the parameters \emph{(log(a), log(var(G)), beta0, beta1, ...)}.
    \code{NULL} is accepted.
}
  \item{IPRT}{
A logical, passed to Iprt of function \code{optim}. 
If \code{TRUE} then print iterations. 
}
  \item{deps}{
A small positive real number for the stopping criterion: the algorithm stops iteration when \code{max(p.new - p.old) < deps}
}
\item{maxit}{
The maximum number of iterations.
}
}

\details{
The algorithm repeats the following four steps until a stoping criterion is satisfied:

Step 1) Estimate the coefficients of covariates by the method of generalized Least Squares.

Step 2) Approximate the distribution of the random effect \emph{G_i} by \emph{gamma_i}.

Step 3) Estimate \emph{a} by minimizing the negative psudo-profile likelihood. The numerical minimization is carried out using \code{optimize} and the numerical integration is carried out using adaptive quadrature.

Step 4) Estimate \emph{var(G)} by the medhod of moment and update the weights.


}

\value{

\item{opt}{
The values returned by \code{optimize} to minimize the negative of the psudo-profile log-likelihood with respect to \code{a}.
}

\item{V}{
If the number of covariate is nonzero, \code{vcm} returns a naive estimate of the variance covariance matrix of \code{b0}, \code{b1},...
}

\item{est}{
A vector of length # covariates + 2 containing the estimates of the parameters \emph{(log(a), log(var(G)), b0, b1, ...)}
}

\item{gi}{
A vector of length \emph{N}, containing the approximation of \emph{g_i}s
}

\item{gtb}{
The relative frequency table of \emph{g_i}, (\emph{i=1,...,N}).
\code{gh1} contains the unique values in ascending order and \code{ghw} contains relative frequencies.
}



\item{mod}{
\code{"NoN"}, denoting that the fitted model is a semi-parametric mixed-effect model.
}

% \item{idat}{
% A dataframe, containing ID, CEL, x.1, x.2, ...  
% %Vcode contains the scan number where -1,0,...,ni correspond to screening, baseline, and follow-up scans 1,...,ni.
% The column labeled as CEL contains the response counts.
% }


\item{cor}{\code{"ind"}, indicating that the model assumes independent structure of the count responses given the random effects.}

\item{maxit}{The maximum number of iterations.}
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
The functions to fit the other models:
\code{\link{mle.fun}},
\code{\link{mle.ar1.fun}},
%\code{\link{mle.a3.fun}},
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
## generate a simulated dataset from the negative binomial mixed-effect independent model:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i is from unknown distribution
## For the simulation purpose, G_i's are from the mixture of the gamma and the log-normal distributions.


sn = 5 ## the number of repeated measures of each subject
n = 80 ## the number of subjects
logtheta = 1.3
th = exp(logtheta) ## the parameter for the gamma distribution of g_i
loga = -0.5
a = exp(loga) ## the parameter for the failure probability of the negative binomial distribution
b0 = 0.5
u1 = rep(exp(b0),sn) ## the mean vector for group 1 at time point 1,...,sn
u2 = rep(exp(b0),sn) ## the mean vector for group 2 at time point 1,...,sn

DT3= rNBME.R(gdist="GN", n=n, a=a, th=th, u1=u1, u2=u2, sn=sn,
	    othrp=list(p.mx=0.1,u.n=3,s.n=1,sh.mx = NA) ## 0  < p.mx < 1
	    )

ID = DT3$id
dt3=data.frame(CEL=DT3$y)

Vcode=rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
new=Vcode>0         # new scans: 1,2,3


## 1) Fit the negative binomial mixed-effect AR(1) model where random effects is from the gamma distribution

logitd = -0.2
re.gamma.ar1=mle.ar1.fun(formula=CEL~1,data=dt3,ID=ID,
		      Vcode=Vcode, 
		      p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
		       model="G", 
		      IPRT=TRUE)


## compute the estimates of the conditional probabilities with sum of the new repeated measure as a summary statistics 
Psum=index.batch(olmeNB=re.gamma.ar1,data=dt3,ID=ID, Vcode=Vcode, labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE)  
		 


## 2) Fit the negative binomial mixed-effect AR(1) model where random effects is from the log-normal distribution


re.logn.ar1=mle.ar1.fun(formula=CEL~1,data=dt3,ID=ID,
		      Vcode=Vcode, model="N", IPRT=TRUE)
\dontrun{
## REQUIRES SOME TIME..
Psum=index.batch(olmeNB=re.logn.ar1, data=dt3,ID=ID,Vcode=Vcode,
	         labelnp=new,qfun="sum", IPRT=TRUE)#,i.se=FALSE) 
}


## 3) Fit the negative binomial independent model where random effects is from the gamma distribution
re.gamma.ind=mle.fun(formula=CEL~1,data=dt3,ID=ID, 
                   model="G", 
		   IPRT=TRUE)

Psum=index.batch(olmeNB=re.gamma.ind, data=dt3,ID=ID,Vcode=Vcode,
                 labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE) 



## 4) Fit the negative binomial independent model where random effects is from the lognormal distribution
re.logn.ind=mle.fun(formula=CEL~1,data=dt3,ID=ID, 
                   model="N", 			   	
		   p.ini=c(loga,logtheta,b0), 		
		   IPRT=TRUE)

Psum=index.batch(olmeNB=re.logn.ind,data=dt3,ID=ID,labelnp=new,qfun="sum", IPRT=TRUE) 

## 5) Fit the semi-parametric negative binomial AR(1) model 

logvarG = -0.4

re.semi.ar1=mle.ar1.non3(formula=CEL~1,data=dt3,ID=ID,Vcode=Vcode)
Psum=index.batch(olmeNB=re.semi.ar1,data=dt3,ID=ID,Vcode=Vcode, labelnp=new,qfun="sum", IPRT=TRUE,iMC=TRUE,i.se=FALSE) 


## 6) Fit the semi-parametric negative binomial independent model 
## This is closest to the true model
re.semi.ind=mle.a3.fun(formula=CEL~1,data=dt3,ID=ID, p.ini=c(loga, logvarG, b0))

## compute the estimates of the conditional probabilities with sum of the new repeated measure as a summary statistics 
Psum=index.batch(olmeNB=re.semi.ind,data=dt3,ID=ID, labelnp=new, qfun="sum", IPRT=TRUE,i.se=FALSE) 
## compute the estimates of the conditional probabilities with max of the new repeated measure as a summary statistics 
Pmax=index.batch(olmeNB=re.semi.ind, data=dt3,ID=ID,labelnp=new, qfun="max", IPRT=TRUE,i.se=FALSE) 


## Which patient's estimated probabilities based on the sum and max statistics disagrees the most?
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
