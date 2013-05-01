\name{mle.ar1.fun}
\alias{mle.ar1.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs the maximum likelihood estimation for the negative binomial mixed-effect AR(1) model 
}

%
\description{
%
This function fits a negative binomial mixed-effect AR(1) model in the formulation described Zhao et al.
%
Under this model, the count response measure \emph{Y_ij} from a subject \emph{i} at time point \emph{j} (\emph{i=1,...,N,j=1,...,n_i}) is from the negative binomial distribution:
%
 
\emph{ Y_ij | G_i =g_i ~ NB(r_ij,p_i)},

%
and the parametrization of \emph{r_ij, p_i} and the distributional assumption on \emph{G_i} are the same as the independence model (\code{\link{mle.fun}}). 
%
Given emph{ G_i=g_i}, \emph{Y_ij} depends on \emph{Y_i(j-1)} through the beta binomial thinning and is conditionally independent on \emph{Y_ij'} given \emph{Y_i(j-1)} for all \emph{j' < j-1 }.
%
The beta binomial thinning operator depends on a  parameter \emph{d} which indicates the strength of the  positive AR(1) association between repeated measures of a subject: the larger \emph{d}, the stronger the positive correlations between the repeated measures of the same subject are.
%
This interpretation depends on the result:
% 
\emph{cov(Y_ij,Y_ij'|G_i=g_i)=d^{j-j'} E(Y_ij')*(a*g_i^2+g_i)}.
%
See Zhao et al., for more details.
}

\usage{
mle.ar1.fun(formula, data, ID, Vcode,
            p.ini=NULL, IPRT = FALSE, model = "G", 
            i.tol = 1e-75, o.tol = 0.001)
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
Each row must contain the data corresponding to the repeated measure j of subject and the rows \emph{(i,j)}s must be ordered as 
\emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}. Missing data are currently not accepted.
}
%

\item{ID}{
A vector of length \emph{ sum_i^N n_i }, containing the subject IDs. i.e., \code{c(rep(ID_1,n_1),rep(ID_2,n_2),...,rep(ID_N,n_N))}
}
%

\item{Vcode}{
A vector of length the total number of repeated measures, containing the indices of time point.
For example, 
there are three subjects and the first two subjects do not have missing visits and completed five visits while the last  subject missed the third visit  and have four visits in total, then
\code{Vcode=c(1,2,3,4,5,1,2,3,4,5,1,2,4,5)} 
}

  \item{p.ini}{
Initial values for the parameters
 (log(\emph{a}), log(\emph{theta}),\emph{logit(d)} \emph{beta0}, \emph{beta1}, ...)
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

%to add
\details{
\code{mle.ar1.fun}  calls \code{optim} to minimize the negative log-likelihood of the negative binomial model with respect to the model parameters: c(log(\emph{a}), log(\emph{theta}), logit(\emph{d}), \emph{beta0}, \emph{beta1}, ...).
%
The Nelder-Mead algorithm is employed.
% 
The log-likelihood is obtained by marginalizing out the random effects. 
%
The numerical integration is carried out using adaptive quadrature.

When missing visits are present, an approximation of the likelihood is used (see Web Appendix II in Zhao et al for details.)
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
The (4 + # covariates) by 2 matrix. The first column contains the estimates of the fixed effects, 
(log(\emph{a}), log(\emph{theta}),logit(\emph{d}), \emph{beta0},\emph{beta1},...)
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

%add
\item{Vcode}{A vector containing the input indices of time point.}

\item{cor}{\code{"ar1"}, indicating that the model assumes AR1 correlation structure of the count responses given the random effects.}
}

\references{
 Zhao, Y., Li, D.K.B., Petkau, A.J., Riddehough, A. & Traboulsee, A.  Detection of unusually large increases in 
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
%\code{\link{mle.ar1.fun}},
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

## ==========================================================================================
## generate a data based on the negative binomial mixed-effect AR(1) model.
## Under this model, the response counts follows the negative binomial:
## Y_ij | G_i = g_i ~ NB(r_ij,p_i) where r_ij = exp(X^T beta)/a , p_i =1/(a*g_i+1)
## with G_i ~ Gamma(scale=th,shape=1/th)
## 
## The adjacent repeated measures of the same subject are correlated with correlation structure:
## cov(Y_ij,Y_ij'|G_i=g_i)=d^{j-j'} E(Y_ij')*(a*g_i^2+g_i)  

# log(a) = -0.5, log(th)=1.3, logit(delta) = -0.2
# b0 =  0.5, no covariates; 
loga = -0.5
logtheta= 1.3
logitd = -0.2
b0 = 0.5
# 80 subjects each with 5 scans
n = 80
sn = 5


DT2 =  rNBME.R(gdist = "G",
               n = n, ## 	the total number of subjectss	       
	       sn = sn,
	       th=exp(logtheta),
               u1 = rep(exp(b0),sn),
	       u2 = rep(exp(b0),sn),
	       a = exp(loga),
	       d = exp(logitd)/(1+exp(logitd))
	      )
Vcode=rep(-1:(sn-2),n) # scan number -1, 0, 1, 2, 3
ID = DT2$id
new = Vcode > 0
dt2=data.frame(CEL=DT2$y)

## ==========================================================================================

## 1) Fit the negative binomial mixed-effect AR(1) model where the random effects are from the gamma distribution
## This is the true model

re.gamma.ar1=mle.ar1.fun(formula=CEL~1,data=dt2,ID=ID,
		      Vcode=Vcode, 
		      p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
		       model="G", 
		      IPRT=TRUE) 



## compute the estimates of the conditional probabilities with sum of the new repeated measure as a summary statistics 
Psum=index.batch(olmeNB=re.gamma.ar1,data=dt2,ID=ID,Vcode=Vcode,
		 labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE)  
		 

## compute the estimates of the conditional probabilities with max of the new repeated measure as a summary statistics 
Pmax=index.batch(olmeNB=re.gamma.ar1,data=dt2,ID=ID,Vcode=Vcode, labelnp=new,qfun="max", IPRT=TRUE)

## Which patient's estimated probabilities based on the sum and max statistics disagrees the most?
( IDBigDif <- which(rank(abs(Pmax$condProbSummary[,1]-Psum$condProbSummary[,1]))==80) )
## Show the patient's CEL counts  
dt2$CEL[ID==IDBigDif]
## Show the estimated conditional probabilities based on the sum summary statistics
Psum$condProbSummary[IDBigDif,]
## Show the estimated conditional probabilities based on the max summary statistics
Pmax$condProbSummary[IDBigDif,]


## 2) Fit the negative binomial mixed-effect AR(1) model where random effects is from the log-normal distribution

re.logn.ar1=mle.ar1.fun(formula=CEL~1,data=dt2,ID=ID,
		      Vcode=Vcode, 
		      p.ini=c(loga,logtheta,logitd,b0), ## log(a), log(theta), logit(d), b0
		       model="N", 
		      IPRT=TRUE)

## Require some time
Psum=index.batch(olmeNB=re.logn.ar1,data=dt2,ID=ID,Vcode=Vcode, labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE) 
re.logn.ar1$Psum=Psum




## 3) Fit the negative binomial independent model where random effects are from the gamma distribution
re.gamma.ind=mle.fun(formula=CEL~1,data=dt2,ID=ID, 
                   model="G", 
		   p.ini=c(loga,logtheta,b0), 
		   IPRT=TRUE)

Psum=index.batch(olmeNB=re.gamma.ind,data=dt2,ID=ID, 
                 labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE)  




## 4) Fit the negative binomial independent model where random effects are from the lognormal distribution
re.logn.ind=mle.fun(formula=CEL~1,data=dt2,ID=ID, 
                   model="N", 			   	
		   p.ini=c(loga,logtheta,b0), 		
		   IPRT=TRUE)

Psum=index.batch(olmeNB=re.logn.ind, data=dt2,ID=ID, 
                 labelnp=new,qfun="sum", IPRT=TRUE,i.se=TRUE) 


## 5) Fit the semi-parametric negative binomial AR(1) model 

logvarG = -0.5

re.semi.ar1=mle.ar1.non3(formula=CEL~1,data=dt2,ID=ID, p.ini=c(loga, logvarG, logitd,b0),Vcode=Vcode)

Psum=index.batch(olmeNB=re.semi.ar1,data=dt2,ID=ID, Vcode=Vcode,
		 labelnp=new,qfun="sum", IPRT=TRUE,i.se=FALSE)  ## No option of i.se=TRUE



## 6) Fit the semi-parametric negative binomial independent model 
re.semi.ind=mle.a3.fun(formula=CEL~1,data=dt2,ID=ID, p.ini=c(loga, logvarG, b0))
Psum=index.batch(olmeNB=re.semi.ind,data=dt2,ID=ID,  labelnp=new, qfun="sum", IPRT=TRUE,i.se=FALSE) ## No option of i.se=TRUE



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
main="Simulated from the AR(1) model \n with random effect ~ gamma")

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



%This function fits a negative binomial mixed-effect AR(1) model in the formulation described Zhao et al.
%
%Denote a repeated measure from a subject \emph{i} at time point \emph{j} as \emph{Y_ij} where \emph{i=1,...,N,j=1,...,n_i}.
%
%The model assumes that the repeated measure at time point \emph{1} given the subject specific random effect \emph{G_i=g_i} to follow the negative binomial distribution:
%

%\emph{Y_i1| G_i=g_i ~ NB(r_i1,p_i),} 

%where the failure probability \emph{ p_i } at each time point \emph{j} of subject \emph{i} is parametrized as: 


%\emph{ p_i= 1/(g_i*a+1),}  


%and \emph{ a > 0 }.
%
%For \emph{j > 1}, the model assumes the first-order Markov chain:
%

%\emph{ Y_ij = Z(Y_i(j-1),d) + e_ij, }

%
%with \emph{0 < delta < 1} and \emph{r_ij-d*r_i(j-1) > 0}.
% 
%Given Y_i(j-1) and d, Z is from the beta binomial distribution:


%\emph{ Z(Y_i(j-1),d) ~ betabinom(n=Y_i(j-1),alpha=d*r_i,(j-1),beta=(1-delta)*r_i(j-1))} 


%and \emph{e_ij} is distributed as the negative binomial. 


%\emph{e_ij ~ NB(r_ij-d*r_i(j-1),p)i)}


%These assumptions lead that the response counts to be from \emph{ Y_ij | G_i =g_i ~ NB(r_ij,p_i)} as the independence model (\code{\link{mle.fun}}) and \emph{Y_ij} is conditionally independent on \emph{Y_ij'} given \emph{Y_i(j-1)} for all \emph{j' < j-1 }.


%Finally, the parameter \emph{d} can be interpreted as a measure of associations between repeated measures of a subject: the larger \emph{d}, the stronger the positive correlations are.
%This interpretation depends on the result:
% \emph{cov(Y_ij,Y_ij'|G_i=g_i)=d^{j-j'} E(Y_ij')*(a*g_i^2+g_i)},



%The distributional assumptions on \emph{G_i} are the same as the independent model. See \code{\link{mle.fun}} for details.
