
\name{lmeNB}
\alias{lmeNB}

%- Also NEED an '\alias' for EACH other topic documented here.
%
\title{
%
Performs the maximum likelihood estimation for the negative binomial
mixed-effect model.
This function is a wrapper for \code{\link{mle.fun}},
\code{\link{mle.ar1.fun}}, \code{\link{mle.a3.fun}} and \code{\link{mle.ar1.non3}}.
}

\description{
  %
  The model assumes that given the random effect \emph{G_i=g_i}, the count responses \emph{ Y_ij}s of subject \emph{i}, (\emph{i = 1, ..., N}), at time points \emph{j (=1,...,n_i)} follow the negative binomial distribution:
  %
  
  \emph{ Y_ij | G_i =g_i ~ NB(r_ij,p_i),} 
  
  %
  where \emph{ r_{ij} = exp(X_ij^T beta)/a} and \emph{beta} is fixed effect coefficients.
  The failure probability \emph{ p_i } of subject \emph{i}, is parametrized as: 
  
  \emph{ p_i= 1/(g_i*a+1),}  
  
  and \emph{ a > 0 }.
  % 
  The model assumes \emph{ E(G_i) = 1 } so that \emph{ E(Y_ij|G_i=g_i)=r_ij*g_i*a } and \emph{ E(Y_ij)= r_ij*g_i = exp(X_ij^T beta)}.
  %
  This assumption allows the interpretation of the latent random variable \emph{ G_i } as the subject i's activity rate relative to the overall cohort.
  %
  
  
  Regarding the dependence structures of \emph{Y_ij} and \emph{Y_ij'} conditional on the random effect, we consider two scenarios.
  
  [1]: \emph{Y_ij} and \emph{Y_ij'} are independent conditionally on \emph{Gi}. This assumption leads us to 
  %
  \emph{cov(Y_ij,Y_ij'|G_i=g_i)=0}
  \emph{cov(Y_ij,Y_ij')=E(Y_ij)^2*Var(Y_ij)+E(Y_ij)*(1+(Var(Y_ij)+1)*a)}
  
  [2]: AR(1) structures for \emph{Y_ij} and \emph{Y_ij'} conditionally on \emph{Gi}.
  
  Given emph{ G_i=g_i}, \emph{Y_ij} depends on \emph{Y_i(j-1)} through the beta binomial thinning and is conditionally independent on \emph{Y_ij'} given \emph{Y_i(j-1)} for all \emph{j' < j-1 }.
  %
  The beta binomial thinning operator depends on a  parameter \emph{d} which indicates the strength of the  positive AR(1) association between repeated measures of a subject: the larger \emph{d}, the stronger the positive correlations between the repeated measures of the same subject are.
  %
  This interpretation depends on the result:
  % 
  \emph{cov(Y_ij,Y_ij'|G_i=g_i)=d^{j-j'} E(Y_ij')*(a*g_i^2+g_i)}.
  %
  
  
  Regarding the random effect \emph{G_i}s distribution, \code{lmNB} allows
  three scenarios:
  
  (1) The log-normal distribution with \emph{E(G_i)=1} and \emph{Var(G_i)=theta}
  
  (2) The gamma distribution  with \emph{E(G_i)=1} and \emph{Var(G_i)=theta}.
  
  (3) No distributional assumption and the random effect distribution is approximated by the estimated values of the quantity:
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
  
  
}


\usage{
lmeNB(formula,data,ID,p.ini=NULL,IPRT=FALSE,AR=FALSE,RE=c("G","N","semipara"),deps=1e-04,Vcode,i.tol=1e-75,o.tol=1.e-3,maxit=100)
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
    The initial values of the parameters:
    
    If \code{AR=0} and (\code{RE="G"} or \code{RE="N"}),\code{p.ini}=\emph{(log(a), log(var(G)), beta0, beta1, ...)}.
    If \code{AR=0} and \code{RE="semipara"},            \code{p.ini}=\emph{(log(a), log(var(G)),beta0, beta1, ...)}.
    If \code{AR=1} and (\code{RE="G"} or \code{RE="N"}),\code{p.ini}=\emph{(log(a), log(var(G)),logit(d),beta0,beta1,...)}.
    If \code{AR=1} and \code{RE="semipara"},            \code{p.ini}=
    \emph{(log(a), log(var(G)),logit(d), b0, b1, ...)}.

    \code{NULL} is accepted.
  }
  \item{IPRT}{
    A logical, passed to Iprt of function \code{optim}. 
    If \code{TRUE} then print iterations. 
    
  }
  \item{AR}{WRITE SOMETHING}
  
  \item{RE}{
  The distribution of random effects \emph{G_i}.
  If \code{model="G"} then the random effects are assumed to be from the gamma distribution.
  If \emph{model="N"} then they are assumed to be form the log-normal distribution.
}
\item{deps}{aasdf}

\item{Vcode}{
A vector of length the total number of repeated measures, containing the indices of time point.
For example, 
there are three subjects and the first two subjects do not have missing visits and completed five visits while the last  subject missed the third visit  and have four visits in total, then
\code{Vcode=c(1,2,3,4,5,1,2,3,4,5,1,2,4,5)}. 
Necessary only if the AR(1) model is fit.
}
\item{i.tol}{
  A real number to determine the tolerance for \code{\link{integrate}}.
  Necessary only for semiparametric random effect model.
}
\item{o.tol}{
  A real number to determine the tolerance for \code{\link{optim}}.
  Necessary only for semiparametric random effect model.
}

\item{maxit}{
  The maximum number of iterations.
  Necessary only for semiparametric random effect model.
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

\item{idat}{
A dataframe, containing ID, CEL, x.1, x.2, ...  
%Vcode contains the scan number where -1,0,...,ni correspond to screening, baseline, and follow-up scans 1,...,ni.
The column labeled as CEL contains the response counts.
}


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
The subroutines of this function is:
\code{\link{mle.fun}},
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
## See the examples in help files of the subroutines.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
