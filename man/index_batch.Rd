\name{index.batch}
\alias{index.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the point estimate and its 95% confidence interval (for a parametric model) 
of the conditional probability \emph{ Pr(q(Y_i,new)>=q(y_i,new)| Y_i,pre=y_i,pre)}
}

\description{
%
Given an output of \code{\link{mle.fun}}, \code{\link{mle.ar1.fun}},
\code{\link{mle.a3.fun}}, \code{\link{mle.ar1.non3}} or \code{\link{lmeNB}}.
This function computes the probability of observing the response counts as large as those new observations of subject \emph{i} \emph{y_i,new=(y_i,m_i+1,...,y_i,ni)} conditional on the subject's previous observations \emph{y_i,pre(y_i,1,...,y_i,m_i)}.
%
%
That is, this function returns a point estimate and its asymptotic 95% confidence interval (for a parametric model) of the conditional probability for each subject:
%

\emph{ Pr(q(Y_i,new)>=q(y_i,new)| Y_i,pre=y_i,pre)}.

The standard error is not produced when the semi-parametric approach is employed.
%
A scalar statistic to summarize the new response counts can be either the total count, \emph{q(Y_i,new)=sum_{j=m_i+1}^ni Y_ij}, or the maximum, \emph{q(Y_i,new)=max \{ Y_ij;j=m_i+1,...,n_i \} }.
%
See Zhao et al., for more details.

}



\usage{
index.batch(data,labelnp, ID, Vcode,olmeNB,subset=NULL, 
            qfun = "sum", IPRT = TRUE,  i.se = TRUE, iMC = FALSE)	
}


\arguments{
  \item{data}{
    A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
    The each row must contains the data corresponding to the repeated measure j of subjects and the rows \emph{(i,j)}s must be ordered as 
    \emph{(1,1),...,(1,n_1),(2,1),...,(2,n_2),...,(N,n_N)}.
    This dataset does not have to be the same as the one used in the
  computations of negative binomial mixed effect regression (\code{\link{mle.fun}},\code{\link{mle.ar1.fun}},\code{\link{mle.a3.fun}},\code{\link{mle.ar1.non3}} or \code{\link{lmeNB}}).
  }
  \item{labelnp}{
    A vector of length the total number of repeated measures, indicating new measures by \code{TRUE}.
    For examples, if patient \emph{i} has a \emph{n_i} repeated measures and the last \emph{n_i-m_i+1} measures are new, then 
    then \code{labelnp=c(rep(rep(FALSE,m_1),rep(TRUE,n_1-m_1+1)),rep(rep(FALSE,m_2),rep(TRUE,n_2-m_2+1)),...,rep(rep(FALSE,m_N),rep(TRUE,n_N-m_N+1)))}.}
  \item{ID}{
   A vector of length \emph{ sum_i^N n_i }, containing the subject IDs. i.e., \code{c(rep(1,n_1),rep(2,n_2),...,rep(N,n_N))}
 }
 \item{Vcode}{
   A vector of length the total number of repeated measures, containing the indices of time point.
   For example, 
   there are three subjectss and two subjects do not have missing visits and completed five visits while the other subjects missed a visit at the third time and there are four visits in total, then
   \code{Vcode=c(1,2,3,4,5,1,2,3,4,5,1,2,4,5)}.
   Necessary only if the olmeNB is an output of AR(1) models from \code{\link{mle.ar1.fun}},\code{\link{mle.ar1.non3}}.
}
  \item{olmeNB}{
    Output of \code{\link{mle.fun}},\code{\link{mle.ar1.fun}},\code{\link{mle.a3.fun}},\code{\link{mle.ar1.non3}}
    or \code{\link{lmeNB}}.
  }
 \item{subset}{
   An optional expression indicating the subset of the subjects of that the index should be computed. 
 }
 
 \item{qfun}{
   If \code{qfun="sum"},
   a scalar statistic to summarize the new response counts is the total count.
   
   If \code{qfun="max"},
   a scalar statistic to summarize the new response counts is the maximum.
  
 }
 \item{IPRT}{
   print control.
 }
 
 
 \item{ i.se}{
   If \code{i.se=TRUE} then the standard errors of the estimator of the conditional probability are returned for the output of 
   \code{\link{mle.fun}} or \code{\link{mle.ar1.fun}}.
   The semi-parametric approach, \code{\link{mle.a3.fun}} or \code{\link{mle.ar1.non3}}, do not return the standard errors.
 }
 \item{iMC}{
   Necessary for the AR(1) outputs,\code{\link{mle.ar1.fun}} and \code{\link{mle.ar1.non3}}.
   If \code{iMC=TRUE} then the function \code{MCCP.ar1} is called and the Monte carlo integration is performed. Fast but could be unreliable; not recommended for computing the confidence intervals
   If \code{iMC=FALSE} then the function \code{CP1.ar1} is called and the adaptive quadrature is performed. Slow but reliable. 
 }
}


\details{
%
The standard error of the point estimate on the logit scale is constructed using the delta method for the parametric model, where distributional assumption was made for random effects. %put this in Details
%
}
\value{
The \emph{N} by 4 (3, if hide the SE) numeric matrix, containing the point estimate of the conditional probability, and the lower and the upper bounds of the 95 % confidence intervals. 
% hide the SE?
}

\references{
Zhao, Y., Li, D.K.B., Petkau, J.A., Riddehough, A. & Traboulsee, A.  Detection of unusual increases in MRI lesion counts in multiple sclerosis patients.
}
\author{
 Zhao, Y. and Kondo, Y.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
A wrapper function to fit negative binomial mixed effect model:
  \code{\link{lmeNB}}
  
The functions to fit the relevant models:
\code{\link{mle.fun}},
\code{\link{mle.ar1.fun}},
\code{\link{mle.a3.fun}},
\code{\link{mle.ar1.non3}},

The subroutines of \code{\link{index.batch}}:
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
## See the examples in linked R-functions. 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
