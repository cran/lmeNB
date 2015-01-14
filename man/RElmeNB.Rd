
\name{RElmeNB}
\alias{RElmeNB}
\alias{int.numRE}
\alias{int.denRE}
\alias{dens_Yi.gY}
\alias{densYijGivenYij_1AndGY}
\alias{dbb}

%- Also NEED an '\alias' for EACH other topic documented here.
%
\title{
Calculate predicted values of E(Gi|Yi) given the estimates of parameters
}

\description{
 
  Compute predicted values of random effects for each patient
  
}


\usage{
RElmeNB(theta, alpha, betas, delta, formula, ID, Vcode = NULL,
        data, AR, RE, rel.tol = .Machine$double.eps^0.8, expG = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
 A scalar containing the estimated variance of the random effect distribution, \eqn{\theta}{theta}.
  }

  \item{alpha}{
  A scalar containing the estimated dispersion parameter, \eqn{\alpha}{alpha}.
  }
  %
  
  \item{betas}{
  A vector containing the estimated regression coefficients, \eqn{\boldsymbol{\beta}}{beta}.
  }
  %
  \item{delta}{
   AR(1) parameter, \eqn{\delta}{delta}
  }

  \item{ID}{
   See \code{\link{lmeNB}}.
  }

  \item{Vcode}{
   Necessary only if the AR(1) model is used.
   See \code{\link{lmeNB}}.
  }
  
  \item{RE}{
  The distribution of random effects \eqn{G_i}{G[i]}.
  If \code{RE="G"} then the random effects are assumed to be from the gamma distribution.
  If \code{RE="N"} then they are assumed to be form the log-normal distribution.

  The current version of RElmeNB only accept parametric model.
  }
  \item{AR}{
  See \code{\link{lmeNB}}.
  }  

\item{formula}{ See \code{\link{lmeNB}}.}
\item{data}{ See \code{\link{lmeNB}}.}
\item{rel.tol}{relative tolerance for the integration of the random effect. passed to \code{integrate} function.}
\item{expG}{ Internal use only}
}


\value{
return the predicted RE of each patient. 
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

The subroutines of this function is:
\code{\link{fitParaIND}},
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
## See the examples in help files of rNBME.R.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
