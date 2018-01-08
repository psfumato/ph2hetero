\name{design.tournoux}
\alias{design.tournoux}

\title{Tournoux 2-stage Phase II design}
       
       \description{
Calculates Optimal 2-stage Phase II designs given by Caroline Tournoux-Facon (2011).}
       
\usage{
  design.tournoux(alpha=0.05,beta=0.2,p0n,p0p,p1n,p1p,w=1,gamma=0.6)
}
\arguments{
       \item{alpha}{Type I error.}
       
       \item{beta}{Type II error.}
       
       \item{p0n}{Max unacceptable response rate for "negative" subgroup.}
       
       \item{p0p}{Max unacceptable response rate for "positive" subgroup.}
       
       \item{p1n}{Min acceptable response rate for "negative" subgroup.}

       \item{p1p}{Min acceptable response rate for "positive" subgroup.}
       
       \item{w}{Ratio between the two subgroups \eqn{\frac{N_{pos}}{N_{neg}}}.}

       \item{gamma}{Level of the heterogeneity test.}
       }

\value{
  Returns a data.frame object which components are number of patients, error-rates and power-levels for each stage. The first line of the data.frame corresponds to the arguments for an One-stage Fleming design and the others are two-stages Fleming designs with : 
  
  - Non-stratified heterogeneous Fleming two-stage design when psi=0.
  
  - Stratified adaptive Fleming two-stage with identification of heterogeneity of responses in favour of "negative" subgroup when psi=1. 
  
  - Stratified adaptive Fleming two-stage with identification of heterogeneity of responses in favour of "positive" subgroup when psi=2. 
}

\examples{
out.tournoux<-design.tournoux(p0n=0.15,
                              p0p=0.15,
                              p1n=0.30,
                              p1p=0.25,
                              w=2,
                              gamma=0.6,
                              alpha=0.05,
                              beta=0.1)
}
\author{
  Patrick Sfumato and Bastien Cabarrou.
}

\references{
 Tournoux-Facon C, De Rycke Y and Tubert-Bitter P (2011).\emph{Targeting population entering phase III trials: A new stratified adaptive phase II design.} Statistics in Medecine 30(8) 801-11.
 
 Tournoux-Facon C, De Rycke Y and Tubert-Bitter P (2011).\emph{How a new stratified adaptive phase II design could improve targeting population.} Statistics in Medecine 30(13) 1555-62.
 }