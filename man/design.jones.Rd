\name{design.jones}
\alias{design.jones}

\title{Jones 2-stage Phase II design}
       
       \description{
Calculates Optimal 2-stage Phase II designs given by Cheryl L. Jones (2007).}
       
\usage{
  design.jones(alpha=0.05,beta=0.2,p0,p1n,p1p,Nmax=100,Ppos=0.5,NumThreads=2)
}
\arguments{
       \item{alpha}{Type I error.}
       
       \item{beta}{Type II error.}
       
       \item{p0}{Max unacceptable response rate.}
       
       \item{p1n}{Min acceptable response rate for "negative" subgroup.}
       
       \item{p1p}{Min acceptable response rate for "positive" subgroup.}
       
       \item{Nmax}{Maximum sample size.}
       
       \item{Ppos}{Prevalence of "biomarker-positive" subjects.}
       
       \item{NumThreads}{Number of threads used for parallel compilation.}
       }

\value{
  Returns a data.frame object which components are :
  
  \item{alpha}{Type I error probability of the optimal design.}

  \item{power}{Power of the optimal design.}

  \item{p0}{Max unacceptable response rate.}
       
  \item{p1n}{Min acceptable response rate for "negative" subgroup.}
       
  \item{p1p}{Min acceptable response rate for "positive" subgroup.}

  \item{PET}{Overall probability of early termination.}
  
  \item{EN}{Expected sample size.}
  
  \item{k1n}{Minimum number of responses in marker negative subjects to go into the unselected population during stage 2.}

  \item{k1p}{Minimum number of responses in the marker positive subjects.}
  
  \item{N1n}{Number of marker negative subjects enrolled during the first stage.}

  \item{N1p}{Number of marker positive subjects enrolled during the first stage.}

  \item{k2p}{Total number of response required during stage 1 and stage 2 combined in the amplified when preliminary efficacity is in the amplified subgroup.}
  
  \item{N2p}{Number of amplified subjects enrolled during the second stage if preliminary evidence of efficacy leads to this path.}
  
  \item{kn}{Total number of response required during stage 1 and stage 2 combined in the negative sample to conclude efficacy in the unselected population.}
  
  \item{kp}{Total number of response required during stage 1 and stage 2 combined in the amplified sample to conclude efficacy in this population when preliminary efficacy is unselected.}
    
  \item{N2un}{Number of unselected subjects enrolled during the second stage if preliminary evidence of efficacy leads to this path.}
  
}

\examples{
out.jones<-design.jones(alpha=0.05,
                        beta=0.2,
                        p0=0.03,
                        p1n=0.10,
                        p1p=0.15,
                        Nmax=150,
                        Ppos=0.40)
}
\author{
  Patrick Sfumato and Bastien Cabarrou.
}

\references{

 Jones CL, Holmgren E (2007).\emph{An adaptive Simon Two-Stage Design for Phase 2 studies of targeted therapies.} Contemporary Clinical Trials 28 654-661.}