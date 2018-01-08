\name{design.parashar}
\alias{design.parashar}

\title{Parashar's 2-stage Phase II design}
       
       \description{
Calculates Optimal 2-stage Phase II designs given by Deepak Parashar (2016).}
       
\usage{
  design.parashar(alpha=0.05,beta=0.2,p0,p1n,p1p,Nmax=100,NumThreads=2)
}

\arguments{
       \item{alpha}{Type I error.}
       
       \item{beta}{Type II error.}
       
       \item{p0}{Max unacceptable response rate.}
       
       \item{p1n}{Min acceptable response rate for "negative" subgroup.}
       
       \item{p1p}{Min acceptable response rate for "positive" subgroup.}
       
       \item{Nmax}{Maximum sample size.}

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

  \item{kep}{Minimum number of responses in the enrichment population.}
  
  \item{Nep}{Number of subjects in the enrichment population.}
  
  \item{kn}{Total number of response required during stage 1 and stage 2 combined in the negative sample to conclude efficacy in the unselected population.}
  
    \item{kp}{Total number of response required during stage 1 and stage 2 combined in the amplified sample to conclude efficacy in this population when preliminary efficacy is unselected.}
  
  \item{Nn}{Number of marker negative subjects enrolled during the first stage and the second stage.}
  
  \item{N2p}{Number of marker positive subjects enrolled during the first stage and the second stage.}
}

\examples{
out.parashar<-design.parashar(alpha=0.05,
                              beta=0.2,
                              p0=0.03,
                              p1n=0.25,
                              p1p=0.40,
                              Nmax=30,
                              NumThreads=2)

}
\author{
  Patrick Sfumato and Bastien Cabarrou.
}

\references{

 Parashar D, Bowden J, Starr C, Wernisch L and Mander A (2016).\emph{An optimal stratified Simon two-stage design.} Pharmaceutical statistics 15(4) 333-40.}