design.parashar <- function(alpha=0.05,beta=0.2,p0,p1n,p1p,Nmax=100,NumThreads=2) {
  
  if (alpha<=0||alpha>=1){stop("alpha value must be included in a ]0,1[ range")}
  if (beta<=0||beta>=1){stop("beta value must be included in a ]0,1[ range")}
  if ((p0<=0||p0>=1)||(p1n<=0||p1n>=1)||(p1p<=0||p1p>=1)){stop("Probabilities values must be included in a ]0,1[ range")}
  if ((p0>=p1n)||(p0>=p1p)){stop("p0 value must be strictly less than p1n and p1p values")}
  if (p1n>p1p){stop("p1n value must be less than p1p value")}
  if (NumThreads<=0){stop("NumThreads value must be strictly greater than 0")}
  if (Nmax<=0){stop("Nmax value must be strictly greater than 0")}
          
          
    output <- .C("Dparashar",
               alpha_in=as.double(alpha),
               beta_in=as.double(beta),
               RR0_=as.double(p0),
               RRn_=as.double(p1n),
               RRp_=as.double(p1p),
               Nmax_=as.integer(Nmax),
               k1n_=as.integer(0),
               k1p_=as.integer(0),
               N1n_=as.integer(0),
               N1p_=as.integer(0),
               kep_=as.integer(0),
               Nep_=as.integer(0),
               kn_=as.integer(0),
               kp_=as.integer(0),
               Nn_=as.integer(0),
               Np_=as.integer(0),
               alpha_out=as.double(0),
               power_out=as.double(0),
               PET=as.double(0),
               EN=as.double(0),  
               NumThreads_=as.integer(NumThreads),
               pkg="ph2hetero")
  
  if(output$EN==0){warning("No optimal design using Parashar's method were found for theses parameters. Increase maximum sample size (Nmax).")}
  
  res<-data.frame(design="Parashar optimal design",
             alpha=output$alpha_out,
             power=output$power_out,
             p0=output$RR0_,
             p1n=output$RRn_,
             p1p=output$RRp_,
             PET=output$PET,
             EN=output$EN,
             k1n=output$k1n_,
             k1p=output$k1p_,
             N1n=output$N1n_,
             N1p=output$N1p_,
             kep=output$kep_,
             Nep=output$Nep_,
             kn=output$kn_,
             kp=output$kp_,
             Nn=output$Nn_,
             Np=output$Np_)
  
  colnames(res)[1]=""
  return(res)
}