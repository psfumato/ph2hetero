design.jones <-function(alpha=0.05,beta=0.2,p0,p1n,p1p,Nmax=100,Ppos=0.5,NumThreads=2) {
  
  if (alpha<=0||alpha>=1){stop("alpha value must be included in a ]0,1[ range")}
  if (beta<=0||beta>=1){stop("beta value must be included in a ]0,1[ range")}
  if ((p0<=0||p0>=1)||(p1n<=0||p1n>=1)||(p1p<=0||p1p>=1)){stop("Probabilities values must be included in a ]0,1[ range")}
  if ((p0>=p1n)||(p0>=p1p)){stop("p0 value must be strictly less than p1n and p1p values")}
  if (p1n>p1p){stop("p1n value must be less than p1p value")}
  if (NumThreads<=0){stop("NumThreads value must be strictly greater than 0")}
  if (Nmax<=0){stop("Nmax value must be strictly greater than 0")}
  if (Ppos<=0||Ppos>=1){stop("Ppos value must be included in a ]0,1[ range")}
      
  
  out.An<-ph2simon(pu=p0, pa=p1n, ep1=alpha/2, ep2=beta, nmax=Nmax)
  out.Ap<-ph2simon(pu=p0, pa=p1p, ep1=alpha/2, ep2=beta, nmax=Nmax)
  nopt.An=((1:nrow(out.An$out))[out.An$out[,5]==min(out.An$out[,5])])[1]
  nopt.Ap=((1:nrow(out.Ap$out))[out.Ap$out[,5]==min(out.Ap$out[,5])])[1]
  
  N1n<- out.An$out[c(nopt.An,1),][1,2]
  N1p<- out.Ap$out[c(nopt.Ap,1),][1,2]
  k1n<- out.An$out[c(nopt.An,1),][1,1]+1
  k1p<- out.Ap$out[c(nopt.Ap,1),][1,1]+1

  output <- .C("Djones",
               alpha_in=as.double(alpha),
               beta_in=as.double(beta),
               RR0_=as.double(p0),
               RRn_=as.double(p1n),
               RRp_=as.double(p1p),
               Ppos_=as.double(Ppos),
               Nmax_=as.integer(Nmax),
               k1n_=as.integer(k1n),
               k1p_=as.integer(k1p),
               N1n_=as.integer(N1n),
               N1p_=as.integer(N1p),
               N2u_=as.integer(0),
               kn_=as.integer(0),
               kp_=as.integer(0),
               N2p_=as.integer(0),
               k2p_=as.integer(0),
               alpha_out=as.double(0),
               power_out=as.double(0),
               PET=as.double(0),
               EN=as.double(0),
               NumThreads_=as.integer(NumThreads),
               pkg="ph2hetero")
  
  if(output$EN==0){warning("No optimal design using Jones method were found for theses parameters. Increase maximum sample size (Nmax).")}
  
  res<-data.frame(design="Jones optimal design",
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
                  k2p=output$k2p_,
                  N2p=output$N2p_,
                  kn=output$kn_,
                  kp=output$kp_,
                  N2un=output$N2u_
                  )
  
  colnames(res)[1]=""
  return(res)
}