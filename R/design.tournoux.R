design.tournoux <- function(alpha=0.05,beta=0.2,p0n,p0p,p1n,p1p,w=1,gamma=0.6){

if (alpha<=0||alpha>=1){stop("alpha value must be included in a ]0,1[ range")}
if (beta<=0||beta>=1){stop("beta value must be included in a ]0,1[ range")}
if ((p0n<=0||p0n>=1)||(p0p<=0||p0p>=1)||(p1n<=0||p1n>=1)||(p1p<=0||p1p>=1)){stop("Probabilities values must be included in a ]0,1[ range")}
if (p0n>=p1n){stop("p0n value must be strictly less than p1n value")}
if (p0p>=p1p){stop("p0p value must be strictly less than p1p value")}
if (gamma<=0||gamma>=1){stop("gamma value must be included in a ]0,1[ range")}
if (w<=0){stop("w value must be strictly positive")}

#Two stage Fleming function
  
fleming.two.stage.design<- function(p0,p1,n1,n2,alpha,beta){
pa=(sqrt((n1+n2)*p0)+sqrt(1-p0)*qnorm(1-alpha))^2/((n1+n2+qnorm(1-alpha)^2))
r1=round(n1*p0+qnorm(1-alpha)*sqrt((n1+n2)*p0*(1-p0)))+1
a1=round(n1*pa-qnorm(1-alpha)*sqrt((n1+n2)*pa*(1-pa)))
    
r2=round((n1+n2)*p0+qnorm(1-alpha)*sqrt((n1+n2)*p0*(1-p0)))+1
a2=round((n1+n2)*pa-qnorm(1-alpha)*sqrt((n1+n2)*pa*(1-pa)))
prod1=0
prod2=0
for (m in (a1+1):(r1-1)){
    prod1=prod1+dbinom(m,n1,p0)*(1-pbinom(r2-m-1,n2,p0))
    prod2=prod2+dbinom(m,n1,p1)*pbinom(a2-m,n2,p1)}
    alpha.out=prod1+1-pbinom(r1-1,n1,p0)
    beta.out=prod2+pbinom(a1,n1,p1)
    
return(list(r1=r1,a1=a1,r2=r2,a2=a2,alpha.out=alpha.out,beta.out=beta.out))}
  
#Function to calculate heterogeneity parameter
  
c.psi.func <- function(N1,N2,p0n,p0p,n=10000,gamma){
    
R1 <- rbinom(n,N1,p0n)
R2 <- rbinom(n,N2,p0p)
    
d1 <- (R1/N1) - p0n
d2 <- (R2/N2) - p0p
d <- abs(d1) + abs(d2)
S <- sign(d1) + sign(d2)
    
c <- 0
g <- 1
while(g>(gamma*gamma)/2){
  c <- c +0.01
  g <- sum(d>c & S==0)/n  
}
    
return(c)
}
  
  #One stage Fleming
  
  p0 <- (p0n+w*p0p)/(1+w)
  p1 <- (p1n+w*p1p )/(1+w)
  
  N.norm.prep <- ((qnorm(1-beta)*sqrt(p1*(1-p1))+ qnorm(1-alpha)*sqrt(p0*(1-p0))) / ((p1-p0)))**2
  N.norm=trunc(N.norm.prep)+1
  N.range=trunc(1/4*N.norm):trunc(7/4*N.norm+1)
  N=N.range[which(N.range%%2==0)]
  R<-round((N*p0+qnorm(1-alpha)*sqrt(N*p0*(1-p0))))+1
  alpha.one.stage.fleming.design=1-pbinom(R-1,N,p0)
  beta.one.stage.fleming.design=pbinom(R-1,N,p1)
  
  n1=NULL
  for (i in 1:length(N)){
    if (N[i]%%2==0){n1[i]=N[i]/2}
    if (N[i]%%2!=0){n1[i]=(N[i]+1)/2}}
  n2=N-n1
  
  one.stage.fleming.design<-data.frame(N,alpha.one.stage.fleming.design,beta.one.stage.fleming.design,n1=n1,n2=n2)
  
  #Two stage Fleming
  
  i=1
  alpha.two.stage.fleming.design=NULL
  beta.two.stage.fleming.design=NULL
  r1=NULL
  r2=NULL
  a1=NULL
  a2=NULL
  for (i in 1:length(N)){
    n1=one.stage.fleming.design$n1[i]
    n2=one.stage.fleming.design$n2[i]
    out<-fleming.two.stage.design(p0,p1,n1,n2,alpha,beta)
    alpha.two.stage.fleming.design[i]=out$alpha.out
    beta.two.stage.fleming.design[i]=out$beta.out
    r1[i]=out$r1
    r2[i]=out$r2
    a1[i]=out$a1
    a2[i]=out$a2
    i=i+1
  }
  
  two.stage.fleming.design=data.frame(one.stage.fleming.design,a1,a2,r1,r2,alpha.two.stage.fleming.design,beta.two.stage.fleming.design)
  
  fleming.select=two.stage.fleming.design[which((two.stage.fleming.design$alpha.one.stage.fleming.design<=alpha)
                                                &(two.stage.fleming.design$beta.one.stage.fleming.design<=beta)
                                                &(two.stage.fleming.design$alpha.two.stage.fleming.design<=alpha)
                                                &(two.stage.fleming.design$beta.two.stage.fleming.design<=beta)
                                                &((two.stage.fleming.design$n1*w)%%(w+1)==0)),]
  
  if (dim(fleming.select)[1]==0){stop("Design not found for these parameters")}
  
  fleming.select.min=fleming.select[1,]
  
  n1=fleming.select.min$n1
  
  N=fleming.select.min$N
  n1p=n1*w/(w+1)
  n1n=n1-n1p
  a1=fleming.select.min$a1
  b1=fleming.select.min$r1
  
  #If psi=0
  
  n2n=n1n
  n2p=N-(n1n+n1p+n2n)
  a2=fleming.select.min$a2
  b2=fleming.select.min$r2
  
  #If psi=1
  
  n2Fn=0
  alpha.F1=1
  beta.F1=1
  while((alpha.F1>alpha)|(beta.F1>beta)){
    n2Fn=n2Fn+1
    out<-fleming.two.stage.design(p0n,p1n,n1n,n2Fn,alpha,beta)
    a2Fn<-out$a2
    b2Fn<-out$r2
    alpha.F1<-out$alpha.out
    beta.F1<-out$beta.out
  }
  
  #If psi=2
  
  n2Fp=0
  alpha.F2=1
  beta.F2=1
  while((alpha.F2>alpha)|(beta.F2>beta)){
    n2Fp=n2Fp+1
    out<-fleming.two.stage.design(p0p,p1p,n1p,n2Fp,alpha,beta)
    a2Fp<-out$a2
    b2Fp<-out$r2
    alpha.F2<-out$alpha.out
    beta.F2<-out$beta.out
  }
  
  #Calculate c1, c2, c2F1, c2F2
  
  c1<-c.psi.func(N1=n1n,N2=n1p,p0n,p0p,n=10000,gamma)
  c2<-c.psi.func(N1=n1n+n2n,N2=n1p+n2p,p0n,p0p,n=10000,gamma)
  
  Tournoux.design=data.frame(Fleming=c("One stage","Two stage","","",""),
                             Stage=c("","Stage 1","Stage 2","",""),
                             psi=c("","",0,1,2),
                             Neg=c("",n1n,n2n,n2Fn,0),
                             Pos=c("",n1p,n2p,0,n2Fp),
                             a=c("",a1,a2,a2Fn,a2Fp),
                             b=c("",b1,b2,b2Fn,b2Fp),
                             c=c("",c1,c2,"",""),
                             N=c(N,n1n+n1p,n1n+n1p+n2n+n2p,n1n+n1p+n2Fn,n1n+n1p+n2Fp),
                             alpha=c(round(fleming.select.min$alpha.one.stage.fleming.design,4),"",
                                     round(fleming.select.min$alpha.two.stage.fleming.design,4),
                                     round(alpha.F1,4),round(alpha.F2,4)),
                             power=c(round(1-fleming.select.min$beta.one.stage.fleming.design,4),"",
                                     round(1-fleming.select.min$beta.two.stage.fleming.design,4),
                                     round(1-beta.F1,4),round(1-beta.F2,4)),row.names=NULL)
  
  return(Tournoux.design)}
