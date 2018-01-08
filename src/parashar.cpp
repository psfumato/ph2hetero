// STL
#include <iostream>
#include <algorithm> 
#include <cmath>

// SCYTHE
#include "mersenne.h"
#include "rng.h"
#include "distributions.h"
#include "ide.h" 
#include "la.h"
#include "matrix.h" 
#include "stat.h" 
#include "smath.h" 
// R interface
#include <R.h>           
#include <R_ext/Utils.h> 
#include <Rdefines.h>
#include <Rinternals.h>
#include "usefunc.h"
#include "extra.h"

using namespace scythe;
using namespace std;


extern "C" {
  void Dparashar(const double *alpha_in, 
     const double *beta_in, 
     const double *RR0_, 
     const double *RRn_, 
     const double *RRp_, 
     const int *Nmax_,
     int *k1n_,
		 int *k1p_,
     int *N1n_,
     int *N1p_,
     int *kep_, 
		 int *Nep_, 
     int *kn_,
		 int *kp_,
     int *Nn_,
     int *Np_,
     double *alpha_out , 
     double *power_out, 
  	 double *PET, 
     double *EN, 
     const int *NumThreads_
		 ) {
  
  const double RR0=*RR0_;
  const double RRn=*RRn_;
  const double RRp=*RRp_;
  const double alpha=*alpha_in;
  const double beta=*beta_in;
  const int Nmax=*Nmax_;
  int Nmax1=Nmax;
  int X1n;
  int X2n;
  int k1n;
  int kn;
  int X1p;
  int X2p;
  int kp;
  int X2ep;
  int k1p;
  int kep;

#ifdef _OPENMP
omp_set_num_threads(*NumThreads_);
#pragma omp parallel shared(Nmax1) private(X1n,X2n,k1n,kn,X1p,X2p,k1p,kp,X2ep,kep)
{
#pragma omp for schedule(dynamic)
#endif
  
  for(X1n=0; X1n<=Nmax1; X1n++){
    int k1nmin=1;
    if (round(X1n*RR0)-1>1){k1nmin=round(X1n*RR0)-1;}
    int k1nmax=round(X1n*RRn)+1;
    for (X2n=0; X2n<=Nmax1; X2n++){
      int Xn=X1n+X2n;
      for(k1n=k1nmin; k1n<=k1nmax; k1n++){
        int knmin=k1n;
        if (round(Xn*RR0)-1>knmin){knmin=round(Xn*RR0)-1;}
        int knmax=round(Xn*RRn)+1;
        for(kn=knmin; kn<=knmax; kn++){  
          Matrix<double>R1=R1_calc_parashar(X1n,X2n,k1n,kn,RR0,RRn);
          if((R1[0]>alpha)||(R1[1]<1-beta)){continue;}
          double EN1_sum=EN1_sum_calc_parashar(X1n,k1n,kn,RR0);
          for(X1p=0; X1p<=Xn; X1p++){
            if (X1n+X1p>Nmax1){continue;}
            if ((X1n+X1p>EN[0])&(EN[0]>0)){continue;}
             int k1pmin=1;
             if (round(X1p*RR0)-1>1){k1pmin=round(X1p*RR0)-1;}
             int k1pmax=round(X1p*RRp)+1;
            for(X2p=0; X2p<=Xn; X2p++){
              int Xp=X1p+X2p;
              if (Xp>Xn){continue;}
              if ((X1n+X1p+(Xn+Xp-(X1n+X1p))*EN1_sum>EN[0])&(EN[0]>0)){continue;}
              for (k1p=k1pmin; k1p<=k1pmax; k1p++){
                int kpmin=k1p;
                if (round(Xp*RR0)-1>kpmin){kpmin=round(Xp*RR0)-1;}
                int kpmax=round(Xp*RRp)+1;
                for(kp=kpmin; kp<=kpmax;kp++){
                Matrix<double>R2=R2_calc_parashar(X1n,X2n,Xp,k1n,kp,kn,RR0,RRp);
                if (R1[0]+R2[0]>alpha){continue;}
                  for (X2ep=0; X2ep<=Xn; X2ep++){
                  int Nmax2=max(X1n+X1p+X2ep,Xn+Xp);
                  if (Nmax2>Nmax1){continue;}
                  int Xep=X1p+X2ep;
                  int kepmin=k1p;
                  if (round(Xep*RR0)-1>kepmin){kepmin=round(Xep*RR0)-1;}
                  int kepmax=round(Xep*RRp)+1;
                    for (kep=kepmin; kep<=kepmax; kep++){
                    Matrix<double> R3=R3_calc_parashar(X1n,X1p,X2ep,k1n,k1p,kep,RR0,RRp);
                    if ((R1[0]+R2[0]+R3[0]>alpha)||(min(R1[1],R2[1]+R3[1])<1-beta)){continue;}
                    double EN2_sum=EN2_sum_calc_parashar(X1p,k1p,kep,RR0);
                    const double k1n1d=k1n-1;
                    if ((X1n+X1p+(Xn+Xp-(X1n+X1p))*EN1_sum+(Xep-X1p)*pbinom(k1n1d,X1n,RR0)*EN2_sum>EN[0])&(EN[0]>0)){continue;}
                    EN[0]=X1n+X1p+(Xn+Xp-(X1n+X1p))*EN1_sum+(Xep-X1p)*pbinom(k1n1d,X1n,RR0)*EN2_sum;
                    k1n_[0]=k1n;
                    k1p_[0]=k1p;
                    N1n_[0]=X1n;
                    N1p_[0]=X1p;
                    kep_[0]=kep;
                    Nep_[0]=Xep;
                    kn_[0]=kn;
                    kp_[0]=kp;
                    Nn_[0]=Xn;
                    Np_[0]=Xp;
                    alpha_out[0]=R1[0]+R2[0]+R3[0];
                    power_out[0]=min(R1[1],R2[1]+R3[1]);    
                    const double kn1d=kn-1;
                    const double kep1d_=kep-1;
                    const double k1p1d=k1p-1;
                    PET[0]=1-pbinom(kn1d,X1n,RR0) + pbinom(k1n1d,X1n,RR0)*(1-pbinom(kep1d_,X1p,RR0)+pbinom(k1p1d,X1p,RR0));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
#ifdef _OPENMP
}
#endif
  
	}         
} // extern "C"
