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
  void Djones(const double *alpha_in, 
     const double *beta_in, 
     const double *RR0_, 
     const double *RRn_, 
     const double *RRp_, 
     const double *Ppos_,
     const int *Nmax_,
     const int *k1n_,
  	 const int *k1p_,
     const int *N1n_,
     const int *N1p_,
     int *N2u_, 
		 int *kn_, 
     int *kp_,
		 int *N2p_, 
     int *k2p_,
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
  const double Ppos=*Ppos_;
  const int N1n=*N1n_;
  const int N1p=*N1p_;
  const int k1n=*k1n_;
  const int k1p=*k1p_;
  const int Nmax=*Nmax_;
  int Nmax1=Nmax;
  int X2u;
  int kn;
  int kp;
  int X2p;
  int k2p;
            
            
PET[0]=PET_calc_jones(N1n,N1p,k1n,k1p,RR0);          
double EN1_sum = EN1_sum_calc_jones(N1n,k1n,N1p,RR0);           
double EN2_sum = EN2_sum_calc_jones(N1p,k1p,N1n,k1n,RR0);

#ifdef _OPENMP
omp_set_num_threads(*NumThreads_);
#pragma omp parallel shared(Nmax1) private(X2u,X2p,kn,kp,k2p) 
{
#pragma omp for schedule(dynamic)
#endif
  for(X2u=0; X2u<=Nmax1; X2u++){
    if ((PET[0]*(N1n+N1p) + EN1_sum*(N1n+N1p+X2u)>EN[0])&(EN[0]>0)){continue;}
    int e2n=round(X2u*(1-Ppos));
    int e2p=X2u-e2n;
    int knmin=k1n+1;
    if (round((N1n+e2n)*RR0)>knmin){knmin=round((N1n+e2n)*RR0);}
    int knmax=round((N1n+e2n)*RRn)+1;
    for(X2p=0; X2p<=Nmax1; X2p++){
      int Nmax2=max(N1n+N1p+X2u,N1n+N1p+X2p);
      if (Nmax2>Nmax1){continue;}
      for (kn=knmin; kn<=knmax; kn++){
        Matrix<double> R1=R1_calc_jones(N1n,e2n,k1n,kn,RR0,RRn);
        if((R1[0]>alpha)||(R1[1]<1-beta)){continue;}
        int kpmin=k1p+1;
        if (round((N1p+e2p)*RR0)>kpmin){kpmin=round((N1p+e2p)*RR0);}
        int kpmax=round((N1p+e2p)*RRp)+1;
        for (kp=kpmin; kp<=kpmax; kp++){
          Matrix<double> R2=R2_calc_jones(N1n,N1p,e2n,e2p,k1n,kp,kn,RR0,RRp);
          if (R1[0]+R2[0]>alpha){continue;}
          if ((PET[0]*(N1n+N1p) + EN1_sum*(N1n+N1p+X2u) + EN2_sum*(N1n+N1p+X2p)>EN[0])&(EN[0]>0)){continue;}
          int k2pmin=k1p+1;
          if (round((N1p+X2p)*RR0)>k2pmin){k2pmin=round((N1p+X2p)*RR0);}
          int k2pmax=round((N1p+X2p)*RRp)+1;
          for (k2p=k2pmin; k2p<=k2pmax; k2p++){
            Matrix<double> R3=R3_calc_jones(N1n,N1p,X2p,k1n,k1p,k2p,RR0,RRp);
            if ((R1[0]+R2[0]+R3[0]>alpha)||(min(R1[1],R2[1]+R3[1])<1-beta)){continue;}
            EN[0]=PET[0]*(N1n+N1p) + EN1_sum*(N1n+N1p+X2u) + EN2_sum*(N1n+N1p+X2p);
            alpha_out[0]=R1[0]+R2[0]+R3[0];
            power_out[0]=min(R1[1],R2[1]+R3[1]);    
            N2u_[0]=X2u;
            kn_[0]=kn;
            kp_[0]=kp;
            N2p_[0]=X2p;
            k2p_[0]=k2p;
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
