// STL
#include <iostream>
#include <algorithm> // random_shuffle, reverse, sort, ...
#include <cmath>
// SCYTHE
#include "matrix.h" 
#include "distributions.h"
#include "ide.h" 
#include "la.h"
#include "mersenne.h"
#include "rng.h"
#include "stat.h" 
#include "smath.h" 
// R interface
#include <R.h>           //  Rprintf()
#include <R_ext/Utils.h> //  user interrupts
#include <Rdefines.h>
#include <Rinternals.h>
// 
#include "extra.h"
#include "usefunc.h"
#ifdef _OPENMP
#include <omp.h>
#endif


extern const double SDtol = 1e-9;
//*************************************************


using namespace scythe;
using namespace std;

/*Jones Design*/

Matrix<double> R1_calc_jones(const unsigned &X1n,          
            const unsigned &e2n,
            const unsigned &k1n,
            const unsigned &kn,            
            const double &RR0,
            const double &RRn) {
      
  Matrix<double> R1_out(1,2);
  double R1_a_prep=0;
  double R1_b_prep=0;

  for (unsigned i=k1n;i<=X1n;i++){
    const double id=i;
    const double knid=kn-id-1;
    R1_a_prep=R1_a_prep+(1-pbinom(knid,e2n,RR0))*dbinom(id,X1n,RR0);
    R1_b_prep=R1_b_prep+(1-pbinom(knid,e2n,RRn))*dbinom(id,X1n,RRn);
  }
  
  R1_out[0]=R1_a_prep;
  R1_out[1]=R1_b_prep;

  
  return(R1_out);
}

Matrix<double> R2_calc_jones(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &e2n,
            const unsigned &e2p,
            const unsigned &k1n,
            const unsigned &kp,
            const unsigned &kn,
            const double &RR0,
            const double &RRp) {
               
  Matrix<double> R2_out(1,2);
  double R2_ab_prep=0;


  for (unsigned i=k1n;i<=X1n;i++){
    const double id=i;
    const double knid=kn-id-1;
    R2_ab_prep=R2_ab_prep+(pbinom(knid,e2n,RR0))*dbinom(id,X1n,RR0);
  }
  
  const double kp1d=kp-1;

  R2_out[0]=(1-pbinom(kp1d,X1p+e2p,RR0))*R2_ab_prep;
  R2_out[1]=(1-pbinom(kp1d,X1p+e2p,RRp))*R2_ab_prep;
  
  return(R2_out);
}

Matrix<double> R3_calc_jones(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &X2p,
            const unsigned &k1n,
            const unsigned &k1p,
            const unsigned &k2p,
            const double &RR0,
            const double &RRp) {
      
  Matrix<double> R3_out(1,2);
  double R3_a_prep=0;
  double R3_b_prep=0;

  for (unsigned i=k1p;i<=X1p;i++){
    const double id=i;
    const double kepid=k2p-id-1;
    R3_a_prep=R3_a_prep+(1-pbinom(kepid,X2p,RR0))*dbinom(id,X1p,RR0);
    R3_b_prep=R3_b_prep+(1-pbinom(kepid,X2p,RRp))*dbinom(id,X1p,RRp);
  }
  
  const double k1n1d=k1n-1;

  R3_out[0]=pbinom(k1n1d,X1n,RR0)*R3_a_prep;
  R3_out[1]=pbinom(k1n1d,X1n,RR0)*R3_b_prep;
  
  return(R3_out);
}

double PET_calc_jones(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &k1n,
            const unsigned &k1p,
            const double &RR0) {
  
  double PET1_sum_out=0;

  for (unsigned k=0;k<k1n;k++){
    const double kd=k;
    PET1_sum_out=PET1_sum_out+dbinom(kd,X1n,RR0);
  }
  
  double PET2_sum_out=0;

  for (unsigned k=0;k<k1p;k++){
    const double kd=k;
    PET2_sum_out=PET2_sum_out+dbinom(kd,X1p,RR0);
  }
  
  double PET_out=PET1_sum_out*PET2_sum_out;
  return(PET_out);
}

double EN1_sum_calc_jones(const unsigned &X1n,
            const unsigned &k1n,
            const unsigned &X1p,
            const double &RR0) {
  
  double EN11_sum_out=0;

  for (unsigned k=k1n;k<=X1n;k++){
    const double kd=k;
    EN11_sum_out=EN11_sum_out+dbinom(kd,X1n,RR0);
  }
  
  double EN12_sum_out=0;

  for (unsigned k=0;k<=X1p;k++){
    const double kd=k;
    EN12_sum_out=EN12_sum_out+dbinom(kd,X1p,RR0);
  }
  
  double EN1_sum_out=EN11_sum_out*EN12_sum_out;
  
  return(EN1_sum_out);
}

double EN2_sum_calc_jones(const unsigned &X1p,
            const unsigned &k1p,
            const unsigned &X1n,
            const unsigned &k1n,
            const double &RR0) {
  
  double EN21_sum_out=0;

  for (unsigned l=0;l<k1n;l++){
    const double ld=l;
    EN21_sum_out=EN21_sum_out+dbinom(ld,X1n,RR0);
  }
  
  double EN22_sum_out=0;

  for (unsigned l=k1p;l<=X1p;l++){
    const double ld=l;
    EN22_sum_out=EN22_sum_out+dbinom(ld,X1p,RR0);
  }
  
  double EN2_sum_out=EN21_sum_out*EN22_sum_out;
  
  return(EN2_sum_out);
}

/*Parashar design*/

Matrix<double> R1_calc_parashar(const unsigned &X1n,          
            const unsigned &X2n,
            const unsigned &k1n,
            const unsigned &kn,            
            const double &RR0,
            const double &RRn) {
      
  Matrix<double> R1_out(1,2);
  double R1_a_prep=0;
  double R1_b_prep=0;
  const unsigned isupp=min(X1n,kn-1);

  for (unsigned i=k1n;i<=isupp;i++){
    const double id=i;
    const double knid=kn-id-1;
    R1_a_prep=R1_a_prep+(1-pbinom(knid,X2n,RR0))*dbinom(id,X1n,RR0);
    R1_b_prep=R1_b_prep+(1-pbinom(knid,X2n,RRn))*dbinom(id,X1n,RRn);
  }
  
  const double kn1d=kn-1;

  R1_out[0]=1-pbinom(kn1d,X1n,RR0)+R1_a_prep;
  R1_out[1]=1-pbinom(kn1d,X1n,RRn)+R1_b_prep;

  
  return(R1_out);
}

Matrix<double> R2_calc_parashar(const unsigned &X1n,
            const unsigned &X2n,
            const unsigned &Xp,
            const unsigned &k1n,
            const unsigned &kp,
            const unsigned &kn,
            const double &RR0,
            const double &RRp) {
               
  Matrix<double> R2_out(1,2);
  double R2_ab_prep=0;
  const unsigned isupp=min(X1n,kn-1);

  for (unsigned i=k1n;i<=isupp;i++){
    const double id=i;
    const double knid=kn-id-1;
    R2_ab_prep=R2_ab_prep+(pbinom(knid,X2n,RR0))*dbinom(id,X1n,RR0);
  }
  
  const double kp1d=kp-1;

  R2_out[0]=(1-pbinom(kp1d,Xp,RR0))*R2_ab_prep;
  R2_out[1]=(1-pbinom(kp1d,Xp,RRp))*R2_ab_prep;
  
  return(R2_out);
}

Matrix<double> R3_calc_parashar(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &X2ep,
            const unsigned &k1n,
            const unsigned &k1p,
            const unsigned &kep,
            const double &RR0,
            const double &RRp) {
      
  Matrix<double> R3_out(1,2);
  double R3_a_prep=0;
  double R3_b_prep=0;
  const unsigned isupp=min(X1p,kep-1);

  for (unsigned i=k1p;i<=isupp;i++){
    const double id=i;
    const double kepid=kep-id-1;
    R3_a_prep=R3_a_prep+(1-pbinom(kepid,X2ep,RR0))*dbinom(id,X1p,RR0);
    R3_b_prep=R3_b_prep+(1-pbinom(kepid,X2ep,RRp))*dbinom(id,X1p,RRp);
  }
  
  const double k1n1d=k1n-1;
  const double kep1d=kep-1;;

  R3_out[0]=pbinom(k1n1d,X1n,RR0)*((1-pbinom(kep1d,X1p,RR0))+R3_a_prep);
  R3_out[1]=pbinom(k1n1d,X1n,RR0)*((1-pbinom(kep1d,X1p,RRp))+R3_b_prep);
  
  return(R3_out);
}

double EN1_sum_calc_parashar(const unsigned &X1n,
            const unsigned &k1n,
            const unsigned &kn,
            const double &RR0) {
  
  double EN1_sum_out=0;

  for (unsigned k=k1n;k<kn;k++){
    const double kd=k;
    EN1_sum_out=EN1_sum_out+dbinom(kd,X1n,RR0);
  }
  return(EN1_sum_out);
}

double EN2_sum_calc_parashar(const unsigned &X1p,
            const unsigned &k1p,
            const unsigned &kep,
            const double &RR0) {
  
  double EN2_sum_out=0;

  for (unsigned l=k1p;l<kep;l++){
    const double ld=l;
    EN2_sum_out=EN2_sum_out+dbinom(ld,X1p,RR0);
  }
  return(EN2_sum_out);
}





