#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrix.h" 

using namespace scythe;
using namespace std;

Matrix<double> R1_calc_jones(const unsigned &X1n,          
            const unsigned &e2n,
            const unsigned &k1n,
            const unsigned &kn,            
            const double &RR0,
            const double &RRn);
            
Matrix<double> R2_calc_jones(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &e2n,
            const unsigned &e2p,
            const unsigned &k1n,
            const unsigned &kp,
            const unsigned &kn,
            const double &RR0,
            const double &RRp);
            
Matrix<double> R3_calc_jones(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &X2p,
            const unsigned &k1n,
            const unsigned &k1p,
            const unsigned &k2p,
            const double &RR0,
            const double &RRp);
            

double PET_calc_jones(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &k1n,
            const unsigned &k1p,
            const double &RR0);
            
double EN1_sum_calc_jones(const unsigned &X1n,
            const unsigned &k1n,
            const unsigned &X1p,
            const double &RR0);
            
double EN2_sum_calc_jones(const unsigned &X1p,
            const unsigned &k1p,
            const unsigned &X1n,
            const unsigned &k1n,
            const double &RR0);

Matrix<double> R1_calc_parashar(const unsigned &X1n,          
            const unsigned &X2n,
            const unsigned &k1n,
            const unsigned &kn,            
            const double &RR0,
            const double &RRn);
            
Matrix<double> R2_calc_parashar(const unsigned &X1n,
            const unsigned &X2n,
            const unsigned &Xp,
            const unsigned &k1n,
            const unsigned &kp,
            const unsigned &kn,
            const double &RR0,
            const double &RRp);
            
Matrix<double> R3_calc_parashar(const unsigned &X1n,
            const unsigned &X1p,
            const unsigned &X2p,
            const unsigned &k1n,
            const unsigned &k1p,
            const unsigned &kep,
            const double &RR0,
            const double &RRp);
            
double EN1_sum_calc_parashar(const unsigned &X1n,
            const unsigned &k1n,
            const unsigned &kn,
            const double &RR0);

double EN2_sum_calc_parashar(const unsigned &X1p,
            const unsigned &k1p,
            const unsigned &kep,
            const double &RR0);
