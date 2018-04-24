#ifdef _OPENMP
#include <omp.h>
#endif

#include <complex>
#include <silo.h>
#include <fftw3.h>
#include "Array.h"

using Array::array1;
using Array::array2;
using Array::array3;

typedef std::complex<double> Complex;

const static int DIM = 2;  
const static size_t align = sizeof(Complex);

//fft variables
int Nx, Ny, HNy, Ksize;
double fftnorm;

fftw_plan fw_pbc, bw_pbc;
fftw_plan fw_nopbc, bw_nopbc;
array1<double>  kx, ky;
array2<double>  K2;
array2<Complex> xout, xout2;

//silo variables
array1<double> xval, yval;

//system vars
double dt;
double dx;

template<typename T> inline T SQ(const T &a){ return a*a ;}
template<typename T> inline T MIN(const T &a, const T &b){ return ( a <= b ? a : b);}
template<typename T> inline T MAX(const T &a, const T &b){ return ( a >= b ? a : b);}
template<typename T> inline T ABS(const T &a){return ( a >= 0 ? a : -a);}

inline void fft(fftw_plan &fw_plan, double* fr, Complex* fk){
  fftw_execute_dft_r2c(fw_plan, fr, reinterpret_cast<fftw_complex*>(fk));
}
inline void ifft(fftw_plan &bw_plan, Complex *fk, double* fr){
  fftw_execute_dft_c2r(bw_plan, reinterpret_cast<fftw_complex*>(fk), fr);
}
inline int fftfreq(const int&i, const int&n){return (i >= n/2 ? -(n-i) : i);}

inline void PeriodicBC(array2<double> &X){
  // bottom and upper
#pragma omp parallel for
  for(auto i=1;i<=Nx;i++){
    X(i,0   ) = X(i,Ny);
    X(i,Ny+1) = X(i,1 );
  }
  // right and left
#pragma omp parallel for  
  for(auto j=1;j<=Ny  ;j++){
    X(0   ,j) = X(Nx,j);
    X(Nx+1,j) = X(1 ,j);
  }
}

inline void initialize(const int n1, const int n2, const double _dx = 1.0){
  dx = _dx;
  Nx = 1 << n1;
  Ny = 1 << n2;
  HNy = (Ny/2 + 1);
  Ksize = Nx*HNy;
  fftnorm = 1.0 / static_cast<double>(Nx*Ny);

  xout.Allocate(Nx, HNy, align);            // reciprocal space variable (global memory to perform ffts)
  xout2.Allocate(Nx, HNy, align);
#ifdef _OPENMP
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  array2<double>  xin(Nx+2,Ny+2, align);    // real space variable      
  { // pbc plans
    {
      fftw_iodim64 dims[DIM] = {{.n = Nx, .is = (Ny+2), .os = HNy},
				{.n = Ny, .is = 1,      .os = 1}};
      fw_pbc = fftw_plan_guru64_dft_r2c(DIM, dims,
					0, (fftwf_iodim64 *) 0,
					&xin(1,1), (fftw_complex*)(xout()), FFTW_MEASURE);
    }
    {
      fftw_iodim64 dims[DIM] = {{.n = Nx, .is = HNy, .os = (Ny+2)},
				{.n = Ny, .is = 1,   .os = 1}};
      bw_pbc = fftw_plan_guru64_dft_c2r(DIM, dims,
					0, (fftwf_iodim64 *) 0,
					(fftw_complex*)(xout()), &xin(1,1), FFTW_MEASURE);
    }
  }
  { // no pbc plans
    int dims[DIM] = {Nx, Ny};
    fw_nopbc = fftw_plan_dft_r2c(DIM, dims, xin(), (fftw_complex*)(xout()), FFTW_MEASURE);
    bw_nopbc = fftw_plan_dft_c2r(DIM, dims, (fftw_complex*)(xout()), xin(), FFTW_MEASURE);
  }
  

  { // real-space variables
    xval.Allocate(Nx+2, align);
    yval.Allocate(Ny+2, align);
    for(auto i=0; i < Nx+2; i++) xval(i) = i*dx;
    for(auto j=0; j < Ny+2; j++) yval(j) = j*dx;
  }

  { // reciprocal-space variables
    kx.Allocate(Nx, align);
    ky.Allocate(HNy, align);
    K2.Allocate(Nx, HNy, align);
    double twopi     = 2.0*M_PI;
    double length[DIM] = {Nx*dx, Ny*dx};
    double deltak[DIM] = {twopi/length[0], twopi/length[1]};
    for(auto i = 0; i < Nx; i++){
      double qx = kx(i) = fftfreq(i,Nx)*deltak[0];
      for(auto j = 0; j < HNy; j++){
	double qy = ky(j) = j*deltak[1];
	
	K2(i,j) = qx*qx + qy*qy;
      }
    }
  }
}


int main(int arg, char** argv){

  initialize(6, 7);
  array2<double> rho_pbc(Nx+2,Ny+2,align);
  array2<double> rho_nopbc(Nx,Ny,align);
  rho_pbc   = 0.0;
  rho_nopbc = 0.0;
#pragma omp parallel for
  for(auto i = 1; i <= Nx; i++) for(auto j = 1; j <= Ny; j++)
				  if(sqrt(SQ(i-Nx/2) + SQ(j-Ny/2)) < 10) rho_pbc(i,j) = rho_nopbc(i-1,j-1) = 1.0;
  PeriodicBC(rho_pbc);

  fft(fw_pbc, &rho_pbc(1,1), xout());
  fft(fw_nopbc, rho_nopbc(), xout2());
  for(auto i = 0; i < Nx; i++) for(auto j = 0; j < HNy; j++)
				 fprintf(stderr, "%20.16lf + i %20.16lf\n", real(xout(i,j)), imag(xout(i,j)));
  for(auto i = 0; i < Nx; i++) for(auto j = 0; j < HNy; j++)
				 fprintf(stdout, "%20.16lf + i %20.16lf\n", real(xout(i,j)), imag(xout(i,j)));				 

  return 0;
}
