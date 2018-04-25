#ifdef _OPENMP
#include <omp.h>
#endif

#include <complex>
#include <silo.h>
#include <fftw3.h>
#include "Array.h"
#include "utils.hxx"

using Array::array1;
using Array::array2;
using Array::array3;
using Array::array4;
typedef std::complex<double> Complex;

const static int DIM = 3;  
const static size_t align = sizeof(Complex);
const static int BUFFER_SIZE = 256;

//fft variables
int Nx, Ny, Nz, HNz, Ksize;
double fftnorm;

fftw_plan fw, bw;
array1<double>  kx, ky, kz;
array3<double>  K2;
array3<Complex> xout, hL;

//silo variables
array1<double> xval, yval, zval;

//system vars
double dt;
double dx;
const static double Drho = 1.0;

inline void fft(const array3<double> &fr, const array3<Complex> &fk){
  fftw_execute_dft_r2c(fw, &fr(1,1,1), reinterpret_cast<fftw_complex*>(fk()));
}
inline void ifft(const array3<Complex> &fk, const array3<double> &fr){
  fftw_execute_dft_c2r(bw, reinterpret_cast<fftw_complex*>(fk()), &fr(1,1,1));
}
inline int fftfreq(const int&i, const int&n){return (i >= n/2 ? -(n-i) : i);}

inline void PeriodicBC(array3<double> &X){
  // bottom and upper
#pragma omp parallel for
  for(auto i=1;i<=Nx;i++){
    for(auto j=1;j<=Ny;j++){
      X(i,j,0   ) = X(i,j,Nz);
      X(i,j,Nz+1) = X(i,j,1 );
    }
  }
  // right and left
#pragma omp parallel for  
  for(auto j=1;j<=Ny  ;j++){
    for(auto k=0;k<=Nz+1 ;k++){
      X(0   ,j,k) = X(Nx,j,k);
      X(Nx+1,j,k) = X(1 ,j,k);
    }
  }

  // front and back
#pragma omp parallel for    
  for(auto i=0;i<=Nx+1;i++){
    for(auto k=0;k<=Nz+1;k++){
      X(i,0   ,k) = X(i,Ny,k);
      X(i,Ny+1,k) = X(i,1 ,k);
    }
  }
}

inline void initialize(const int n1, const int n2, const int n3, const double _dx = 1.0){
  dx = _dx;
  Nx = 1 << n1;
  Ny = 1 << n2;
  Nz = 1 << n3;
  HNz = (Nz/2 + 1);
  Ksize = Nx*Ny*HNz;
  fftnorm = 1.0 / static_cast<double>(Nx*Ny*Nz);
  
  array3<double>  xin(Nx+2,Ny+2,Nz+2, align);    // real space variable
  xout.Allocate(Nx, Ny, HNz, align);             // reciprocal space variable (global memory to perform ffts)
#ifdef _OPENMP
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  { // fft
    fftw_iodim64 dims[DIM] = {{.n = Nx, .is = (Ny+2)*(Nz+2), .os = Ny*HNz},
			    {.n = Ny, .is = (Nz+2),        .os = HNz},
			    {.n = Nz, .is = 1,             .os = 1}};
    fw = fftw_plan_guru64_dft_r2c(DIM, dims,
				  0, (fftwf_iodim64 *) 0,
				  &xin(1,1,1), (fftw_complex*)(xout()), FFTW_MEASURE);
  }
  {
    fftw_iodim64 dims[DIM] = {{.n = Nx, .is = Ny*HNz, .os = (Ny+2)*(Nz+2)},
			    {.n = Ny, .is = HNz,    .os = (Nz+2)},
			    {.n = Nz, .is = 1,      .os = 1}};
    bw = fftw_plan_guru64_dft_c2r(DIM, dims,
				  0, (fftwf_iodim64 *) 0,
				  (fftw_complex*)(xout()), &xin(1,1,1), FFTW_MEASURE);
  }

  { // real-space variables
    xval.Allocate(Nx+2, align);
    yval.Allocate(Ny+2, align);
    zval.Allocate(Nz+2, align);
    for(auto i=0; i < Nx+2; i++) xval(i) = i*dx;
    for(auto j=0; j < Ny+2; j++) yval(j) = j*dx;
    for(auto k=0; k < Nz+2; k++) zval(k) = k*dx;
  }

  { // reciprocal-space variables
    kx.Allocate(Nx, align);
    ky.Allocate(Ny, align);
    kz.Allocate(HNz, align);
    K2.Allocate(Nx, Ny, HNz, align);
    hL.Allocate(Nx, Ny, HNz, align);
    double twopi     = 2.0*M_PI;
    double length[DIM] = {Nx*dx, Ny*dx, Nz*dx};
    double deltak[DIM] = {twopi/length[0], twopi/length[1], twopi/length[2]};
    double maxk2     = 0.0;
    for(auto i = 0; i < Nx; i++){
      double qx = kx(i) = fftfreq(i,Nx)*deltak[0];
      for(auto j = 0; j < Ny; j++){
	double qy = ky(j) = fftfreq(j,Ny)*deltak[1];
	for(auto k = 0; k < HNz; k++){
	  double qz = kz(k) = k*deltak[2];

	  K2(i,j,k) = qx*qx + qy*qy + qz*qz;
	  maxk2 = MAX(K2(i,j,k), maxk2);
	}
      }
    }
    dt = 1.0 / (Drho*maxk2);
    double hDrho = -Drho*dt;    
    for(auto im = 0; im < Ksize; im++) hL(im) = exp(hDrho*K2(im))*fftnorm;
  }
}

void output(double time, const array3<double> &rho){
  static int frame_id = 0;
  static array4<double>  jr;
  static array3<Complex> jk;
  static char* jnames[] = {"j_z", "j_y", "j_x"}; // SILO HACK: silo requires transposed data
  static char* anames[] = {"z", "y", "x"};       // SILO HACK
  if(frame_id == 0){
    jr.Allocate(DIM, Nx+2, Ny+2, Nz+2, align);
    jk.Allocate(Nx, Ny, HNz, align);
  }

  char buffer[BUFFER_SIZE];
  snprintf(buffer, BUFFER_SIZE, "diff3d_%05d.silo", frame_id);
  DBfile *dbfile = DBCreate(buffer, DB_CLOBBER, DB_LOCAL, "diffusion", DB_HDF5);
  int dims[]   = {Nz+2, Ny+2, Nx+2}; // SILO HACK
  DBoptlist *optlist = DBMakeOptlist(5);  
  DBAddOption(optlist, DBOPT_DTIME, &time);
  DBAddOption(optlist, DBOPT_CYCLE, &frame_id);
  DBAddOption(optlist, DBOPT_XLABEL, anames[0]);
  DBAddOption(optlist, DBOPT_YLABEL, anames[1]);
  DBAddOption(optlist, DBOPT_ZLABEL, anames[2]);
  { //Grid
    double *coords[] = {zval(), yval(), xval()}; // SILO HACK
    DBPutQuadmesh(dbfile, "grid", NULL, coords, dims, DIM, DB_DOUBLE, DB_COLLINEAR, optlist);
  }
  {// scalar variables
    DBPutQuadvar1(dbfile, "rho", "grid", rho(), dims, DIM, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);        
  }
  {// vector variables
    fft(rho, xout);
    auto inorm = -Complex(0.0, 1.0)*fftnorm;

    for(auto i = 0; i < Nx; i++) for(auto j = 0; j < Ny; j++) for(auto k = 0; k < HNz; k++)
								jk(i,j,k) = inorm*kx(i)*xout(i,j,k);
    auto jrx = jr[0];    
    ifft(jk, jrx);
    PeriodicBC(jrx);    
    
    for(auto i = 0; i < Nx; i++) for(auto j = 0; j < Ny; j++) for(auto k = 0; k < HNz; k++)
								jk(i,j,k) = inorm*ky(j)*xout(i,j,k);
    auto jry = jr[1];    
    ifft(jk, jry);
    PeriodicBC(jry);

    for(auto i = 0; i < Nx; i++) for(auto j = 0; j < Ny; j++) for(auto k = 0; k < HNz; k++)
								jk(i,j,k) = inorm*kz(k)*xout(i,j,k);
    auto jrz = jr[2];    
    ifft(jk, jrz);
    PeriodicBC(jrz);

    double* ptr[DIM] = {jrz(), jry(), jrx()}; // SILO HACK
    DBPutQuadvar(dbfile, "flux", "grid", DIM, jnames, ptr, dims, DIM, NULL, 0, DB_DOUBLE, DB_NODECENT, optlist);
  }
  DBFreeOptlist(optlist);
  DBClose(dbfile);
  frame_id++;
}

inline void update(array3<double> &rho){
  fft(rho, xout);

#pragma omp parallel for
  for(auto im = 0; im < Ksize; im++) xout(im) *= hL(im);

  ifft(xout, rho);
  PeriodicBC(rho);
}

int main(int arg, char** argv){
  const int FRAMES = 100;
  const int GTS    = 20;
  init_threads();
  initialize(7, 7, 7);

  const int r0     = MIN(MIN(Nx,Ny), Nz)/8;  
  array3<double> rho(Nx+2,Ny+2,Nz+2,align);
  rho = 0.0;
#pragma omp parallel for
  for(auto i = 1; i <= Nx; i++) for(auto j = 1; j <= Ny; j++) for(auto k = 1; k <= Nz; k++)
								if(sqrt(SQ(i-Nx/2) + SQ(j-Ny/2) + SQ(k-Nz/2)) < r0) rho(i,j,k) = 1.0;
#pragma omp parallel for
  for(auto i = 1; i <= Nx; i++) for(auto j = 1; j <= Ny; j++) for(auto k = 1; k <= Nz; k++)
								if(ABS(i-Nx/6) < r0 && ABS(j-Ny/6) < r0 && ABS(k-Nz/6) < r0) rho(i,j,k) += 1.0;
#pragma omp parallel for
  for(auto i = 1; i <= Nx; i++) for(auto j = 1; j <= Ny; j++) for(auto k = 1; k <= Nz; k++)
								if(ABS(i-5*Nx/6) < r0 && ABS(j-5*Ny/6) < r0 && ABS(k-5*Nz/6) < r0) rho(i,j,k) += 1.0;
  
  double time = 0.0;
  output(time, rho);
  for(auto frame = 1; frame <= FRAMES; frame++){
    fprintf(stderr, "# frame = %8d\n", frame);
    for(auto gts = 0; gts < GTS; gts++) update(rho);
    time += GTS*dt;
    output(time, rho);
  }
  return 0;
}
