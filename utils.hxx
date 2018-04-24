#ifndef UTILS_HXX_
#define UTILS_HXX_
#include <cstdio>

inline void init_threads(){
  using std::cerr;
  using std::endl;
#ifdef _OPENMP
  int nthreads, tid, procs, maxt, inpar, dynamic, nested;
#pragma omp parallel private(nthreads, tid)
  {
    tid = omp_get_thread_num();
    if(tid == 0){
      procs   = omp_get_thread_num();
      nthreads= omp_get_num_threads();
      maxt    = omp_get_max_threads();
      inpar   = omp_get_dynamic();
      nested  = omp_get_nested();

      cerr << "# " << endl;
      cerr << "# OMP RUNTIME :" << endl;
      cerr << "# Number of processors        = " << procs    << endl;
      cerr << "# Number of threads           = " << nthreads << endl;
      cerr << "# Max threads                 = " << maxt     << endl;
      cerr << "# In parallel ?               = " << inpar    << endl;
      cerr << "# Dynamic thread enabled ?    = " << dynamic  << endl;
      cerr << "# Nested parallelism enabled? = " << nested   << endl;
      cerr << "# " << endl;
    }
  }
#endif
}
template<typename T> inline T SQ(const T &a){ return a*a ;}
template<typename T> inline T MIN(const T &a, const T &b){ return ( a <= b ? a : b);}
template<typename T> inline T MAX(const T &a, const T &b){ return ( a >= b ? a : b);}
template<typename T> inline T ABS(const T &a){return ( a >= 0 ? a : -a);}

#endif
