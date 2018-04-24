OMP  = NO
EXE  = sopt
ifeq ($(ENV), GCC)
     CC        = gcc-7
     CXX       = g++-7
     FFT_DIR   = /opt/fftw-3.3.4
     HDF_DIR   = /opt/hdf5
     SILO_DIR  = /opt/silo

     ifeq ($(OMP), YES)
          EXE      = popt
          CFLAGS   = -fopenmp 
          FFTLINKS = -L$(FFT_DIR)/lib -lfftw3_omp -lfftw3
     else
          FFTLINKS = -L$(FFT_DIR)/lib -lfftw3
     endif
     CFLAGS  += -Ofast -Wall -Wno-unknown-pragmas -Wfatal-errors -Wmisleading-indentation\
              -DWITHOUT_MPI -D__ArrayExtensions -DNDEBUG -lc++ -lz \
              -fomit-frame-pointer -fstrict-aliasing -Wstrict-aliasing=2 -ffast-math \
              -msse2 -mfpmath=sse \
	      -I$(FFT_DIR)/include -I$(HDF_DIR)/include -I$(SILO_DIR)/include
     CXXFLAGS = -std=c++11 $(CFLAGS)
     LINKS = $(FFTLINKS) -L$(SILO_DIR)/lib -L$(HDF_DIR)/lib -lsiloh5 -lhdf5_hl -lhdf5 -lm
endif
ifeq ($(ENV), CLANG)
     CC       = clang
     CXX      = clang++
     FFT_DIR  = /usr/local
     HDF_DIR  = /usr/local
     SILO_DIR = /usr/local

     FFTLINKS = -L$(FFT_DIR)/lib -lfftw3
     CFLAGS  += -Ofast -Wall -Wno-unknown-pragmas -Wfatal-errors -Wno-c++1z-extensions\
              -DWITHOUT_MPI -D__ArrayExtensions -DNDEBUG \
              -fomit-frame-pointer -fstrict-aliasing -Wstrict-aliasing=2 -ffast-math \
              -msse2 -mfpmath=sse -Wno-deprecated-declarations \
	      -I/usr/local/include
     CXXFLAGS = -std=c++11 --stdlib=libc++ $(CFLAGS)
     LINKS = $(FFTLINKS) -L$(SILO_DIR)/lib -L$(HDF_DIR)/lib -lsiloh5 -lhdf5_hl -lhdf5 -lm
endif
ifeq ($(ENV), ICC)
     CC        = icc
     CXX       = icpc
     FFT_DIR   = /usr/local/fftw-3.3.4
     VISIT_DIR = /opt/visit/thirdparty/2.10.2/shared/icc
     HDF_DIR   = $(VISIT_DIR)/hdf5/1.8.14/linux-x86_64_icc
     SILO_DIR  = $(VISIT_DIR)/silo/4.10.1/linux-x86_64_icc

     ifeq ($(OMP), YES)
         EXE      = popt
         CFLAGS   = -qopenmp -parallel 
         FFTLINKS = -L$(FFT_DIR)/lib -lfftw3_omp -lfftw3
     else
         FFTLINKS = -L$(FFT_DIR)/lib -lfftw3
     endif
     CFLAGS += -O3 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 \
              -Wall -DWITHOUT_MPI -D__ArrayExtensions -DNDEBUG \
              -fomit-frame-pointer \
	      -I$(FFT_DIR)/include -I$(HDF_DIR)/include -I$(SILO_DIR)/include
     CXXFLAGS = -std=c++11 $(CFLAGS)
     LINKS = $(FFTLINKS) -L$(SILO_DIR)/lib -L$(HDF_DIR)/lib -lsiloh5 -lhdf5_hl -lhdf5 -lm
endif
ifeq ($(ENV), KUDPC)
     LDIR      = /LARGE0/gr10331
     CC        = icc
     CXX       = icpc
     FFT_DIR   = /opt/app/fftw/3.3.5/intel-17.0-impi-2017.1
     HDF_DIR   = /opt/app/hdf5/1.8.17/intel-17.0
     SILO_DIR  = $(LDIR)/app/silo/4.10.2/intel-17.0

     ifeq ($(OMP), YES)
         EXE      = popt
         CFLAGS   = -qopenmp -parallel 
         FFTLINKS = -L$(FFT_DIR)/lib -lfftw3_omp -lfftw3
     else
         FFTLINKS = -L$(FFT_DIR)/lib -lfftw3
     endif
     CFLAGS += -xavx2 -O3 \
              -Wall -D__ArrayExtensions -DNDEBUG \
              -fomit-frame-pointer\
	      -I$(FFT_DIR)/include -I$(HDF_DIR)/include -I$(SILO_DIR)/include
     CXXFLAGS = -std=c++11 $(CFLAGS)
     LINKS = $(FFTLINKS) -L$(SILO_DIR)/lib -L$(HDF_DIR)/lib -lsiloh5 -lhdf5_hl -lhdf5 -lm
endif

.SUFFIXES: .c .cc .cxx .o

OBJS_DIFF3 = diffusion3d.o
OBJS_DIFF2 = diffusion2d.o
OBJS_FFT2  = fft2d.o
OBJS = $(OBJS_DIFF3) $(OBJS_DIFF2) $(OBJS_FFT2)

all: diffusion3d diffusion2d fft2d

diffusion3d: $(OBJS_DIFF3)
	$(CXX) $(OBJS_DIFF3) -o $@.$(EXE) $(CFLAGS) $(LINKS)
	$(SIGN)

diffusion2d: $(OBJS_DIFF2)
	$(CXX) $(OBJS_DIFF2) -o $@.$(EXE) $(CFLAGS) $(LINKS)
	$(SIGN)

diffusion2dk: $(OBJS_DIFF2K)
	$(CXX) $(OBJS_DIFF2K) -o $@.$(EXE) $(CFLAGS) $(LINKS)
	$(SIGN)

fft2d: $(OBJS_FFT2)
	$(CXX) $(OBJS_FFT2) -o $@.$(EXE) $(CFLAGS) $(LINKS)
	$(SIGN)


.cxx.o:
	$(CXX) -c $< $(CXXFLAGS) -o $@

.cc.o:
	$(CXX) -c $< $(CXXFLAGS) -o $@

.c.o:
	$(CC) -c $< $(CFLAGS) -o $@

clean:
	rm -f *~ $(OBJS)
purge:
	rm -f *~ $(OBJS) *.sopt *.popt
