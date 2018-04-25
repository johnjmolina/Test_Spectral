Test code showing how to solve diffusion equation using spectral method with the fft.
Uses fftw and silo/hdf5 libraries.

- diffusion2d : solve 2d diffusion equation
- diffusion3d : solve 3d diffusion equation
- fft2d	      : sample code showing how to define forward/backward transforms

Output Silo files can be read by VisIt/Paraview.
They can also be imported/exported using Python.
For this we need to the pyublas/pyvisfile packages.

How to install:
# install silo using homebrew
  > brew install silo
# create separate conda environment (just to be safe)
  > conda create --name silo --channel conda-forge python=3.4 boost h5py numpy matplotlib ipython-notebook
  > source active silo
# install pyublas
  > git clone https://github.com/inducer/pyublas.git
  > cd pyublas
  > ./install_pyublas.sh   [edit to set variables for your system]
# install pyvisfile
  > git clone https://github.com/inducer/pyvisfile.git
  > cd pyvisfile
  > ./install_pyvisfile.sh [edit to set variables for your system]
