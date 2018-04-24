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
