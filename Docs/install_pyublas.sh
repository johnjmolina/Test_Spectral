#!/bin/bash
export MACOSX_DEPLOYMENT_TARGET=10.10
export CONDA=/opt/anaconda/envs/silo
./configure.py \
    --python-exe=$CONDA/bin/python \
    --boost-inc-dir=$CONDA/include \
    --boost-lib-dir=$CONDA/lib \
    --boost-compiler=clang \
    --boost-python-libname=boost_python3 \
    --cxxflags="-stdlib=libc++"
make
make install

