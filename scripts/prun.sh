#!/bin/bash
#QSUB -q  gr10331b
#QSUB -ug gr10331
#QSUB -W  0:10
#QSUB -o  std.out
#QSUB -e  std.err
#QSUB -M  john@cheme.kyoto-u.ac.jp
#QSUB -A  p=1:t=18:c=18:m=2G
#QSUB -ug gr10331
#QSUB -m  be
#================ PBS Options ====================

#==================Shell Script===================
cd $QSUB_WORKDIR
set -x
export LD_LIBRARY_PATH=/LARGE0/gr10331/app/silo/4.10.2/intel-17.0/lib:$LD_LIBRARY_PATH
# automatically
./diffusion2d.popt
