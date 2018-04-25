#!/bin/bash
#QSUB -q  gr10331b
#QSUB -ug gr10331
#QSUB -W 12:0
#QSUB -oo std.out
#QSUB -eo std.err
#QSUB -M john@cheme.kyoto-u.ac.jp
#QSUB -A p=1:t=1:c=1:m=2G
#QSUB -ug gr10331
#QSUB -m be
#================ PBS Options ====================

#==================Shell Script===================
cd $QSUB_WORKDIR
set -x

# automatically
./diffusion2d.sopt
