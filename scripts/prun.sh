#!/bin/bash
#================ LSF Options ====================
#QSUB -q gr10331d
#QSUB -ug gr10331
#QSUB -W 12:0
#QSUB -oo std.out
#QSUB -eo std.err
#QSUB -u john@cheme.kyoto-u.ac.jp
#QSUB -B
#QSUB -N
#QSUB -A p=1:t=28:c=28:m=2G

#==================Shell Script===================
set -x
# automatically
export OMP_NUM_THREADS=$LSB_THREADS
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=STATIC

export KMP_AFFINITY="verbose,compact,1"
export MKL_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS="MKL_ALL=1, MKL_FFT=$OMP_NUM_THREADS"
export MKL_DYNAMIC=FALSE

aprun -n $LSB_PROCS -N $LSB_PPN -d $LSB_CPUS -cc none ./kapsel -Iinput.udf -Ddefine.udf -Rrestart.udf -Ooutput.udf
