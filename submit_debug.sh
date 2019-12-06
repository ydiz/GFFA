#!/bin/bash
#COBALT -t 00:30:00
#COBALT -n 4
#COBALT --attrs mcdram=cache:numa=quad
#COBALT -A CSC249ADSE03
#COBALT -o output-$COBALT_JOBID.log
###COBALT -e ./errors/$COBALT_JOBID.err
#COBALT -M yz3210@columbia.edu
#COBALT --jobname GFFA_theta
#COBALT -q debug-flat-quad

# mkdir -p logs
# mkdir -p errors

echo "Starting Cobalt job script"
export n_nodes=4 # $COBALT_JOBSIZE
export n_mpi_ranks_per_node=1
export n_mpi_ranks=$(($n_nodes * $n_mpi_ranks_per_node))

export n_openmp_threads_per_rank=64 # zyd: this number should be equal to number of cores or cores * hyperthreads per core
export n_hyperthreads_per_core=1
export n_hyperthreads_skipped_between_ranks=64 # zyd: this number should be equal to n_openmp_threads_per_rank 

source ~/theta_modules.sh

program=Test_gauge_force

aprun -n $n_mpi_ranks -N $n_mpi_ranks_per_node -d $n_hyperthreads_skipped_between_ranks -j $n_hyperthreads_per_core -cc depth -d $n_hyperthreads_skipped_between_ranks --env OMP_NUM_THREADS=$n_openmp_threads_per_rank  \
  /projects/CSC249ADSE03/yidizhao/GFFA/$program /lus/theta-fs0/projects/CSC249ADSE03/yidizhao/GFFA/Oct25/beta10_M0.25_traj1/ckpoint_lat.5000 --mpi 1.1.2.2 --grid 16.16.16.16 > output-${COBALT_JOBID}.log 2>&1 &

