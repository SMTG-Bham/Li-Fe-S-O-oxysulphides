#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -m n
#$ -N aiida-551524
#$ -V
#$ -o _scheduler-stdout.txt
#$ -e _scheduler-stderr.txt
#$ -pe mpi 20
#$ -l h_rt=36:00:00

# Nothing

module unload default-modules/2018
module load rcps-core/1.0.0
module unload compilers mpi
module load compilers/intel/2017/update1
module load mpi/intel/2017/update1/intel

module list

'mpirun' '-np' '20' '/shared/ucl/apps/vasp/5.4.4-18apr2017/intel-2017/bin/vasp_std'  > 'vasp_output' 2>&1

# Nothing
