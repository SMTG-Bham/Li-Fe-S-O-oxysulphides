#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -m n
#$ -N aiida-512731
#$ -o _scheduler-stdout.txt
#$ -e _scheduler-stderr.txt
#$ -pe mpi 189
#$ -l h_rt=24:00:00

module unload default-modules/2018
module load rcps-core/1.0.0

module load compilers/intel
module load mpi/intel

'gerun' '/home/uccabz1/vasp/6.2.0/bin/vasp_std'  > 'vasp_output' 2>&1
