#!/bin/bash -l
#$ -cwd
#$ -S /bin/bash
#$ -m n
#$ -N aiida-548912
#$ -o _scheduler-stdout.txt
#$ -e _scheduler-stderr.txt
#$ -pe mpi 136
#$ -l h_rt=24:00:00

#$ -P Gold
#$ -A Faraday_FCAT

#$ -ac allow=A

'gerun' '/home/uccabz1/vasp/6.2.0/bin-a-node/vasp_std'  > 'vasp_output' 2>&1
