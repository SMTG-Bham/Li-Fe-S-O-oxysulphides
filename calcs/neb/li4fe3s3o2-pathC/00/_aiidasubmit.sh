#!/bin/bash
#SBATCH --no-requeue
#SBATCH --job-name="aiida-568927"
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --partition=standard
#SBATCH --account=e05-power-dos
#SBATCH --qos=standard
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:00:00


module load epcc-job-env
export OMP_NUM_THREADS=1

'srun' '-n' '256' '--cpu-bind=rank' '--hint=nomultithread' '--distribution=block:block' '/work/e05/e05/bz1/vasp/6.2.0/bin/vasp_std'  > 'vasp_output' 2>&1
