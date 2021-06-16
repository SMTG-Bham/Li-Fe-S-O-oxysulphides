#!/bin/bash
#SBATCH --no-requeue
#SBATCH --job-name="aiida-567431"
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --partition=standard
#SBATCH --account=e05-power-dos
#SBATCH --qos=standard
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=128
#SBATCH --time=12:00:00


module load epcc-job-env
export OMP_NUM_THREADS=1

'srun' '-n' '1280' '--cpu-bind=rank' '--hint=nomultithread' '--distribution=block:block' '/work/y07/shared/vasp5/vasp.5.4.4.pl2-vtst-gcc10/bin/vasp_std'  > 'stdout' 2>&1
