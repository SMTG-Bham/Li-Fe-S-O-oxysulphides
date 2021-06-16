#!/bin/bash
#SBATCH --no-requeue
#SBATCH --job-name="aiida-569017"
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

module restore /etc/cray-pe.d/PrgEnv-gnu
module load cray-fftw
module load cpe/21.03
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH:${LD_LIBRARY_PATH}"
echo "ldd output"
ldd /work/e05/e05/bz1/vtst-code/bin/6.2.0-cpe-21.03/vasp_std 
echo "module list output: "
module list

'srun' '-n' '1280' '--cpu-bind=rank' '--hint=nomultithread' '--distribution=block:block' '/work/e05/e05/bz1/vtst-code/bin/6.2.0-cpe-21.03/vasp_std'  > 'stdout' 2>&1
