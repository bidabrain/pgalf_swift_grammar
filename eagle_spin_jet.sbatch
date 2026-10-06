#!/bin/bash
#SBATCH --job-name=pgalf
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH -N 1 
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --exclude=grammar082
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --output=pgalf_%j.log       # 标准输出+错误都写这里
#SBATCH --error=pgalf_%j.log


module purge
module load gnu14/14.2.0
module load openmpi5/5.0.7
module load hwloc/2.12.0
module load libfabric/1.18.0
module load cmake/3.31.10 
module load gsl/2.8
module load hdf5/1.14.6
module unload ucx
module load ucx/1.18.0-mt
module list

export OMP_NUM_THREADS=1
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/fftw3/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/GKlib/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/METIS/lib
export LIBRARY_PATH=$LIBRARY_PATH:/home/dbi224/local/omp5/ParMETIS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/fftw3/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/GKlib/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/METIS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/dbi224/local/omp5/ParMETIS/lib



#--- 2) 路径 ---------------------------------------------------
SIM=/gpfs/dbi224/swift_eagle_flamingo/eagle_spin_jet
OUT=/gpfs/dbi224/swift_eagle_flamingo/pgalf/eagle_spin_jet
PGALF=/home/dbi224/opt/pgalf_swift

#--- 3) 跑流水线 ----------------------------------
# SLURM 上用 srun 启动 MPI (比 mpirun 稳); 若 srun 报 PMI/PMIx 错, 换回默认 mpirun
cd "$PGALF"
./run_pgalf_eagle.sh -i "$SIM" -o "$OUT" -s all -n 32 -N 32 -z 15 -m "srun -n"
