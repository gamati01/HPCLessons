#!/bin/bash
# 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p g100_usr_interactive
#SBATCH -A tra23_cfd
#SBATCH -t 0:10:00


echo $HOSTNAME > hostname.dat
export OMP_NUM_THREADS=1

./mm.x  > out.4096.C.$SLURM_JOBID.dat


