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
#
./OpenCL-Wrapper.x > out.$SLURM_JOBID.dat

echo "That's all folks!!!"



