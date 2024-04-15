#!/bin/bash
# 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p g100_usr_interactive
#SBATCH -A tra23_cfd
#SBATCH -t 00:10:00


echo $HOSTNAME > hostname.dat
for size in 1024 2048 4096 8192; do
    echo $size > size.dat
    for ver in 0 1; do
        ./mm.$ver.x < size.dat > out.$size.$ver.$SLURM_JOBID.dat
    done
done

echo "That's all Folks!!!"



