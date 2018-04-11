#!/bin/bash

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --job-name=adsorption
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --time=02:00:00
#SBATCH --output=slurm.%N.%j.out
#SBATCH --error=slurm.%N.%j.err
#SBATCH --export=ALL

cd  /home/hw412/adsorption_in_slit_pore_1DMC_v1/prepare/
srun /home/hw412/adsorption_in_slit_pore_1DMC_v1/prepare/a.out < input_pressure.dat 
