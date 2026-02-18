#!/bin/bash
#SBATCH --job-name=G500_128C
#SBATCH --account=r250142
#SBATCH --partition=instant
#SBATCH --constraint=x64cpu
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=200G
#SBATCH --output=g500_128cores_%j.out

module purge
module load  openmpi

echo "Lancement Graph500 - 128 coeurs - Scale 27"
# --bind-to core assure que chaque processus reste sur son coeur physique
mpirun -n 128 ./graph500_reference_bfs_sssp 27 16
