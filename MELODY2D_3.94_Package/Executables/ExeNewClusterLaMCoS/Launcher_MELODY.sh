#!/bin/bash
#SBATCH -J MELODY_MySim
#SBATCH --account=MyProject
#SBATCH --partition=medium
#SBATCH --time='48:00:00'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=4000
#SBATCH --mail-type=end
#SBATCH --mail-user=surname.name@insa-lyon.fr
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt

echo "Starting at `date`"
echo "Running on hosts : $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."

export OMP_NUM_THREAD=24
./MELODY_2D_X.XX 0