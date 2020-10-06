#!/bin/bash -e
#SBATCH -A PROJECT_ID
#SBATCH -J admixv1
#SBATCH --time 1:00:00 #
#SBATCH -N 1
#SBATCH -c 36
#SBATCH -n 1
#SBATCH --mem=5GB
#SBATCH --array=1-10
#SBATCH --partition=large
#SBATCH --mail-user=natalie.forsdick@postgrad.otago.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --output admix_DATE.%j.out # CHANGE number for new run
#SBATCH --error admix_DATE.%j.err #  CHANGE number for new run

# Make directories for the outputs
mkdir /nesi/nobackup/PROJECT_ID/natalie/admixture_${SLURM_ARRAY_TASK_ID}
cd /nesi/nobackup/PROJECT_ID/natalie/admixture_${SLURM_ARRAY_TASK_ID}

# Call admixture program
admixture=/nesi/nobackup/PROJECT_ID/bin/admixture_linux-1.3.0/admixture ## to update

# Designate how many hypothetical populations/groups to test:
for k in {1..10};
do

# Run admixture. -C = designate termination criteria, --cv = 10-fold cross-validation. 
# -s = seed, -j = number of threads (multithreaded mode).
$admixture -C 0.0001 --cv=10 ../INPUT_FILE.bed $k \
-s time -j24 | tee log_FILE_NAME_${k}.out 

done
