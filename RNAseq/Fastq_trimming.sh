#!/bin/bash
#SBATCH --partition=kamiak         # Partition/Queue to use
#SBATCH --job-name=Trout_mapping_all         # Job name
#SBATCH --output=myJob_%j.out      # Output file (stdout)
#SBATCH --error=myJob_%j.err       # Error file (stderr)
#SBATCH --mail-type=ALL            # Email notification: BEGIN,END,FAIL,ALL
#SBATCH --mail-user=your.bradford.dimos@wsu.edu  # Email address for notifications
#SBATCH --time=1-8:00:00          # Wall clock time limit Days-HH:MM:SS
#SBATCH --nodes=1                  # Number of nodes (min-max)
#SBATCH --ntasks-per-node=1        # Number of tasks per node (max)
#SBATCH --ntasks=1                 # Number of tasks (processes)
#SBATCH --cpus-per-task=8         # Number of cores per task (threads)
module load trimgalore/0.6.6
samples=`cat sample_list.txt`
for i in $samples;do
        /opt/apps/trimgalore/0.6.6/trim_galore /scratch/user/bradford.dimos/20220414_103438/reads/${i}/${i}.fastq 
done