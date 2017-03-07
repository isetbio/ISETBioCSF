#!/bin/bash
# Additional qsub options here . . .
# Make sure we only use 1 matlab instance
#$ -l matlab=1

# Request 8 GB of RAM. This is per processor slot
#$ -l h_vmem=64G

# Request 8 cores, 
#$ -pe onenode 8 

#$ -M cottaris@upenn.edu

#$ -N P&W96Job1 

#$ -e /home1/c/cottaris/documents/MATLAB/PoirsonAndWandellJob1.log 

#$ -o /home1/c/cottaris/documents/MATLAB/PoirsonAndWandellJob1.log


# Run the job
matlab -nodesktop < /home1/c/cottaris/documents/MATLAB/runPoirsonAndWandell96ClusterJob1.m

