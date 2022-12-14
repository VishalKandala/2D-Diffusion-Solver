#!/bin/sh
#ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=2diff       #Set the job name to "Example1"
#SBATCH --time=0-03:00:00          #Set the wall clock limit to 0 Days and 5hrs
#SBATCH --ntasks=25                   #Request 50 tasks
#SBATCH --ntasks-per-node=25          #Request 245tasks/cores per node
#SBATCH --mem=4G                  #Request 4096MB (4GB) per node 
#SBATCH --output=2diff.%j      #Send stdout/err to "Example3Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --account=122759699519             #Set billing account
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=vishalkandala@tamu.edu    #Send all emails to email_address 

cd /scratch/user/vishalkandala/PETSC/2DIFF/

./strongscale
