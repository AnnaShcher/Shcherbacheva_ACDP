#!/bin/bash -l
# created: Mar 20, 2019 12:10 PM
# author: shcherba
#SBATCH -J MCMC_ACDC
#SBATCH --constraint="snb|hsw"
#SBATCH -o std.out
#SBATCH -e std.err
#SBATCH -p serial
#SBATCH -n 1
#SBATCH -t 01:20:00
#SBATCH --mem-per-cpu=64000
#SBATCH --mail-type=END
#SBATCH --mail-user=anna.shcherbacheva@helsinki.fi

# commands to manage the batch script
#   submission command
#     sbatch [script-file]
#   status command
#     squeue -u shcherba
#   termination command
#     scancel [jobid]

# For more information
#   man sbatch
#   more examples in Taito guide in
#   http://research.csc.fi/taito-user-guide

# example run commands
# srun ./my_serial_program
./mcmcrun

# This script will print some usage statistics to the
# end of file: std.out
# Use that to improve your resource request estimate
# on later jobs.
used_slurm_resources.bash
