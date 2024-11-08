#!/bin/bash

#PBS -N ERF
#PBS -l ncpus=8
#PBS -l mem=16GB
#PBS -J 5000-5020
#PBS -o logs
#PBS -j oe

cd ${PBS_O_WORKDIR}

apptainer run image.sif estimate.R
