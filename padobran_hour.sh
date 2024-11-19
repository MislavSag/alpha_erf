
#!/bin/bash

#PBS -N ERFHOUR
#PBS -l ncpus=8
#PBS -l mem=8GB
#PBS -l walltime=120:00:00
#PBS -J 1-829
#PBS -o logs_hour
#PBS -j oe

cd ${PBS_O_WORKDIR}
apptainer run image.sif estimate_hour.R

