
#!/bin/bash

#PBS -N ERFHOUR
#PBS -l ncpus=8
#PBS -l mem=10GB
#PBS -J 1-829
#PBS -o logs_hour
#PBS -j oe

cd ${PBS_O_WORKDIR}
apptainer run image.sif estimate.R

