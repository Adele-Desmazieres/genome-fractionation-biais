#!/bin/sh 


################################ Slurm options #################################

### Job name
#SBATCH --job-name=iadhore

### Limit run time "days-hours:minutes:seconds"
##SBATCH --time=15:00:00

### Requirements
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --mem-per-cpu=15
#SBATCH --cpus-per-task=35

### Email
#SBATCH --mail-user=adele.desmazieres@inrae.fr
#SBATCH --mail-type=NONE

### Output
#SBATCH --output=log/%A_out_iadhore.log
################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

i-adhore scripts/iadhore/iadhore.ini

echo "////////// job done ///////////"
