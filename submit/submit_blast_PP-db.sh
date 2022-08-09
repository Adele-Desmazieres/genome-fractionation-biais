#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=blast_prunus

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
#SBATCH --mail-type=ALL

### Output
#SBATCH --output=out_blast_prunus_%A.log
###############################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

FASTA='../data/Malus-domestica-proteome.fasta'
NAME="malus_vs_prunus"
DB="../database/PP-db"

OUT="$NAME.blastp"

blastp \
 -query $FASTA -db $DB -out $OUT \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 35 -outfmt 5
 
echo "////////// job done ///////////"
