#!/bin/bash

################################ Slurm options #################################

### Job name
#SBATCH --job-name=blast

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
#SBATCH --output=log/%A_out_blast_all_vs_all.log
################################################################################

echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
echo '########################################'

#FASTA="data"
#DB="database"
#OUT="results/blast"


blastp \
 -query "$DATA/Malus-domestica-proteome.fasta" -db "$DB/MD-db" -out "$OUT/MD_vs_MD.txt" \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 35 -outfmt 6

blastp \
 -query "$DATA/Prunus-persica-proteome.fasta" -db "$DB/PP-db" -out "$OUT/PP_vs_PP.txt" \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 35 -outfmt 6

blastp \
 -query "$DATA/Malus-domestica-proteome.fasta" -db "$DB/PP-db" -out "$OUT/MD_vs_PP.txt" \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 35 -outfmt 6


echo "////////// job done //////////"
