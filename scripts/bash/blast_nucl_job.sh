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
#SBATCH --mem-per-cpu=15
#SBATCH --cpus-per-task=10

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


blastn \
 -query "$DATA/Malus-domestica-${nucleique_flag}.fasta" -db "$DB/MD-db/MD-db" -out "$OUT/MD_vs_MD.txt" \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6

blastn \
 -query "$DATA/Prunus-persica-${nucleique_flag}.fasta" -db "$DB/PP-db/PP-db" -out "$OUT/PP_vs_PP.txt" \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6

blastn \
 -query "$DATA/Malus-domestica-${nucleique_flag}.fasta" -db "$DB/PP-db/PP-db" -out "$OUT/MD_vs_PP.txt" \
 -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6


echo "////////// job done //////////"
