#!/bin/bash

TEST=$1

DATA="data"
DB="database"
SUBMIT="submit"
RES="results"

if [[ TEST -eq 1 ]]
then
    DATA="${DATA}_test"
    DB="${DB}_test"
    SUBMIT="${SUBMIT}_test"
    RES="${RES}_test"
fi

# TODO : ajouter une commande pour créer l'arborescence de dossiers si nécessaire

makeblastdb -in $DATA/Prunus-persica-proteome.fasta -dbtype prot -out $DB/PP-db
makeblastdb -in $DATA/Malus-domestica-proteome.fasta -dbtype prot -out $DB/MD-db

sbatch --export=ALL,DATA=$DATA,DB=$DB,OUT="$RES/blast" $SUBMIT/blast_job.sh 

# TODO : s'assurer que le job donné à slurm termine avant de lancer les commandes suivantes
#MD_vs_PP.txt results/blast/MD_vs_MD.txt results/blast/PP_vs_PP
cat $RES/blast/*.txt > $RES/blast/all_vs_all.txt

# TODO : déplacer les fichiers d'entrée de DATA dans data_tmp pour que iadhore les choppe là-bas
# TODO s'assurer que tmp existe ?
cp $RES/blast/* tmp/blast/
sbatch $SUBMIT/iadhore_job.sh
cp tmp/iadhore/* $RES/iadhore/

# TODO : s'assurer que le job donné à slurm termine avant de lancer les commandes suivantes
sed 's:\t\t*:\t:g' multiplicon_pairs.txt > multiplicon_pairs_modified.txt 

cd scripts/python
python3 main2.py