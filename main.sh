#!/bin/bash

TEST=$1

# variables de chemins
DATA="data"
DB="database"
SUBMIT="submit"
RES="results"

# si arg1 = 1 alors tous les dossiers doivent finir par _test
if [[ TEST -eq 1 ]]
then
    DATA="${DATA}_test"
    DB="${DB}_test"
    RES="${RES}_test"
fi

# TODO : ajouter une commande pour créer l'arborescence de dossiers si nécessaire

# crée les bases de données nécessaires au blast
printf "\n### MAKEBLASTDB ###"
makeblastdb -in $DATA/Prunus-persica-proteome.fasta -dbtype prot -out $DB/PP-db
makeblastdb -in $DATA/Malus-domestica-proteome.fasta -dbtype prot -out $DB/MD-db

# soumet le blast à slurm
printf "\n### BLAST ###\n"
sbatch --wait --export=ALL,DATA=$DATA,DB=$DB,OUT="$RES/blast" $SUBMIT/blast_job.sh 
# --wait permet d'exit seulement quand le job termine, fait attendre le script en attendant
# --export permet d'exporter des variables

# regroupe les 3 blast PP-PP MD-MD et MD-PP dans un même fichier
cat $RES/blast/*.txt > $RES/blast/all_vs_all.txt

# déplace les fichiers d'entrée de DATA dans data_tmp pour que iadhore les choppe là-bas sans avoir à modifier son path
# TODO : s'assurer que tmp existe ?
cp $RES/blast/* tmp/blast/
cp $DATA/MD_lst/* tmp/data/MD_lst
cp $DATA/PP_lst/* tmp/data/PP_lst

# soumet le travail iadhore à slurm
printf "\n### IADHORE ###\n"
sbatch --wait $SUBMIT/iadhore_job.sh

# déplace les résultats de tmp vers le dossier de résultats
cp tmp/iadhore/* $RES/iadhore/

# coupe la tabulation en trop dans le fichier multiplicon_pairs.txt
sed 's:\t\t*:\t:g' $RES/iadhore/multiplicon_pairs.txt > $RES/iadhore/multiplicon_pairs_modified.txt 

# lance le script python
printf "\n### PYTHON ###\n"
cd scripts/python
python3 main2.py > ../../$RES/python/fractionation_stat.txt

printf "\ndone\n"