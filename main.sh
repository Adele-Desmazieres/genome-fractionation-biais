#!/bin/bash

TEST=$1

# variables de chemins
DATA="data"
DB="database"
RES="results"
SUBMIT="submit"
# si arg1 = 1 alors tous les dossiers doivent finir par _test sauf submit
if [[ TEST -eq 1 ]]
then
    DATA="${DATA}_test"
    DB="${DB}_test"
    RES="${RES}_test"
fi

# crée l'arborescence de dossiers si elle n'existe pas
mkdir -p $DB $RES/{blast,iadhore,python}

# crée les bases de données de blast de PP et MD si elles n'existent pas (si leur dossier est inexistant ou vide)
# base de données de PP
if [[ ! -d "$DB/PP-db" || ! "$(ls -A $DB/PP-db)" ]]
then
    printf "\n### MAKEBLASTDB PP ###"
    makeblastdb -in "$DATA/Prunus-persica-proteome.fasta" -dbtype prot -out "$DB/PP-db"
fi
# base de données de MD
if [[ ! -d "$DB/MD-db" || ! "$(ls -A $DB/MD-db)" ]]
then
    printf "\n### MAKEBLASTDB MD ###"
    makeblastdb -in "$DATA/Prunus-persica-proteome.fasta" -dbtype prot -out "$DB/MD-db"
fi

# soumet le blast à slurm
if [[ ! "$(ls -A $RES/blast)" ]]
then
    printf "\n### BLAST ###\n"
    sbatch --wait --export=ALL,DATA=$DATA,DB=$DB,OUT="$RES/blast" $SUBMIT/blast_job.sh 
    # --wait permet d'exit seulement quand le job termine, fait attendre le script en attendant
    # --export permet d'exporter des variables
fi


# regroupe les 3 blast PP-PP MD-MD et MD-PP dans un même fichier
cat $RES/blast/*.txt > $RES/blast/all_vs_all.txt

# déplace les fichiers d'entrée de DATA dans data_tmp pour que iadhore les choppe là-bas sans avoir à modifier son path
# TODO : s'assurer que tmp existe et est vide
rm -r tmp/blast/* tmp/data/*
cp $RES/blast/* tmp/blast/
cp $DATA/MD_lst/* tmp/data/MD_lst
cp $DATA/PP_lst/* tmp/data/PP_lst

# soumet le iadhore job à slurm si dossier de résultat d'iadhore est vide
if [[ ! "$(ls -A $RES/iadhore)" ]]
then
    printf "\n### IADHORE ###\n"
    sbatch --wait $SUBMIT/iadhore_job.sh
fi

# déplace les résultats de tmp vers le dossier de résultats
cp tmp/iadhore/* $RES/iadhore/

# coupe la tabulation en trop dans le fichier multiplicon_pairs.txt
sed 's:\t\t*:\t:g' $RES/iadhore/multiplicon_pairs.txt > $RES/iadhore/multiplicon_pairs_modified.txt 

# lance le script python
printf "\n### PYTHON ###\n"
python3 scripts/python/main2.py > $RES/python/fractionation_stat.txt $TEST

printf "\n### SCRIPT DONE ###\n"