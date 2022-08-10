#!/bin/bash

# TODO : ajouter une commande pour créer l'arborescence de dossiers

makeblastdb -in data/Prunus-persica-proteome.fasta -dbtype prot -out database/PP-db
makeblastdb -in data/Malus-domestica-proteome.fasta -dbtype prot -out database/MD-db

# TODO : modifier le script du blast pour que le path soit valide et pour qu'il réalise tous les blasts MD->PP MD->MD et PP->PP
sbatch submit/blast_job.sh

# TODO : s'assurer que le job donné à slurm termine avant de lancer les commandes suivantes
cat results/iadhore/all/malus_vs_prunus.txt results/iadhore/all/malus_vs_malus.txt results/iadhore/all/prunus_vs_prunus.txt > results/iadhore/all/blast_all.txt

sbatch submit/iadhore_job.sh

# TODO : s'assurer que le job donné à slurm termine avant de lancer les commandes suivantes
sed 's:\t\t*:\t:g' multiplicon_pairs.txt > multiplicon_pairs_modified.txt 

cd script/python
python3 main2.py