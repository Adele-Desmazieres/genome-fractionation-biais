#!/bin/bash

# TODO : ajouter une commande pour créer l'arborescence de dossiers

makeblastdb -in data/Prunus-persica-proteome.fasta -dbtype prot -out database/PP-db
makeblastdb -in data/Malus-domestica-proteome.fasta -dbtype prot -out database/MD-db

sbatch submit/blast_job.sh

# TODO : s'assurer que le job donné à slurm termine avant de lancer les commandes suivantes
#MD_vs_PP.txt results/blast/MD_vs_MD.txt results/blast/PP_vs_PP
cat results/blast/*.txt > results/blast/blast_all.txt

sbatch submit/iadhore_job.sh

# TODO : s'assurer que le job donné à slurm termine avant de lancer les commandes suivantes
sed 's:\t\t*:\t:g' multiplicon_pairs.txt > multiplicon_pairs_modified.txt 

cd script/python
python3 main2.py