# Genome fractionation biais


## Description
The goal of this project is to measure the fractionation biais between homologous chromosomes of Malus domestica. It uses Prunus persica genome as reference. 

## Prérequis
- Python 3.9.2
  - numpy 1.23.1
  - pandas 1.4.3
  - plotly 5.9.0
  - plotly.express 0.4.1
  - scipy 1.8.1
- i-ADHoRe 3.0.01
- blastp 2.11.0

## Installation

## Input files
Dans le dossier data :
- le protéome de Malus domestica au format .fasta
- le protéome de Prunus persica au format .fasta
- les listes les gènes de Malus domestica et Prunus persica, sous forme d'un fichier .lst par chromosome


## Usage
- créer l'arborescence de dossiers



- produire la base de données des gènes de PP avec makedb

        makeblastdb -in PP.fasta -dbtype prot -out PP-db

- soumettre à slurm le job blastp du protéome Malus domestica contre la base de données Prunus persica

        blastp -db databse/PP-db -query data/Malus-domestica-proteome.fasta -max_target_seqs=5 -evalue=1 -outfmt 6 -out results/blast/malus_vs_prunus.txt

- (concatener malus_vs_malus.txt + prunus_vs_prunus.txt + prunus_vs_malus2.txt dans blast_all.txt)

        cat results/iadhore/all/malus_vs_malus.txt > results/iadhore/all/blast_all.txt
        cat results/iadhore/all/prunus_vs_prunus.txt >> results/iadhore/all/blast_all.txt
        cat results/iadhore/all/prunus_vs_malus2.txt >> results/iadhore/all/blast_all.txt

- soumettre à slurm le job iadhore_job.sh

        sbatch submit/iadhore_job.sh

- supprimer les tabulations en double dans le fichier multiplicon_pairs.txt > multiplicon_pairs_modified.txt
  
        $ sed 's:\t\t*:\t:g' multiplicon_pairs.txt > multiplicon_pairs_modified.txt 

- lancer le script Python


## Author
Adèle DESMAZIERES












