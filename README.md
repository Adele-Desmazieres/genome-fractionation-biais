# Genome fractionation biais


## Description
The goal of this project is to measure the fractionation biais between homologous chromosomes of Malus domesticae.

## Prérequis
- Python v3.9.2
- modules Python ?
- iadhore v ?
- blast v ?

## Installation

## Input files
Dans le dossier data :
- le proteome de Malus domestica .fasta
- le proteome de Prunus persica .fasta
- les listes les gènes de Malus domestica et Prunus persica, un fichier .lst par chromosome


## Usage
- reproduire l'arborescence de dossiers
- faire la database de PP avec makedb
- soumettre les jobs blastp MD-MD PP-PP et PP-db2
- concatener malus_vs_malus.txt + prunus_vs_prunus.txt + prunus_vs_malus2.txt dans blast_all.txt
- soumettre le job iadhore_job.sh
- supprimer la tabulation en trop dans le fichier multiplicon_pairs.txt > multiplicon_pairs_modified.txt
- lancer le script Python


## Author
Adèle DESMAZIERES












