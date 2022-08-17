# Genome fractionation biais


## Description
The goal of this project is to measure the fractionation biais between homologous chromosomes of Malus domestica. It uses Prunus persica genome as a reference. 

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
- la liste des gènes de Malus domestica et Prunus persica, sous forme d'un fichier .lst par chromosome

## Usage
Dans un terminal, se placer dans le dossier genome-fractionation-biais. Lancer la commande suivante : 
```bash
./main.sh
```
Pour réaliser un test le premier arg doit être égal à 1 :
```bash
./main.sh 1
```
Pour forcer la réexécution des jobs déjà faits, le 2e arg doit être égal à 1 :
```bash
./main.sh 0 1
```

## Author
Adèle DESMAZIERES












