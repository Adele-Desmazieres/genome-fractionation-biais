# Genome Fractionation Biais


## Description
L'objectif de ce projet est de mesurer le biais de fractionnement entre les choromosomes homologues de Malus domestica. Le génome de Prunus persica sert de référence pour la comparaison interchromosomique de Malus domestica. 

The goal of this project is to measure the fractionation biais between homologous chromosomes of Malus domestica. Prunus persica genome is the reference for the interchomosomic comparison of Malus domestica. 

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
Télécharger le projet [Genome Fractionation Biais - IRHS](https://forgemia.inra.fr/irhs-bioinfo/genome-fractionation-biais) dans votre espace de stockage local. 

Copier ce projet de votre espace de stockage local sur le cluster de calcul. 
```bash
scp -r [path/to/]genome-fractionation-biais [user]@cluster-irhs.angers-nantes.inrae.fr:~
```

Création d'un environnement i-ADHoRe 3 grâce à conda. 
```bash 
conda create -c thiesgehrmann -c bioconda -c conda-forge -n iadhore_env iadhore
``` 

## Input files
Dans le dossier **data** :
- le protéome de Malus domestica au format .fasta avec le nom `Malus-domestica-proteome.fasta`
- le protéome de Prunus persica au format .fasta avec le nom `Prunus-persica-proteome.fasta`
- un dossier **MD_lst** contenant :
  - la liste des gènes de Malus domestica, sous la forme `Md01.lst`, `Md02.lst`, ... `Md17.lst`
- un dossie **PP_lst** contenant :
  - la liste des gènes de Prunus persica, sous la forme `Pp01.lst`, `Pp02.lst`, ... `Pp08.lst`

## Utilisation
### Lancement
1. Dans un terminal, se connecter au serveur de calcul par ssh. 
```bash
ssh cluster-irhs.angers-nantes.inrae.fr
```

2. Déplacer le répertoire courant dans le dossier genome-fractionation-biais. 
```bash
cd [path/to/]genome-fractionation-biais
```

3. Activer l'environnement virtuel avec i-ADHoRe. 
```bash
conda activate iadhore_env
```

4. Créer l'arborescence de dossiers. 
```bash
mkdir "data" "log"
```

5. Enregistrer les [input files](#input-files) dans le dossier **data**. 

6. Lancer le programme. 
```bash
./main.sh
```

### Résultats
Les résultats finaux sont dans le dossier **results/python/**. 

Les résultats intermédiaires sont dans **database/**, **results/blast/** et **results/iadhore/**. 

Le premier lancement prendra du temps, notamment pendant le blast. Lors des lancements suivants, le programme ne réalisera pas les tâches qu'il à déjà faites, sauf si les fichiers de résultats sont supprimés, ou si une option est passée pour forcer leur rééxecution. 

### Options
Afficher les options disponibles :
```bash
./main.sh -h
```

## Author
Adèle DESMAZIERES












