#!/bin/bash
date # affiche la date et l'heure actuels

# affiche la documentation
print_usage() {
	printf "Options: \
	\n\t-t : test \
	\n\t-f : forced \
	\n"
	exit
}

# récupère les options du script, source : https://google.github.io/styleguide/shellguide.html 
test_flag=0 # si =1 lance les commandes sur les dossiers_test
forced_flag=0 # si =1 force le lancement des commandes même si les résultats ont été calculés auparavant, écrase les résultats précédents
while getopts 'tf' flag; do
  case "${flag}" in
	t) test_flag=1 ;;
	f) forced_flag=1 ;;
	*) print_usage
	   exit 1 ;;
  esac
done
readonly test_flag
readonly forced_flag

# variables de chemins
DATA="data"
DB="database"
RES="results"
SUBMIT="submit"
TMP="tmp"

# traitement des options
# si arg1 = 1 alors tous les dossiers doivent finir par _test sauf submit
if [[ test_flag -eq 1 ]]
then
	printf "OPTION TEST"
	DATA="${DATA}_test"
	DB="${DB}_test"
	RES="${RES}_test"
fi

# supprime les résultats précédents si l'option forced vaut 1
if [[ forced_flag -eq 1 ]]
then
	printf "OPTION FORCED"
	rm -r $DB $RES
fi


# crée l'arborescence de dossiers si elle n'existe pas
mkdir -p $DB $RES/{blast,iadhore,python}


# crée les bases de données de blast de PP et MD si elles n'existent pas (si leur dossier est inexistant ou vide)
# base de données de PP
if [[ ! -d "$DB/PP-db/" || ! "$(ls -A $DB/PP-db/)" ]]
then
	printf "\n### MAKEBLASTDB PP ###"
	makeblastdb -in "$DATA/Prunus-persica-proteome.fasta" -dbtype prot -out "$DB/PP-db/PP-db"
	printf "makeblastdb PP: done\n"
fi
# base de données de MD
if [[ ! -d "$DB/MD-db/" || ! "$(ls -A $DB/MD-db/)" ]]
then
	printf "\n### MAKEBLASTDB MD ###"
	makeblastdb -in "$DATA/Malus-domestica-proteome.fasta" -dbtype prot -out "$DB/MD-db/MD-db"
	printf "makeblastdb MD: done\n"
fi


# soumet le blast à slurm si son résultat n'existe pas déjà
if [[ ! "$(ls -A $RES/blast)" ]]
then
	printf "\n### BLAST ###\n"
	printf "DATA=$DATA\nDB=$DB\nOUT=$RES/blast\n"
	sbatch --wait --export=ALL,DATA=$DATA,DB=$DB,OUT="$RES/blast" $SUBMIT/blast_job.sh 
	# --wait permet d'exit seulement quand le job termine, fait attendre le script sinon
	# --export permet d'exporter des variables vers le script du job

	# regroupe les 3 blast PP-PP MD-MD et MD-PP dans un même fichier
	cat $RES/blast/*.txt > $RES/blast/all_vs_all0.txt
	# traite le fichier pour ne conserver que les 2 premières colonnes, contenant les paires de gènes
	cut -f1,2 $RES/blast/all_vs_all0.txt > $RES/blast/all_vs_all.txt

	printf "blast: done\n"
fi


# soumet le iadhore job à slurm si son résultat n'existe pas déjà
if [[ ! "$(ls -A $RES/iadhore)" ]]
then
	printf "\n### IADHORE ###\n"
	# déplace les fichiers d'entrée de DATA dans TMP/data pour que iadhore les choppe là-bas sans avoir à modifier son input path
	if [[ -d $TMP ]]
	then
		rm -r $TMP
	fi

	mkdir -p $TMP/iadhore/

	ln -s ../$RES/blast/ $TMP/ # project/tmp/blast/files
	ln -s ../$DATA/ $TMP/ # project/tmp/data/files
	
	if [[ $DATA != "data" ]]
	then
		mv $TMP/$DATA $TMP/data # renomme le dossier
	fi
	#rm -r $TMP/$DATA/
	#cp -r ../$DATA/PP_lst $TMP/data/

	# soumet le job
	sbatch --wait $SUBMIT/iadhore_job.sh

	# déplace les résultats de tmp vers le dossier de résultats
	mv $TMP/iadhore/* $RES/iadhore/

	# si le fichier d'output de iadhore existe alors
	# coupe la tabulation en trop dans le fichier multiplicon_pairs.txt si multiplicon_pairs.txt existe
	if [ -f $RES/iadhore/multiplicon_pairs.txt ]
	then
		sed 's:\t\t*:\t:g' $RES/iadhore/multiplicon_pairs.txt > $RES/iadhore/multiplicon_pairs_modified.txt
	fi

	# supprime le dossier temporaire
	if [[ -d $TMP ]]
	then
		rm -r $TMP
	fi

	printf "iadhore: done\n"
fi


# lance le script python
printf "\n### PYTHON ###\n"
python3 scripts/python/main2.py $test_flag > $RES/python/fractionation_stat.txt 
printf "python: done\n\n"

date # affiche la date et l'heure actuels
printf "script: done\n"