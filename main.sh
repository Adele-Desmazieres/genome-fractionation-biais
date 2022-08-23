#!/bin/bash

printf "script: running...\n"
date # affiche la date et l'heure actuels

# affiche comment obtenir de l'aide
print_usage() {
	printf "Saisissez \" $0 -h \" pour plus d'informations.\n"
}

# affiche la documentation de la commande
print_options() {
	printf "\nUSAGE\n\t$0 [OPTIONS]\n \
	\nOPTIONS \
	\n\t-a : force all reexecution, same as -bip \
	\n\t-b : force blast et blastdb reexecution \
	\n\t-i : force i-ADHoRe reexecution \
	\n\t-p : force Python reexecution \
	\n\t-t : test, execute in directories_test\
	\n\t-h : display options \
	\n\n"
	exit 0
}

# récupère les options du script, source : https://google.github.io/styleguide/shellguide.html 
test_flag=0 # si =1, lance les commandes sur les dossiers_test
blast_flag=0 # si =1, force l'execution de blast
iadhore_flag=0 # si =1, force l'execution de iadhore
python_flag=0 # si =1, force l'execution de python

while getopts 'tbipah' flag; do
  case "${flag}" in
	t) test_flag=1 ;;
	b) blast_flag=1 ;;
	i) iadhore_flag=1 ;;
	p) python_flag=1 ;;
	a) blast_flag=1
	   iadhore_flag=1 
	   python_flag=1 ;;
	h) print_options ;;
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
SUBMIT="scripts/bash"
TMP="tmp"

# traitement des options
# si arg1 = 1 alors les dossiers de data, db et résultats doivent finir par _test
if [[ test_flag -eq 1 ]]
then
	printf "OPTION TEST\n"
	DATA="${DATA}_test"
	DB="${DB}_test"
	RES="${RES}_test"
fi


# crée l'arborescence de dossiers si elle n'existe pas
mkdir -p $DB $RES/{blast,iadhore,python}

# crée les bases de données de blast de PP et MD si leur dossier n'existe pas, ou est vide, ou le flag blast est levé

# base de données MD
if [[ blast_flag -eq 1 || ! -d "$DB/MD-db/" || ! "$(ls -A $DB/MD-db/)" ]]
then
	printf "\n### MAKEBLASTDB MD ###"

	if [[ -d "$DB/MD-db" && "$(ls -A $DB/MD-db)" ]] ; then rm -r $DB/MD-db/* ; fi # supprime l'ancienne db si elle existe
	makeblastdb -in "$DATA/Malus-domestica-proteome.fasta" -dbtype prot -out "$DB/MD-db/MD-db"
	printf "makeblastdb MD: done\n"
fi
# base de données PP
if [[ blast_flag -eq 1 || ! -d "$DB/PP-db/" || ! "$(ls -A $DB/PP-db/)" ]]
then
	printf "\n### MAKEBLASTDB PP ###"

	if [[ -d "$DB/PP-db" && "$(ls -A $DB/PP-db)" ]] ; then rm -r $DB/PP-db/* ; fi # supprime l'ancienne db si elle existe
	makeblastdb -in "$DATA/Prunus-persica-proteome.fasta" -dbtype prot -out "$DB/PP-db/PP-db"
	printf "makeblastdb PP: done\n"
fi


# soumet le blast à slurm si son résultat n'existe pas déjà ou le flag blast est levé
if [[ blast_flag -eq 1 || ! "$(ls -A $RES/blast)" ]]
then
	printf "\n### BLAST ###\n"
	if [[ "$(ls -A $RES/blast)" ]] ; then rm -r $RES/blast/* ; fi
	printf "DATA=$DATA\nDB=$DB\nOUT=$RES/blast\n"

	# soumet le job blast
	sbatch --wait --export=ALL,DATA=$DATA,DB=$DB,OUT="$RES/blast" $SUBMIT/blast_job.sh 
	#    --wait permet d'exit seulement quand le job termine, fait attendre le script sinon
	#    --export permet d'exporter des variables vers le script du job

	# regroupe les 3 blast PP-PP MD-MD et MD-PP dans un même fichier
	cat $RES/blast/*.txt > $RES/blast/all_vs_all0.txt
	# traite le fichier pour ne conserver que les 2 premières colonnes, contenant les paires de gènes
	cut -f1,2 $RES/blast/all_vs_all0.txt > $RES/blast/all_vs_all.txt

	printf "blast: done\n"
fi



# remplit le fichier de config iadhore avec les chromosomes présents dans data
# liste des fichiers de data/MD_lst de la forme :
# filename1 this/is/the/path/filename1.extension
# filename2 this/is/the/path/filename2.extension
# filename3 this/is/the/path/filename3.extension
make_iadhore_config() {
	file="scripts/iadhore/iadhore.ini"
	cat "scripts/iadhore/iadhore1_model.ini" > $file # copie le modèle vide dans le fichier de config
	
	printf "\ngenome = Malus_domestica\n" >> $file 
	ls -1 $TMP/data/MD_lst/* | sed 's/.*/& &/' | sed -e "s/$TMP\/data\/MD_lst\///" -e 's/.lst / /' >> $file
	# lister deux fois les sorties séparées par un espace, puis supprimer le chemin et l'extension de la première sortie

	printf "\ngenome = Prunus_persica\n" >> $file
	ls -1 $TMP/data/PP_lst/* | sed 's/.*/& &/' | sed -e "s/$TMP\/data\/PP_lst\///" -e 's/.lst / /' >> $file

	cat $file # affiche le fichier de config
	printf "\n"
	#ls -R tmp/data
	#printf "\n"
}


# soumet le iadhore job à slurm si son résultat n'existe pas déjà ou le flag iadhore est levé
if [[ iadhore_flag -eq 1 || ! "$(ls -A $RES/iadhore)" ]]
then
	printf "\n### IADHORE ###\n"

	if [[ "$(ls -A $RES/iadhore)" ]] ; then rm -r $RES/iadhore/* ; fi # supprime les précédants resultats iadhore

	# supprime le dossier tmp s'il existe puis en crée un vide
	if [[ -d $TMP ]] ; then rm -r $TMP ; fi
	mkdir -p $TMP/iadhore/

	# déplace les fichiers d'entrée de DATA dans TMP/data pour que iadhore les choppe là-bas sans avoir à modifier son input path
	ln -s ../$RES/blast/ $TMP/ # tmp/blast/files
	ln -s ../$DATA/ $TMP/ # tmp/data_?/files
	
	# renomme le dossier du lien symbolique, afin qu'il s'appelle forcément "data"s
	if [[ $DATA != "data" ]]
	then
		mv $TMP/$DATA $TMP/data # renomme le dossier tmp/data
	fi

	make_iadhore_config # crée le fichier de config iadhore et l'affiche
	
	# soumet le job iadhore
	sbatch --wait $SUBMIT/iadhore_job.sh

	# déplace les résultats de tmp vers le dossier de résultats
	mv $TMP/iadhore/* $RES/iadhore/

	# si le fichier d'output de iadhore existe alors
	# remplace la double tabulation par une simple tabulation, du fichier multiplicon_pairs.txt s'il existe
	if [ -f $RES/iadhore/multiplicon_pairs.txt ]
	then
		sed 's:\t\t*:\t:g' $RES/iadhore/multiplicon_pairs.txt > $RES/iadhore/multiplicon_pairs_modified.txt
	fi

	if [[ -d $TMP ]] ; then rm -r $TMP ; fi # supprime le dossier temporaire

	printf "iadhore: done\n"
fi


# lance le script python
if [[ python_flag -eq 1 || ! "$(ls -A $RES/python)" ]]
then 
	printf "\n### PYTHON ###\n"
	python3 scripts/python/main2.py $test_flag | tee $RES/python/fractionation_stat.txt 
	printf "python: done\n\n"
fi


date # affiche la date et l'heure actuels
printf "script: done\n"