#!/bin/sh
# Before running, ensure it is correctly formatted by using dos2unix script_name.sh
# Make it executable by chmod 755 script_name.sh

# Name: run_association.sh
# commandline: qsub /run_association.sh <phenotypic file name> <phenotypic variable to test>

# ecrit les erreurs dans le fichier de sortie standard 
#$ -j y 

# shell que l'on veut utiliser 
#$ -S /bin/bash 

# indiquer son email pour suivre l'execution : 
#$ -M laurence.albar@ird.fr

# obtenir un message au demarrage (b) , a la fin (e), en cas d'abandon (a) 
#$ -m bea 

# la queue que l'on veut utiliser : 
#$ -q bioinfo.q

#$ -N association

# Export des variables d'environnement : 
#$ -V


###### Definition des variables de chemin

path_to_dir="/data2/projects/irigin/association";
path_to_tmp="/scratch/rice_association-$JOB_ID-$SGE_TASK_ID" 

###### Creation du repertoire temporaire sur noeud

mkdir $path_to_tmp
mkdir $path_to_tmp/$2

scp -rp nas:/$path_to_dir/genotype_data/* $path_to_tmp/
scp -rp nas:/$path_to_dir/phenotype_data/* $path_to_tmp/
scp -rp nas:/$path_to_dir/scripts/* $path_to_tmp/
echo "tranfert donnees master -> noeud";
cd $path_to_tmp

###### Execution du programme
Rscript $path_to_tmp/association_analysis_cluster.R $1 $2

##### Transfert des données du noeud vers master apres compression
scp -rp $path_to_tmp/$2 nas:$path_to_dir

echo "Transfert donnees node -> master";

#### Suppression du repertoire tmp noeud

rm -rf $path_to_tmp

echo "Suppression des donnees sur le noeud";
