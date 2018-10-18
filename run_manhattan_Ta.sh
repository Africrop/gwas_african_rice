#!/bin/sh
# Before running, ensure it is correctly formatted by using dos2unix script_name.sh
# Make it executable by chmod 755 script_name.sh

# Name: run_association.sh
# commandline: qsub /run_manhattan_Ta.sh <variable> <FDR>

# ecrit les erreurs dans le fichier de sortie standard 
#$ -j y 

# shell que l'on veut utiliser 
#$ -S /bin/bash 

# indiquer son email pour suivre l'execution : 
#$ -M philippe.cubry@ird.fr

# obtenir un message au demarrage (b) , a la fin (e), en cas d'abandon (a) 
#$ -m bea 

# la queue que l'on veut utiliser : 
#$ -q bioinfo.q

#$ -N manhattanplot

# Export des variables d'environnement : 
#$ -V

cd /data2/projects/irigin/association/

###### Execution du programme
Rscript ./scripts/manhattan_plots_Ta.R $1 $2
