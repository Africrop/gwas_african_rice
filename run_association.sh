###################################################################################################################################
#
# Copyright 2018 IRD
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD and Grenoble-Alpes University
#
# Written by Philippe Cubry
#
###################################################################################################################################

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
