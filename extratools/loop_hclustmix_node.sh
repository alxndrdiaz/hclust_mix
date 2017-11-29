#!/bin/bash
#PBS -N hopfield_analysis
#PBS -l nodes=1:ppn=10,vmem=64gb,walltime=05:00:00
#PBS -q default
#PBS -V
#PBS -o hopfield.out
#PBS -e hopfield.err

module load Python/2.7.6

	
for inputfile in *.tsv 
do
 #obtains NCBI code + species code only:
 ncbiandspecies=$(echo $inputfile | cut -c1-11)
 specieslabel=$(echo $inputfile | cut -c10-11)
 python hclust_mix_node.py $inputfile -n 
 #moves the results to a new directory:
 foldername=$ncbiandspecies	
 mkdir $foldername
 mv $inputfile $foldername
 mv *.{png,st,trk} $foldername
 cd $foldername
 states='_states'
 track='_track'
 mkdir $specieslabel$states
 mkdir $specieslabel$track
 mv *.st $specieslabel$states
 mv *.trk $specieslabel$track
 cd ..
done 



