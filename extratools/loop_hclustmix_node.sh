#!/bin/bash

#~author: Alexander Ramos Diaz.
#~description: this script takes all the .tsv input files and send them to hclust_mix one by one, then creates a folder for each file with 
#the original input file and the results.
# ~~~~~~You are FREE to use or adapt this code for your own project~~~~~~~~~


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

 
