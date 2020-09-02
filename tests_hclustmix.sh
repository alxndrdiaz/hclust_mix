#!/bin/bash

# runs hclust_mix on 29 data sets in de_souto_datasets/ folder, each matrix has two headers, the first header contains unique sample labels, the second header contains a type label (many samples can have the sampe type label) related to the condition, cell line, etc.  


input_dir='de_souto_datasets/'
input_files=`ls $input_dir*.tsv`
output_dir='all_results_'$input_dir

mkdir $output_dir; echo
for matfile in $input_files
do 
echo 'processing' $matfile
results_dir=$(basename $matfile .tsv)
mkdir $results_dir
echo; python hclust_mix.py $matfile -n -f -p ;  echo
mv attractor_results/ $results_dir
mv *.{ats,png,txt} $results_dir
mv $results_dir $output_dir
done
echo; echo results saved to $output_dir; echo  

