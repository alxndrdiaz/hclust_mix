#!/bin/bash

# removes lines that are not required for processing such as metadata lines and empty lines, also removes lines associated to ERCC transcripts
# this script works specifically for the time-course expression matrices data set from Levin et al. Nature (2016), availabe at NCBI GEO (GSE70185):
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70185
 
# removes unnecessary lines in all .tab input files:
echo
echo "Editing the following files:"
for filename in *.{txt,tab}
do
  # prints file list
  echo $filename  
  # generates a new file name using the NCBI and species code only 
  newfilename=$(echo $filename | cut -c1-11) 
  # reads lines starting at position 2 and generates a new file with the data (excludes first line)
  tail -n +2 $filename > $newfilename.aux
  # in auxiliary file removes all lines that contain ERCC identifiers (spike-in RNAs)
  sed -i '/ERCC/d' $newfilename.aux   
  # removes last five lines containing data set information	
  head -n -5 $newfilename.aux > $newfilename.tmp
done
echo
# removes .aux files
rm *.aux
  
# runs timdataro.py
python timdataro.py

# removes any remaining empty lines:  
for updatedname in *.tmp
do
  # obtains NCBI + species code only
  finalname=$(echo $updatedname | cut -c1-11)
  # removes blank lines avoiding 
  grep -v '^[[:blank:]]*$' $updatedname > $finalname.tsv
done

# removes auxiliary .tmp files and keeps only input .tab files (raw input) and .tsv (files for hclust)
rm *.tmp

# generates general information of the output and separated folders for raw data and edited data: 
# moves all updated files to a folder
mkdir edited_data
echo "Updated files that can be used with hclust are:"
for finalfile in *.tsv
do
  echo $finalfile
  mv $finalfile edited_data
done
echo "Please open the edited_data folder to check updated time-course files."
# moves all input .tab files to a folder
mkdir raw_data
echo
for tabfile in *.{txt,tab}
do
 mv $tabfile raw_data
done  
echo "Raw time-course data moved to raw_data folder."
echo "File editing finished".
echo
