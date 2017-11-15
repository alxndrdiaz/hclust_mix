~~~~~~~~~General information about edit_td.sh~~~~~~~~~

edit_td is a short name for "edit time-course data"

Description:

This small code edits the original gene expression time-course data (for our project the data are time-course gene expression for  ~100  developmental stages of ten metazoan species) available in tables with .tab format, delimiter(separator) in these files is tab or \t which both python and bash recognize. In these files there is a first line [position 1] that tags the developmental stages, in some files there are also control sequences tagged with ERCC, and at the end of all files there are 5 lines summarizing the information of the set; all of these tags and extra information are not required to run hclust or other algorithms, and should be removed. Output data will be .tsv (tab-separated-values) files that hclust requires for processing and generating its own results. 

Features:

0. Runs in bash (Ubuntu Linux and other Linux distributions).

1. Original time-course gene expression data avilable at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70185
   edit_td works for removing tags that are unnecessary for hclust in this specific data set. If you use this code for other 
   time-course gene expression data, it is recommendable to check the format of the text files before processing 
   your data and modify your code accordingly.

3. The script calls a python code called timdataro.py, this code reads each file and selects genes that have a frequency (ocurrence
   in all the stages) greater than or equal to a specific value, the default value used is 1000 (1e3), this means that only the subset
   of all the genes that were detected 1000 times or more will be considered for further data analysis. For more information check
   the timdataro_readme.txt file.


2. Once edit_td finishes, you will have a set of .tsv files that are required to process the time-course data in hclust.


How to run the script?

The easiest way to run edit_td.sh is to have all the files to edit and both edit_td.sh and timdataro.py in the same directory. In this example it is assumed that you already downloaded both files and are in the Desktop directory.

0. Open the Linux terminal and write the following command (then press Enter):

   bash

1. Go to a directory (Desktop) and create a folder (time_course or other name) where you want to edit your 
   time-course data files, write the following commands and press Enter after each one:
  
   cd Desktop
   mkdir time_course
  
2. Now move all the files to edit to the time_course folder, in this example the .tab files are in Desktop and you move them to time_course:
   
   mv *.tab time_course

3. Also move edit_td.sh and timdataro.py to the same folder:

   mv edit_td.sh time_course
   mv timdataro.py time_course

4. Go to time_course directory:
   
   cd time_course

5. Once in time_course folder, type the following commands (press Enter after each one):
   
   chmod +x edit_td.sh
   chmod +x timdataro.py

6. Finally run the script, just write the following command (and then press Enter):
   
   ./edit_td.sh

Results:
In time_course folder you must see a folder called raw_data where all your original .tab files were moved, also there is another folder
called edited_data which contains all the updated .tsv files that hclust requires.


author: Alexander Ramos Diaz.
e-mail: rene.ramos@cinvestav.mx



