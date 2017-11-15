~~~~~~~~~General information about timdataro.py~~~~~~~~~

description: timdataro is a short name for "time-course data reading and output"

~language: Python 2.7.6
~libraries:
 numpy 1.8.2
 pandas 0.17.1

This code loads time-course gene expression data, then selects genes
by their global frequency and finally generates an output file that contains
the new selection of genes. Input and output formats can be changed.
Important: the files to process will be updated, so keep in mind you need to
work with these files and always have a backup of your original data.


Features:

0. edit_td.sh calls timdataro in order to edit time-course gene expression data files. However it can work separately fron Linux terminal:
   
   python timdataro.py 

1. Works if the files to process are in the same folder that the .py file, otherwise the files' directory should be specified.

2. Uses a general strategy: reads current working directory -> reads files in the specified format (default is .tmp) -> process files
   using pandas and numpy -> generates a .tmp updated file with the required format and conditions.

3. pandas is used to read the data as a data frame which simplifies indexing and numpy is used for matrix and array computations.

4. Global frequency means the total number of times a gene was detected considering all the stages, for example if an experiment that 
   considered 3 stages  X, Y, Z if geneA was detected 1, 2 and 3 times at stages X, Y, Z respectively, the global frequency would be
   1 + 2 + 3 = 6. The defaul global frequency used is 1000 (or 1e3), this means timdataro will consider genes that have a global 
   frequency above or equal to 1000 and all other genes will be removed, this updated gene set is available in the .tsv file.

5. Input and output files can be changed depending on how you need to make your data analysis.

author: Alexander Ramos Diaz
e-mail: rene.ramos@cinvestav.mx
   
