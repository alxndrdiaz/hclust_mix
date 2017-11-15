#~language: Python 2.7.6
#~libraries:
#numpy 1.8.2
#pandas 0.17.1
#~~~~> author: Alexander Ramos Diaz
#----> description: timdataro = "time-course data reading and output"
# This code loads time-course gene expression data, then selects genes
# by their global frequency and finally generates an output file that contains
# the new selection of genes. Input and output formats can be changed.
# Important: the files to process will be updated, so keep in mind you need to
#work with these files and always have a backup of your original data.


import os
import glob
import time
import numpy as np
import pandas as pd

#detects current working directory and assigns it to the currentdirectory string
currentdirectory = os.getcwd()

#reads all tmp intermediate files in currentdirectory and assigns them to the list alltmp
alltmp = glob.glob(currentdirectory + "/*.tmp")

#General information
filenumber = len(alltmp)
print "Number of files to process:", filenumber
print

#Data processing in eight steps
for filename in alltmp:
    #~step1: load genes from each file as a dataframe
    loaded_genes =  pd.read_table(filename, index_col=0, header=None)
    print "----------------------------------------------------------------------------------------------------------------------------------------"
    print "Processing:", filename
    #~step2: selects original gene tags only
    gene_tags = loaded_genes.index.tolist()

    #~step3: selects numerical values only
    gene_counts = loaded_genes.values
    print "Initial gene number = ", len(gene_tags)

    #~step4: selects genes(detected transcripts) according to global frequency, 1e3 in this case
    rowsum = np.sum(gene_counts, axis=1) #axis=1 to sum rows
    frequencynumber = 1e3 #= minimum total number of times the gene was detected (condition)
    newset = gene_counts[rowsum >= frequencynumber, :] #select only the rows that satisfy the condition
    newgenenumber = newset.shape[0] #shows which is the new gene number
    print "Removing genes below a total detection frequency = ", frequencynumber
    print "--> New gene number = ", newgenenumber

    #~step5: generates the new gene tags only for the genes that satified the condition in step 3.
    rownumber = rowsum.shape[0] #the first component of "shape"  is the row number
    interval = np.arange(rownumber) #generates the iteration interval [0,1,...n]
    print "Iterarion interval is from 1 to", interval.shape[0] #verifies that it will iterate over all the rows.
    print
    newgenetags = [] #new empty list
    for j in interval: # fills the list 'newgenetags' only when the condition holds
        if rowsum[j] >= frequencynumber: #1e5 must be the same value in step4
           newgenetags.append(gene_tags[j])
    print "Is new gene number equal to new number of gene tags?"
    if len(newgenetags) == newgenenumber:
       print "--> Yes."#verifies that number_of_newsgenetags = newgenenumber.
    else:
       print "No. Remove all input files from the directory, copy your files and try again."

    #~step6: generates a header for the output file
    header1 = [] #generates a number list that will be used as first header
    header2 = [] #generates a letter list that will be used as extra header for the tmp file
    stagenumber = np.arange(loaded_genes.shape[1]) # stagenumer = positions where the tag will be added
    for k in stagenumber:
        header1.append(k) #0, 1, etc..
        header2.append('Sn') #Sn = stage number,it is a generic tag
        headers = [header1, header2]
        headers = pd.MultiIndex.from_arrays(headers, names=None)

    #~step7: generates a dataframe to generate the output file
    new_geneset = pd.DataFrame(newset, index=newgenetags,columns=headers)

    #~step8: generates the output file
    new_geneset.to_csv(filename, sep="\t")
    print "Process finished, updated file:", filename
    print "----------------------------------------------------------------------------------------------------------------------------------------"
print "timdataro total run time in seconds = ", time.clock()
print
