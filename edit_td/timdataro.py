# Python 2.7.6
# libraries: numpy 1.8.2; pandas 0.17.1
# timdataro : time-course data reading and output, loads time-course gene expression data, then selects genes by a minimum 
# number of  counts (mincounts=1e3 by default)  and finally generates an output file that contains these subset of genes.  

import os
import glob
import time
import numpy as np
import pandas as pd

# detects current working directory 
currentdirectory = os.getcwd()

# reads all .tmp intermediate files in currentdirectory
alltmp = glob.glob(currentdirectory + "/*.tmp")

# general information
filenumber = len(alltmp)
print "Number of files to process:", filenumber
print

# minimum number of counts to keep a gene
mincounts = 1e3 

# processing .tmp files
for filename in alltmp:
    loaded_genes =  pd.read_table(filename, index_col=0, header=None)
    print
    print "Processing:", filename
    gene_tags = loaded_genes.index.tolist()

    # selects numerical values only
    gene_counts = loaded_genes.values
    print "Initial gene number = ", len(gene_tags)

    # selects genes with mincounts or more
    rowsum = np.sum(gene_counts, axis=1) 
    newset = gene_counts[rowsum >= mincounts, :]
    newgenenumber = newset.shape[0] 
    print "Removing with minimum number of counts below = ", mincounts
    print "New gene number = ", newgenenumber

    # generates the new gene tags only for the genes that satified the condition in step 3.
    rownumber = rowsum.shape[0]
    interval = np.arange(rownumber) 
    print "Iterarion interval is from 1 to", interval.shape[0] 
    print
    newgenetags = [] 
    for j in interval: 
        if rowsum[j] >= mincounts: 
           newgenetags.append(gene_tags[j])
    print "Is new gene number equal to new number of gene tags?"
    if len(newgenetags) == newgenenumber:
       print "Yes."
    else:
       print "No. Remove all input files from the directory, copy your files and try again."

    # generates the output file
    header1 = []
    header2 = [] 
    stagenumber = np.arange(loaded_genes.shape[1]) 
    for k in stagenumber:
        header1.append(k) 
        header2.append('Sn')
        headers = [header1, header2]
        headers = pd.MultiIndex.from_arrays(headers, names=None)
    new_geneset = pd.DataFrame(newset, index=newgenetags,columns=headers)
    new_geneset.to_csv(filename, sep="\t")
    print "Process finished, updated file:", filename
    print
print "timdataro total run time in seconds = ", time.clock()
print
