#~ language: Python 2.7.6
#~ libraries:
#~ numpy 1.8.2
#~ pandas 0.17.1
#~ scipy 0.13.3

#~ Attractor state analysis: this code takes the output attractors from hclust_mix
#~ in .ats format and analyze overrepresented and underrepresented genes, which
#~ are the predictions of the identified attractors. It also compares attractors
#~ in order to determine if they are equal or if they are different and how they
#~ are structured. Finally it provides output files for further analysis.

import os
import glob
import time
import numpy as np
import pandas as pd

#~~>step0: reading files

#~ detects current working directory and assigns it to the currentdirectory string
currentdirectory = os.getcwd()

#~ reads all ats intermediate files in currentdirectory and assigns them to the list alltmp
all_ats = glob.glob(currentdirectory + "/*.ats")

#~ number of attractors to analyze
filenumber = len(all_ats)

#~ generates empty lists
attlist = [] #~ storages all attractors in table format
attpositions = [] #~ storages stage positions of all attractos
undernumbers = [] #~ storages number of overrepresented genes in each attractor
overnumbers = [] #~ storages number of underrepresented genes in each attractor

#~~> step1: basic analysis of each attractor

#~ generates a list of attractor files by filling the empty list, then
#~ classifies genes according to overrepresentation/underrepresented (+1 or -1)
for at in all_ats:
    attdata = pd.read_table(at, index_col=0, header=0)
    attlist.append(attdata)
    attstage = attdata.columns[0] #~ obtains index name for further selection
    attpositions.append(attstage)
    #~ underrepresentation selection analysis:
    underset = attdata[attdata[attstage] == -1] #~ selects underrepresented genes
    under_gene_list =  underset.index.tolist() #~ generates a list only with gene names
    unumber = len(under_gene_list)
    undernumbers.append(unumber)
    undersetname =  attstage + "_attractor_" + "under.txt"
    #~ overrepresentation selection analysis:
    overset = attdata[attdata[attstage] == 1]  #~ selects overrepresented genes
    over_gene_list = overset.index.tolist()
    onumber = len(over_gene_list)
    overnumbers.append(onumber)
    oversetname =  attstage + "_attractor_" + "over.txt"

#~ using the attractor list to create a file containing all attractors
all_attractors = pd.concat(attlist, axis=1)
all_attractors.to_csv("all_atrractors.aat", sep="\t")

#~~> identifies repeated and complementary attractors and provides a list of the results
#~ it also provides an statistical summary and graphics.

vectors = []
repeated_attractors = []

for att in attlist:
    vector = att.values
    vectors.append(vector)

aux_vectors = vectors
aux_positions = attpositions
pairs_repeated = []

iterate = np.arange(len(vectors))

for i in iterate:
    for j in iterate:
        if i != j:
           if np.array_equal(vectors[i], aux_vectors[j]) == True:
              repeated = "(" + aux_positions[i] + "," + aux_positions[j] + ")"
              pairs_repeated.append(repeated)


itpairs = np.arange(len(pairs_repeated))
updated_pairs = []
for i in itpairs: #~ prints only the pair of attractors by selecting only odd positions
    j = i + 1
    if j % 2 != 0:
       updated_pairs.append(pairs_repeated[j])


genenumbers = np.vstack((undernumbers, overnumbers))
genenumbers = np.transpose(genenumbers)
numberlabels = ['u(-)genes','o(+)genes']

attindex = []
for position in attpositions:
    indexname = 'at stage ' + position
    attindex.append(indexname)
print
att_info = pd.DataFrame(genenumbers, index=attindex, columns=numberlabels)
att_info.to_csv("attractor_information.txt", sep="\t")

print
print "Attractors: "
print att_info
print

print "Repeated attractors:"

krange = np.arange(len(updated_pairs))

for k in krange:
    print updated_pairs[k]
print    
