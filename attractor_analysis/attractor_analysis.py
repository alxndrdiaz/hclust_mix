# Python 2.7.6
# libraries: numpy 1.8.2; pandas 0.17.1; scipy 0.13.3

# takes the output attractors from hclust_mix in .ats format and analyze overrepresented and  underrepresented genes, which are the predictions of the 
# identified attractors. It also compares attractors in order to determine if they are different and how they are structured, provides output files 
# for further analysis (.txt files, including files containing genes per binary state: -1, +1 ; and a summary file).

import os
import glob
import time
import numpy as np
import pandas as pd

# detects current working directory and assigns it to the currentdirectory string
currentdirectory = os.getcwd()

# reads attractor file
att_file = ''.join( glob.glob(currentdirectory + "/*.tab") )
attractors = pd.read_csv(att_file, index_col=0, header=0, sep="\t")
print
print 'genes, samples = ', attractors.shape
print

# finds unique attractors
U_attractors = attractors.T
U_attractors = U_attractors.drop_duplicates(keep='first')
U_attractors = U_attractors.T
print 'total unique attractors = ', U_attractors.shape[1]
print

# iterate to extract genes by state (-1,+1) for each attactor
for column in U_attractors:
    low_states = U_attractors[column][ U_attractors[column] == -1 ]
    high_states = U_attractors[column][ U_attractors[column] == +1 ]
    outname = str(column) + '_sample_attractor'
    low_states.to_csv( outname + '_low.ids', index=True, header=False, sep="\t")
    high_states.to_csv(  outname + '_high.ids', index=True, header=False, sep="\t") 


