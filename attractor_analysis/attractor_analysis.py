# Python 2.7.6
# libraries: numpy 1.8.2; pandas 0.17.1; scipy 0.13.3

# takes the output attractors from hclust_mix in .ats format and analyze overrepresented and  underrepresented genes, which are the predictions of the 
# identified attractors. It also compares attractors in order to determine if they are different and how they are structured, provides output files 
# for further analysis (.txt files, including files containing genes per binary state: -1, +1 ; and a summary file).

import os
import shutil
import glob
import pandas as pd


# detects current working directory and assigns it to the currentdirectory string:
currentdirectory = os.getcwd()

# reads attractor file:
att_file = ''.join( glob.glob(currentdirectory + "/*.ats") )
attractors = pd.read_csv(att_file, index_col=0, header=0, sep="\t")
# total genes and samples used to search attractors
total_genes = attractors.shape[0]
total_samples = attractors.shape[1]

# finds unique attractors:
U_attractors = attractors.T
U_attractors = U_attractors.drop_duplicates(keep='first')
U_attractors = U_attractors.T
# total numbber of unique attractors 
unique_attractors = U_attractors.shape[1]

# generates general summary table:
general_summary = pd.DataFrame({
'total_genes': [total_genes],
'total_samples': [total_samples],
'unique_attractors': [unique_attractors]
})
general_summary.to_csv('general_summary.txt', index=False, header=True, sep="\t")


# extracts genes by state (-1,0,+1) for each attactor:
att_IDS = []
high_genes = []
low_genes = []
zero_genes = []
for column in U_attractors:
    high_states = U_attractors[column][ U_attractors[column] == +1 ]
    low_states = U_attractors[column][ U_attractors[column] == -1 ]
    zero_states = U_attractors[column][ U_attractors[column] == 0 ] 
    outname = 'A'+ str(column)
    att_IDS.append( outname )
    high_genes.append( len(high_states) )
    low_genes.append( len(low_states) )
    zero_genes.append( len(zero_states) )  
    U_attractors[column].to_csv( outname + '.tab', index=True, header=True, sep="\t")  
    high_states.to_csv(  outname + '_genes_high.ids', index=True, header=False, sep="\t")
    low_states.to_csv( outname + '_genes_low.ids', index=True, header=False, sep="\t")
    zero_states.to_csv(  outname + '_genes_zero.ids', index=True, header=False, sep="\t") 

        
# find samples that converged to the same attractor:
atts_samples = []
atts_sorted = [None]*len(attractors.columns)

for column in U_attractors:
    att_by_sample = []
    for nsample in attractors: 
        msample = int(nsample) - 1
        if U_attractors[column].equals( attractors[nsample] ) == True:
           att_by_sample.append( nsample )
           atts_sorted[msample] = 'A' + column
    atts_samples.append(att_by_sample)

print
for atsample in atts_samples: 
    print 'total samples converged to attractor = ', len(atsample); print 
print

samples_to_attractors = pd.DataFrame( 
{'sample': list(attractors.columns),
'attractor': atts_sorted} )

samples_to_attractors = samples_to_attractors.reindex( columns = ['sample','attractor'] )
samples_to_attractors.to_csv('samples_attractors.tab', index=False, header=True, sep="\t")


# move all the results to a directory: 
os.mkdir('attractor_results')
idsfiles = glob.glob(currentdirectory + "/*.ids")
tabfiles = glob.glob(currentdirectory + "/*.tab")
results = idsfiles + tabfiles 
for result in results: 
    shutil.move(result, 'attractor_results')





    


