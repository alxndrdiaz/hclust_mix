# Python 2.7.16

#  evaluates how many unique attractors were detected and generates a file for each one ["A.tab" files]
#  IDs for these files are assigned based on the first observed sample that converged to an attractor ("A1", "A5", etc) 
#  for each attractor separates genes by state (-1,0,+1), ["low", "zero", "high" .ids files] 
#  generates a table containing each sample and its associated attractor ["samples_attractors.tab"]
#  general summary table contains the number of genes, samples and unique attractors ["general_summary.txt"]
#  attractor summary table contains total number of samples and genes by state for each attractor ["attractor_summary.txt"]
#  "attractors_heatmap.png" shows samples clustered by the attractors, each attractor is marked with a color at the top 
#  all .tab, and .ids files are moved to "attractor_results/" directory


import os
import shutil
import glob
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt


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
'unique_attractors': [unique_attractors] })
general_summary.to_csv('general_summary.txt', index=False, header=True, sep="\t")


# extracts genes by state (-1,0,+1) for each attactor:
att_IDS = []
high_genes = []
low_genes = []
zero_genes = []
cols=list(U_attractors.columns.values)
for column in U_attractors:
    high_states = U_attractors[column][ U_attractors[column] == +1 ]
    low_states = U_attractors[column][ U_attractors[column] == -1 ]
    zero_states = U_attractors[column][ U_attractors[column] == 0 ] 
    outname =  'A' + str( 1 + cols.index(column) ) 
    att_IDS.append(outname)
    high_genes.append( len(high_states) )
    low_genes.append( len(low_states) )
    zero_genes.append( len(zero_states) )  
    U_attractors[column].to_csv( outname + '.tab', index=True, header=True, sep="\t")  
    high_states.to_csv(  outname + '_genes_high.ids', index=True, header=False, sep="\t")
    low_states.to_csv( outname + '_genes_low.ids', index=True, header=False, sep="\t")
    zero_states.to_csv(  outname + '_genes_zero.ids', index=True, header=False, sep="\t") 

    
# finds samples that converged to the same attractor: 
atts_samples = []
atts_sorted = [None]*len(attractors.columns)
for column in U_attractors:
    att_by_sample = []
    for nsample in attractors: 
        if U_attractors[column].equals( attractors[nsample] ) == True:
           att_by_sample.append( nsample )
           # gets position of the converged sample 
           msample = attractors.columns.get_loc(nsample)
           # saves the attractor label associated to sample
           atts_sorted[msample] = 'A' + str( 1 + cols.index(column) )
    atts_samples.append(att_by_sample)
# number of samples that converged to each attractor
nsamples_attractor = [] 
for atsample in atts_samples: 
     nsamples_attractor.append( len(atsample) )
# table that contains each sample and its attractor
samples_to_attractors = pd.DataFrame( 
{'sample': list(attractors.columns),
'attractor': atts_sorted} )
samples_to_attractors = samples_to_attractors.reindex( columns = ['sample','attractor'] )
samples_to_attractors.to_csv('samples_attractors.tab', index=False, header=True, sep="\t")


# generates attractor summary table: 
attractor_summary = pd.DataFrame({
'attractor': att_IDS,
'genes_high': high_genes,
'genes_low': low_genes, 
'genes_zero': zero_genes, 
'nsamples': nsamples_attractor })
attractor_summary.to_csv('attractor_summary.txt', index=False, header=True, sep="\t")

# plots attractors:
# colors for column 
samples_palette = sn.color_palette( palette='husl', n_colors= len(att_IDS) )	 
colors_att = dict( zip(att_IDS,samples_palette)  )
colors_list = map(colors_att.get, atts_sorted)
# seaborn clustermap
att_heatmap = sn.clustermap( attractors,  
annot=False, linewidths=.005, 
vmin = -1, vmax=1,  cbar_kws={'ticks':[-1,0,1]},
col_colors = colors_list, cmap='vlag', 
xticklabels=False, yticklabels=False )
# add columns legend (associated to attractors)  
for label in att_IDS:
    att_heatmap.ax_col_dendrogram.bar(0,0,color=colors_att[label],label=label,linewidth=0) 
att_heatmap.ax_col_dendrogram.legend(loc='best', ncol=3)
# title 
plotitle = str(total_samples) + ' samples clustered by N = ' + str( len(att_IDS) ) + ' attractors' + ', genes = ' + str(total_genes) 
att_heatmap.fig.suptitle(plotitle, fontsize=20)
# saves figure
att_heatmap.plot
plt.savefig('attractors_heatmap.png', format='png', dpi=300)


# move all the results to a directory: 
os.mkdir('attractor_results')
idsfiles = glob.glob(currentdirectory + "/*.ids")
tabfiles = glob.glob(currentdirectory + "/*.tab")
results = idsfiles + tabfiles 
for result in results: 
    shutil.move(result, 'attractor_results')

