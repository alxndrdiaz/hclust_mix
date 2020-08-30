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
attractors = pd.read_csv(att_file, index_col=0, header=[0,1], sep="\t")
# gets both column indices and reuse only the first one
sample_labels = attractors.columns.get_level_values(0).tolist() 
type_labels = attractors.columns.get_level_values(1).tolist() 
attractors.columns = sample_labels

# total genes and samples used to search attractors
total_genes = attractors.shape[0]
total_samples = attractors.shape[1]

# finds unique attractors:
U_attractors = attractors.T
U_attractors = U_attractors.drop_duplicates(keep='first')
U_attractors = U_attractors.T
# changes column labels
cols=list(U_attractors.columns.values)
ucolumns = [] 
for column in U_attractors:
    ucolumns.append( 'A' + str( 1 + cols.index(column) ) ) 
U_attractors.columns = ucolumns
# total numbber of unique attractors 
unique_attractors = U_attractors.shape[1]


# generates general summary table:
general_summary = pd.DataFrame({
'total_genes': [total_genes],
'total_samples': [total_samples],
'unique_attractors': [unique_attractors] })
general_summary.to_csv('general_summary.txt', index=False, header=True, sep="\t")


# extracts genes by state (-1,0,+1) for each attactor:
high_genes = []
low_genes = []
zero_genes = []
for column in U_attractors:
    high_states = U_attractors[column][ U_attractors[column] == +1 ]
    low_states = U_attractors[column][ U_attractors[column] == -1 ]
    zero_states = U_attractors[column][ U_attractors[column] == 0 ] 
    high_genes.append( len(high_states) )
    low_genes.append( len(low_states) )
    zero_genes.append( len(zero_states) )  
    U_attractors[column].to_csv( column + '.tab', index=True, header=True, sep="\t")  
    high_states.to_csv(  column + '_genes_high.ids', index=True, header=False, sep="\t")
    low_states.to_csv( column + '_genes_low.ids', index=True, header=False, sep="\t")
    zero_states.to_csv(  column + '_genes_zero.ids', index=True, header=False, sep="\t") 

    
# finds samples that converged to the same attractor: 
atts_samples = []
atts_sorted = [None]*len(attractors.columns)
stypes_sorted = [None]*len(attractors.columns)
for column in U_attractors:
    att_by_sample = []
    for nsample in attractors: 
        if U_attractors[column].equals( attractors[nsample] ) == True:
           att_by_sample.append( nsample )
           # gets position of the converged sample 
           msample = attractors.columns.get_loc(nsample)
           # saves the attractor label and type associated to sample
           atts_sorted[msample] = column
           stypes_sorted[msample] = type_labels[msample]
    atts_samples.append(att_by_sample)
# number of samples that converged to each attractor
nsamples_attractor = [] 
for atsample in atts_samples: 
     nsamples_attractor.append( len(atsample) )
# table that contains each sample and its attractor
samples_to_attractors = pd.DataFrame( 
{'sample': list(attractors.columns),
'type': stypes_sorted,
'attractor': atts_sorted} )
samples_to_attractors = samples_to_attractors.reindex( columns = ['sample','type','attractor'] )
samples_to_attractors.to_csv('samples_attractors.tab', index=False, header=True, sep="\t")


# generates attractor summary table: 
att_IDS=list(U_attractors.columns.values)
attractor_summary = pd.DataFrame({
'attractor': att_IDS,
'genes_high': high_genes,
'genes_low': low_genes, 
'genes_zero': zero_genes, 
'nsamples': nsamples_attractor })
attractor_summary.to_csv('attractor_summary.txt', index=False, header=True, sep="\t")


# plots attractors heatmap:
# colors for column 
samples_palette = sn.color_palette( palette='husl', n_colors= len(att_IDS) )	 
colors_att = dict( zip(att_IDS,samples_palette)  )
colors_list = map(colors_att.get, atts_sorted)
# seaborn clustermap
att_heatmap = sn.clustermap( attractors,
figsize = (8,6),  
annot=False, linewidths=.00005, 
vmin = -1, vmax=1,  cbar_kws={'ticks':[-1,0,1]},
col_colors = colors_list, cmap='vlag', 
xticklabels=False, yticklabels=False )
# title 
plotitle = str(total_samples) + ' samples clustered by N = ' + str( len(att_IDS) ) + ' attractors' + ', genes = ' + str(total_genes) 
att_heatmap.fig.suptitle(plotitle, fontsize=18)
# saves figure
att_heatmap.plot
plt.savefig('attractors_heatmap.png', format='png', dpi=300); plt.clf()


# plots attractors dendrogram:  
from scipy.cluster import hierarchy
# computes distance between samples
datts = hierarchy.linkage(U_attractors.T, metric='euclidean')
# plots the dendrogram
plt.title('Hierarchical Clustering Dendrogram')
plt.ylabel('attractors')	
plt.xlabel('distance [Euclidean]')
hierarchy.set_link_color_palette(None)
hierarchy.dendrogram(datts, labels=U_attractors.columns, leaf_rotation=0, orientation='left')
plt.savefig('attractors_dendrogram.png', format='png', dpi=300); plt.clf()


# plots attractor stacked bar plot for sample type content
att_content =  samples_to_attractors[ ['type','attractor'] ]
att_counts = pd.crosstab(index=att_content['attractor'], columns=att_content['type'])
att_counts.to_csv('attractor_content_summary.txt', index=True, header=True, sep="\t")
att_fracs = att_counts.apply(lambda x: x/sum(x), axis=1)
att_fracs.plot.bar(stacked=True)
plt.title('Distribution of samples in attractors')
plt.ylabel('fraction')
plt.legend(loc='best')
plt.savefig('attractors_barplot.png', format='png', dpi=300); plt.clf()

# move all the results to a directory: 
os.mkdir('attractor_results')
idsfiles = glob.glob(currentdirectory + "/*.ids")
tabfiles = glob.glob(currentdirectory + "/*.tab")
results = idsfiles + tabfiles 
for result in results: 
    shutil.move(result, 'attractor_results')
