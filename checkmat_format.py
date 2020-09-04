
import glob
import pandas as pd

# directory containing example matrix
DIR = 'test_yeoh_reduced'

# reads attractor file:
matfile = ''.join( glob.glob(DIR + "/*.tsv") )
# reads only column indices to check format 
tidx = pd.read_csv(matfile, index_col=0, header=[0,1], nrows=1, skipinitialspace=True, sep="\t")
# loads file with pd.MultiIndex with n=2 level, expected format
if isinstance(tidx.columns, pd.MultiIndex):
   sample_labels = tidx.columns.get_level_values(0).tolist() 
   type_labels = tidx.columns.get_level_values(1).tolist() 
   print; print 'sample labels: \n \n', sample_labels; print
   print; print 'type labels: \n \n', type_labels; print
   print 'Matrix format is correct. \n'   
else: 
   print '\n Format is incorrect. Please provide unique labels for samples. \n'  

