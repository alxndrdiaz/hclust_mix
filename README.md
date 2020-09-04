# What is hclust_mix? 

h: Hopfield networks, clust: clustering, mix: attractor analysis and some tools to process the data.

These tools are adapted for preprocessing and analysis of time-course developmental transcriptomic data using hclust 1.0, 
the algorithm was proposed by Stefan Maetschke and Mark Ragan to characterize cancer subtypes, [Maetschke and Mark Ragan, Bioinformatics (2014)](https://academic.oup.com/bioinformatics/article/30/9/1273/234782). The original version of hclust along with a short tutorial are available at:
[http://bioinformatics.org.au/tools/hclust/](http://bioinformatics.org.au/tools/hclust/). The idea behind their algorithm was to "demonstrate the usage 
of Hopfield networks for clustering, feature selection and network inference" [see tutorial by Mark Ragan](http://bioinformatics.org.au/tools/hclust/)  applied to transcriptomic data. Specifically it aims to model differentiated cell states as attractor states of a Hopfield Network, and was tested in 12 single-cell data sets including Haematopoiesis and human stem cell differentiation, [Fard et al. npj Syst Biol Appl 2, (2016)](https://www.nature.com/articles/npjsba20161).

## 1.  Hopfield network model and hclust workflow

In 1982 John Hopfield proposed a model of neural network to understand how neurons can storage information. Nodes in the network are neurons with binary states (-1,+1) and edges between them are called weights, these can be updated according to a rule called Hebb's rule. Neuron states can be updated using these weights. An important property of  this network is that it can converge from an initial state to a stable state called an attractor, this convergence is achieved by minimizing an "energy function", [Hopfield PNAS (1982)](https://www.pnas.org/content/79/8/2554); [Rojas Springer-Verlag (1996)](http://page.mi.fu-berlin.de/rojas/neural/index.html.html). This type of neural network can be used for retrieving patterns that have been corrupted by noise, a good illustrative example can be found in [this Jupyter notebook by Filippo Galli](https://github.com/philipjk/genetic_algorithm_optimization_sklearn-based/blob/master/hopfield_networks.ipynb). 

The model proposed by Maetschke and Mark Ragan borrows this idea, but instead of neurons, nodes are genes and weights represent co-expresion, the initial pluripotent states can converge to an attractor that represents a differentiated cell state. Additionally they use a ternary set of node states (-1, 0, +1),  [Maetschke and Mark Ragan, Bioinformatics (2014)](https://academic.oup.com/bioinformatics/article/30/9/1273/234782); [Fard et al. npj Syst Biol Appl 2, (2016)](https://www.nature.com/articles/npjsba20161). The basic steps of hclust workflow are: 

1. Normalization of expression matrix and feaure selection (genes with highest variance). 
2. Discretizing normalized expression matrix. 
2. Compute initial weight matrix.
3. Iterate to minimize energy function until convergence to attactor states. 
4. Visualization of energy landscape in PCA space.



## 2. Usage

Describe what each script makes, the output content and how to run each one. 

### 2.2 Requirements
The script requires Python 2.7, and the libraries ......, If you would like to make a virtual environment use: bla bla bla 


### 2.2 Expression matrix format



### 2.3 How to use



### 2.4 Results

![relaxation](test_results_yeoh_reduced/1_relaxation_state_matrix.png)

![pruning](test_results_yeoh_reduced/3_pruning_threesholds.png)

![contour](test_results_yeoh_reduced/6_PCA_contour_plot.png)

![barplot](test_results_yeoh_reduced/attractors_barplot.png)

![dendrogram](test_results_yeoh_reduced/attractors_dendrogram.png)

![heatmap](test_results_yeoh_reduced/attractors_heatmap.png)



## 3. Limitations

While this approach is quite good ilustrating how differentiated cell states tend to stabilize compared to transient cell states, it might not be the the option to choose if you wish to understand or model cell differentiation in general, unless you perform a benchmark analysis to determine that this method outperforms other methods in predicting differentiated cell states. So, if you would like to start learning about cell differentiation models hclust is a good starting point. However, if you are interested in questions related to regulation of cell differentiation using single-cell transcriptomics datasets, there are other packages such as [RNA velocity](http://velocyto.org/), or [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) that provide useful analysis tools. 


## 4. Final remarks: Waddington epigenenetic landscape 

[Conrad Hal Waddington](https://en.wikipedia.org/wiki/C._H._Waddington) in his 1957 book *The strategy of the genes* (London: George Allen & Unwin) proposed a metaphor to explain how a pluripotent cell becomes a differentiated cell. In his metaphor the pluripotent cell is like a ball at the top of a hill, while this cell differentiates it moves down through the rugged landscape of the hill until it reaches the bottom, that is a fully differentiated state. Genes modify this lanscape to allow only certain paths to exist, in such a way that there is only a limited number of possible outcomes.    

#![Waddington](https://jeb.biologists.org/content/jexbio/218/6/816/F2.large.jpg) 
Figure 1. Waddington Epigenetic Landscape. A and B represent different differentiation paths, [Noble J. Exp. Biol. (2015)](https://jeb.biologists.org/content/218/6/816).

Finally, a very important aspect of this metaphor is that if cells differentiated in this way, most cells would be observed to differentiate in a very deterministic way and this process should be regulated by a very stable expression of lineage-specific genes. Some authors have proposed that this might not 
be the case, at least considering new evidence from single-cell transcriptomes, see for example  multilineage priming effect in frog and fish from [Klein et al. Science (2018)](https://science.sciencemag.org/content/360/6392/eaar5780). For a detailed discussion about Waddington's ideas see for example [Huang BioEssays (2011)](https://onlinelibrary.wiley.com/doi/abs/10.1002/bies.201100031) and [Noble J. Exp. Biol. (2015)](https://jeb.biologists.org/content/218/6/816).

## References



