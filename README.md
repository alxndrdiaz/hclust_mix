# What is hclust_mix? 

h: Hopfield networks, clust: clustering, mix: attractor analysis and some tools to process the data.

These tools are adapted for preprocessing and analysis of time-course developmental transcriptomic data using hclust 1.0, 
the algorithm was proposed by Stefan Maetschke and Mark Ragan to characterize cancer subtypes, [Maetschke and Mark Ragan, Bioinformatics (2014)](https://academic.oup.com/bioinformatics/article/30/9/1273/234782). The original version of hclust along with a short tutorial are available at:
[http://bioinformatics.org.au/tools/hclust/](http://bioinformatics.org.au/tools/hclust/). The idea behind their algorithm was to "demonstrate the usage 
of Hopfield networks for clustering, feature selection and network inference" [see tutorial by Mark Ragan](http://bioinformatics.org.au/tools/hclust/)  applied to transcriptomic data. Specifically it aims to model differentiated cell states as attractor states of a Hopfield Network, and was tested in several datasets including Haemopoiesis and human stem cell differentiation [Fard et al. npj Syst Biol Appl 2, (2016)](https://www.nature.com/articles/npjsba20161).


## 1. Waddington epigenenetic landscape 

[Conrad Hal Waddington](https://en.wikipedia.org/wiki/C._H._Waddington) in his 1957 book *The strategy of the genes* (London: George Allen & Unwin) proposed a metaphor to explain how a pluripotent cell becomes a differentiated cell. In his metaphor the pluripotent cell is like a ball at the top of a hill, while this cell differentiates it moves down through the rugged landscape of the hill until it reaches the bottom, that is a fully differentiated state. Genes modify this lanscape to allow only certain paths to exist, in such a way that there is only a limited number of possible outcomes.    

![The Waddington Epigenetic Landscape](https://web.archive.org/web/20050902020936im_/http://zygote.swarthmore.edu/gene33.GIF)

Finally, a very important aspect of this metaphor is that if cells differentiated in this way, most cells would be observed to differentiate in a very deterministic way and this process should be regulated by a very stable expression of lineage-specific genes. Some authors have proposed that this might not 
be the case, at least considering new evidence from single-cell transcriptomes, see for example  multilineage priming effect in frog and fish from [Klein et al. Science (2018)](https://science.sciencemag.org/content/360/6392/eaar5780). Mention references that discuss this idea.


## 2.  What is a Hopfield network and why was used in this context?

In 1982 John Hopfield proposed a model of neural network to understand how neurons can storage information [Hopfield PNAS (1982)](https://www.pnas.org/content/79/8/2554). Nodes in the network are neurons and edges between them are called weights and can be updated according to a rule that defines how they interact to  change their states. These type of network can converge from an initial state to a stable state called an attractor, this convergence is achieved by minimizing an energy function. The model proposed by  [Fard et al. npj Syst Biol Appl 2, (2016)](https://www.nature.com/articles/npjsba20161) borrows this idea, but instead of neurons, nodes are genes and weights represent co-expresion, the initial pluripotent states can converge to an attractor that represents a differentiated cell state, in this way Waddington Epigenetic Landscape is "not just a colourful metaphor".   

### 2.1 Hopfield network model in hclust

### 2.2 Workflow


`
## 3. Limitations

While this approach is quite good ilustrating how differentiated cell states tend to stabilize compared to transient cell states, it might not be the the option to choose if you wish to understand or model cell differentiation in general, unless you perform a benchmark analysis to determine that this method outperforms other methods in predicting differentiated cell states. So, if you would like to start learning about cell differentiation models hclust is a good starting point. However, if you are interested in questions related to regulation of cell differentiation using single-cell transcriptomics datasets, there are other packages such as [RNA velocity](http://velocyto.org/), or [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) that provide useful analysis tools. 



## 4. Usage

Describe what each script makes, the output content and how to run each one. 


