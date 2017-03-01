# fastOC
This package implements an Orthologous Clustering approach to identify conserved and lineage specific co-expressed gene modules across multiple species. In addition, we illustrate various aspects of data manipulation, filtering and visualization of large transcriptomic experiments.

fastOC can be installed in R using the following commands:

````{r}
    #install dependencies
    install.packages(c("igraph", "Matrix", "WGCNA", 
                       "reshape2", "fastcluster", "dynamicTreeCut"),
                     dependencies = TRUE)

    #install fastOC
    require(devtools);
    install_github("mzinkgraf/fastOC");
````

Information about the original program OrthoClust can be found at Yan & Wang et al. (2008). The main differences between fastOC and OrthoClust are that (1) fastOC can work with many species and (2) uses the Louvain cluster method to increase computational efficiency.

Yan K-K, Wang D, Rozowsky J, Zheng H, Cheng C, Gerstein M. 2014. OrthoClust: an orthology-based network framework for clustering data across multiple species. Genome Biology 15(8): R100.

Vincent DB, Jean-Loup G, Renaud L, Etienne L. 2008. Fast unfolding of communities in large networks. Journal of Statistical Mechanics: Theory and Experiment 2008(10): P10008.

