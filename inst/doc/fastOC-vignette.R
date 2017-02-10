## ----results='hide', message=FALSE, warning=FALSE, eval=FALSE------------
#      require(devtools);
#      install_github("mzinkgraf/fastOC");

## ----results='hide', message=FALSE, warning=FALSE------------------------
    require(fastOC);
    options(scipen = 0);
    options(stringsAsFactors = FALSE);
    #change working directory
    setwd(".")

## ------------------------------------------------------------------------
    #load rpkm data sets
      load("../Data/eugra_potri_sapur_rpkm.rdata")
      names(rpkm)

## ---- results = "asis",echo = FALSE--------------------------------------
  pander::pandoc.table(rpkm$potri[1:4,1:2])

## ------------------------------------------------------------------------
  rpkm_var <- filterVariance(rpkm, variance = 0.1)

## ------------------------------------------------------------------------
    GeneMeta <- createGeneMeta(rpkm_var)

## ---- results = "asis",echo = FALSE--------------------------------------
  x<-GeneMeta[c(1,2,28779,64104),]
  row.names(x)<-NULL
  pander::pandoc.table(x)

## ------------------------------------------------------------------------
    ortho = list()
      ortho$potri_eugra <- read.table("../Data/potri_eugra.txt", sep = " ")
      ortho$potri_sapur <- read.table("../Data/potri_sapur.txt", sep = " ")
      ortho$sapur_eugra <- read.table("../Data/sapur_eugra.txt", sep = " ")

## ---- results = "asis",echo = FALSE--------------------------------------
  names(ortho$potri_eugra)[1:2]<-c("GeneA","GeneB")
  pander::pandoc.table(ortho$potri_eugra[1:4,1:2])

## ---- eval=FALSE---------------------------------------------------------
#      gene_edgelist<-getEdgelist(rpkm_var, GeneMeta, top = 5)

## ---- eval=FALSE---------------------------------------------------------
#      ortho_weights <- getOrthoWeights(ortho, GeneMeta, couple_const = 1)
#  

## ---- eval=FALSE---------------------------------------------------------
#    names(ortho_weights) <- names(gene_edgelist)
#    combined_out <- rbind(gene_edgelist, ortho_weights)
#  
#    #write file to disk
#    save(combined_out,file="../Data/Louvain_input.rdata")
#  

## ---- echo = FALSE-------------------------------------------------------
    load("../Data/Louvain_input.rdata");
    load("../Data/MultiSppTrees.rdata");
    load("../Data/Occurance_sparse.rdata");

## ----results='hide', message=FALSE, warning=FALSE, eval=FALSE------------
#      nRuns = 100
#      results <- louvain(combined_out, nRuns)
#      save(results, file="../Data/igraph_louvain.results")

## ---- eval=FALSE---------------------------------------------------------
#      occurance <- filterCommunityAssign(results, minMem = 10)
#      save(occurance, file = "Occurance_sparse.rdata")

## ---- eval=FALSE---------------------------------------------------------
#      MultiSpp_trees <- multiSppHclust(occurance, nRuns, GeneMeta)

## ------------------------------------------------------------------------
    MultiSpp_modules <- multiSppModules(MultiSpp_trees, GeneMeta, 
         minModuleSize = c(800,800,800), 
         cut = c(0.5,0.5,0.5))

    #add module assignemnts to columns 3 of GeneMeta

    GeneMeta[,3]<-do.call(c,lapply(seq_along(MultiSpp_modules), function(y, i) { y[[i]] }, y=MultiSpp_modules))

## ----fig.show='hold',results='hide', message=FALSE, warning=FALSE--------
    #Create heatmap of co-appearance modules across all species
    par(fig=c(0,1,0,1))
    plot_MultiSpp(GeneMeta = GeneMeta, order = MultiSpp_trees$order, 
                  sb=12, CA_keep = occurance, remove_0 = F, lwd = 1, cex = 0.5)

    #Add heatmap color bar
    par(fig=c(0,0.6,0.6,0.9), new=TRUE)
    color.bar(lut = colorRampPalette(c("white", "lightyellow", "red","black"))(n = 100),
              min = 0, max = 1, nticks=5,
              title = "Co-Appearance")

