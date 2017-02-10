#these functions are designed to work with the orthologous clustering methods
#that are being developed for the analysis of multiple angiosperm tree species

#'Filter Louvain community assignemnts
#'
#'This function removes communities from the sparse matrix that have a minimum number of gene members
#'
#' @usage filterCommunityAssign(results,minMem=10)
#' @param results Data frame containing the community assignments for each genes and contains nRuns+1 columns
#' @param minMem Minimum number of genes in a community to be considered significant. Default = 10
#' @import methods
#' @importMethodsFrom Matrix colSums
#' @return Returns a sparse matrix containing to occurance of each gene in each community. 1 = present and 0 = absent
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
filterCommunityAssign <- function(results, minMem=10)
{
  nRuns=ncol(results)-1
  for(j in 1:nRuns)
  {
    tmp=results[,c(1,j+1)]
    tmp[,3]=1
    names(tmp)=c("V1","V2","V3")
    if(j==1)
    {
      d_sparse = sparseMatrix(as.integer(tmp[,1]), as.integer(tmp[,2]), x = tmp$V3)
    } else {
      data.sparse = sparseMatrix(as.integer(tmp[,1]), as.integer(tmp[,2]), x = tmp$V3)
      d_sparse<-cbind(d_sparse,data.sparse)
    }
  }

  #full module assignments
  keep=which(Matrix::colSums(d_sparse)>minMem)
  d_keep<-d_sparse[,keep]
  return(d_keep)
}

#'Filter expression based on variance
#'
#'This function filters rpkm expression data based on user defined variance threshold.
#' @usage filterVariance(rpkm, variance=0.1)
#' @param rpkm List object containing the rpkm values for each species.
#' @param variance Variance threshold.
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
filterVariance = function(rpkm, variance=0.1)
{
  if(is.list(rpkm))
  {
    for(v in 1:length(rpkm))
    {
      vari<-apply(rpkm[[v]],1,var)
      index<-which(vari>variance)
      rpkm[[v]]<-rpkm[[v]][index,]
      return(rpkm)
    }
  } else {
    return("Object is not a list")
  }
}

#'Create gene metadata object
#'
#'This function generates a data frame that contains the species and gene information, and columns to store module results.
#' @usage createGeneMeta(rpkm)
#' @param rpkm List object containing the rpkm values for each species.
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
createGeneMeta = function(rpkm)
{
  n<-names(rpkm)

  nGenes<-cumsum(do.call(c,lapply(seq_along(rpkm), function(y, i) {  length(dimnames(y[[i]])[[1]] ) }, y=rpkm)))

  gene_names<-do.call(rbind,lapply(seq_along(rpkm), function(y, n, i) { cbind(n[[i]], dimnames(y[[i]])[[1]] ) }, y=rpkm, n=names(rpkm)))
  gene_names<-as.data.frame(gene_names)
  gene_names[,3]<-0
  gene_names[,4]<-seq(1:nGenes[length(nGenes)])
  names(gene_names)<-c("species","gene","modules","ID")
  return(gene_names)
}

#'Combine multiple files
#'
#'Read and combine multiple results files from a folder. The results are combined using cbind and the row order in each file must be the same.
#' @usage multMerge(mypath, pattern="*\\.out", header=FALSE, sep=" ")
#' @param mypath Path to folder.
#' @param pattern Unique pattern that matches files to be be read and combined using perl style regular expressions.
#' @param header A logical value indicateing if the file contains a header row. Default = FALSE
#' @param sep The field separator character.
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @seealso  \code{\link[base]{list.files}}
#' @export
multMerge = function(mypath, pattern="*\\.out", header=FALSE, sep=" ")
{
  filenames=list.files(path=mypath, pattern=pattern, full.names=TRUE)
  index<-which(file.info(filenames)$size>0)
  datalist = lapply(filenames[index], function(x) {read.table(file=x,header=header, sep=sep)})
  output=Reduce(function(x,y) {cbind(x,y[,2])}, datalist)
  #output=do.call(rbind.data.frame, datalist)
  return(output)
}

#'Combine multiple htseq files using merge
#'
#'Read and combine tab deliminted htseq results files from a folder using the merge function.
#' @usage multMergeHTseq(mypath, pattern="*\\.htseq", byY="gene", RegEx="(\\w+)\\.htseq\\.txt", Replace="\\1", sep="\t")
#' @param mypath Path to folder.
#' @param pattern Unique pattern that matches files to be read and merged using perl style regular expressions.
#' @param byY Column name that should be used in the merge function.
#' @param RegEx Perl style regular expression that matches pattern in the filename. Used to extract library name from filename. Example "(\\w+)\\.htseq\\.txt" matches Library1 in filename Library1.htseq.txt
#' @param Replace Perl style replacement. Default = "\\1"
#' @param sep The field separator character.
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @seealso  \code{\link[base]{list.files}}
#' @export
#combine htseq results into a single data frame
multMergeHTseq = function(mypath, pattern="*\\.htseq", byY="gene",
                          RegEx="(\\w+)\\.htseq\\.txt", Replace="\\1", sep="\t")
{
  filenames=list.files(path=mypath, pattern=pattern, full.names=TRUE)
  datalist = lapply(filenames, function(x) {read.table(file=x, header=F, sep=sep)})
  output=Reduce(function(x,y) {suppressWarnings(merge(x,y,by.x="V1",by.y="V1",all=T))}, datalist)
  #get short names and create header line
  nm=list.files(path=mypath, pattern=pattern, full.names=F)
  names(output)=c(byY, sub(RegEx, Replace, nm,perl=T))
  return(output)
}

#' Correlation matrix to list of maximum
#'
#'For each row in a matrix determine index values of the top (default=5) most correlated genes. An internal function used by getEdgelist()
#' @usage weighted2rankList(m, top=5)
#' @param m A maxtrix of correlation values
#' @param top An integer specifying the number values to return. Default = 5
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @seealso  \code{\link{getEdgelist}}
#' @export
weighted2rankList<-function(m, top=5)
{
  tmp<-apply(m,1,function(x) order(x,decreasing=TRUE)[1:top])
  tmp<-as.data.frame(tmp)
  names(tmp)<-seq(1,ncol(tmp),1)
  out<-melt(tmp,factorsAsStrings = TRUE)
  out[,1]<-as.numeric(as.character(out[,1]))
  #remove duplicates
  cat1<-paste(out[,1],out[,2],sep="_")
  cat2<-paste(out[,2],out[,1],sep="_")
  dup<-which(cat2 %in% cat1)
  return(out[-dup,])
}

#' Create an edgelist from rpkm expression values
#'
#'Generate an edgelist for each gene from a list of rpkm values and return each gene and its top most correlated neighbors
#' @usage getEdgelist(rpkm, GeneMeta, top=5, weight=1, nThreads = 3)
#' @param rpkm A list object where each element in the list is a data frame of rpkm values for each species
#' @param GeneMeta Data frame that contains the project metadata. See \link{createGeneMeta}
#' @param top An integer specifying the number of neighbors for each gene that should be printed to the edgelist. Default = 5
#' @param weight The edge weight between genes. Default = 1
#' @param nThreads The number of multiple threads that should be used to calculate the correlation matrix. Default = 3
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
getEdgelist<-function(rpkm,GeneMeta,top=5,weight=1,nThreads = 3)
{
  if(is.list(rpkm) & length(rpkm)>0)
  {
    nGenes<-cumsum(do.call(c,lapply(seq_along(rpkm), function(y, i) {  length(dimnames(y[[i]])[[1]] ) }, y=rpkm)))

    full_edgelist<-data.frame(matrix(nrow = 0,ncol = 2))
    for(q in 1:length(rpkm))
    {

      datExpr0<-t(rpkm[[q]])

      Cor<-corFast(datExpr0,nThreads = nThreads)
      diag(Cor)<-0
      collectGarbage()
      edgelist<-weighted2rankList(Cor,top=top)

      if(q==1)
      {
        edgelistB<-edgelist
        edgelistB[,1]<-edgelist[,1]
        edgelistB[,2]<-edgelist[,2]
      } else {
        edgelistB<-edgelist
        edgelistB[,1]<-edgelist[,1]+nGenes[q-1]
        edgelistB[,2]<-edgelist[,2]+nGenes[q-1]
      }

      collectGarbage()

      full_edgelist<-rbind(full_edgelist,edgelistB)
      #write.table(edgelist,file=paste("~/Documents/scripts/louvain-generic/PSE_R_var/",names(rpkm)[q] ,"_edgelist.txt",sep=""),sep="\t",row.names = F,col.names = F,quote = F)

      rm(Cor)
      collectGarbage()
    }

    full_edgelist[,3]<-weight
    Rfull_edgelist<-full_edgelist[,c(2,1,3)]
    names(Rfull_edgelist)<-names(full_edgelist)
    fullE<-rbind(full_edgelist,Rfull_edgelist)
    names(fullE)<-c("variable","value","weight")
    return(fullE[order(fullE[,1],fullE[,2]),])
  } else {
    return(print("Check rpkm to make sure it is a list"))
  }
}


#' Convert orthologous relationships to edgelist
#'
#'Generate a coupling matrix that connects species based on the orthologus relationships between species. Relationships can range from one-to-one to many-to-many
#' @usage getOrthoWeights(ortho, GeneMeta, couple_const=1)
#' @param ortho List object that contains the orthologous gene relationships between species
#' @param GeneMeta Data frame that contains the project metadata. See \link{createGeneMeta}
#' @param couple_const A contant that determines the relative contribution of orthologous gene relationships between species. Default=1
#' @keywords gene ortholog
#' @references Koon-Kiu Yan, Daifeng Wang, Joel Rozowsky, Henry Zheng, Chao Cheng and Mark Gerstein. 2014. OrthoClust: an orthology-based network framework for clustering data across multiple species. Genome Biology. 15:R100
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
getOrthoWeights<-function(ortho,GeneMeta,couple_const=1)
{
  gnames<-as.data.frame(GeneMeta)
  row.names(gnames)<-gnames[,2]
  gnames[,3]<-1:nrow(gnames)

  if(is.list(ortho) & length(ortho)>0)
  {
    ortho_weights<-data.frame(matrix(nrow=0,ncol=3))

    for( k in 1:length(ortho))
    {
      index<-intersect(which(ortho[[k]][,1] %in% gnames[,2]), which(ortho[[k]][,2] %in% gnames[,2]))
      tmp<-ortho[[k]][index,]

      geneA<-as.character(tmp[,1])
      geneB<-as.character(tmp[,2])
      TgeneA<-table(geneA)
      TgeneB<-table(geneB)
      ow1<-TgeneA[tmp[,1]]
      ow2<-TgeneB[tmp[,2]]
      w<-as.data.frame((1/ow1+1/ow2)/2)[,2]*couple_const
      o<-data.frame(cbind(gnames[tmp[,1],3],gnames[tmp[,2],3],w))
      ortho_weights<-rbind(ortho_weights,o)
    }

    Rortho_weights<-ortho_weights[,c(2,1,3)]
    names(Rortho_weights)<-names(ortho_weights)[1:3]
    ortho_sym<-rbind(ortho_weights[,1:3],Rortho_weights)
    names(ortho_sym)<-c("variable","value","weight")
    return(ortho_sym[order(ortho_sym[,1],ortho_sym[,2]),])
  } else {
    return(print("Check rpkm to make sure it is a list"))
  }
}

#' Louvain community detection
#'
#'Run the Louvain clustering method multiple times on a user defined set of edges. For each run, the edgelist is shuffled so the louvain method begins at a random starting point in the network.
#' @usage louvain(edgelist, nruns)
#' @param edgelist A data frame that contains three columns. The first two columns define the edges between gene_A and gene_B. Gene names must be converted to integer values. The third column defines the numeric edge weight.
#' @param nruns An interger defining the number of runs.
#' @keywords louvain
#' @references Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre. 2008. Fast unfolding of communities in large networks. J. Stat. Mech. P10008
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @seealso \code{\link[igraph]{cluster_louvain}}
#' @export
louvain<-function(edgelist,nruns)
{
  if(is.data.frame(edgelist) & ncol(edgelist)==3)
  {
    colnames(edgelist)<-c("V1","V2","weight")

    for(h in 1:nruns)
    {
      if(h==1)
      {
        g<-graph_from_data_frame(edgelist,directed = F)
        community<-cluster_louvain(g)
        memb<-membership(community)
        results<-as.data.frame(as.numeric(names(memb)))
        results[,h+1]<-memb
        print(modularity(community))
      } else {
        rand_edge<-edgelist[order(rnorm(nrow(edgelist))),]
        Rg<-graph_from_data_frame(rand_edge,directed = F)
        Rcommunity<-cluster_louvain(Rg)
        Rmemb<-membership(Rcommunity)

        results[,h+1]<-Rmemb[order(as.numeric(names(Rmemb)))]
        print(modularity(Rcommunity))
      }

    }
    return(results)
  } else {
    retuen("Edges are not a data frame or does not have 3 columns")
  }
}

#Generate weighted edgelist using adjacency and low end cutoff
getEdgelistWeighted<-function(rpkm,nGenes,power=c(6),threshold=0.8,nThreads = 3)
{
  if(is.list(rpkm) & length(rpkm)>0 & length(rpkm)==length(power))
  {
    full_edgelist<-data.frame(matrix(nrow = 0,ncol = 3))
    for(q in 1:length(rpkm))
    {

      datExpr0<-t(rpkm[[q]])
      colnames(datExpr0)<-1:ncol(datExpr0)
      Cor<-corFast(datExpr0,nThreads = nThreads)
      diag(Cor)<-0
      Cor<-Cor^power[q]
      collectGarbage()
      Cor[Cor<threshold]<-NA
      edgelist<-melt(Cor,na.rm = T)

      if(q==1)
      {
        edgelistB<-edgelist
        edgelistB[,1]<-edgelist[,1]
        edgelistB[,2]<-edgelist[,2]
      } else {
        edgelistB<-edgelist
        edgelistB[,1]<-edgelist[,1]+nGenes[q-1]
        edgelistB[,2]<-edgelist[,2]+nGenes[q-1]
      }

      collectGarbage()

      full_edgelist<-rbind(full_edgelist,edgelistB)
      #write.table(edgelist,file=paste("~/Documents/scripts/louvain-generic/PSE_R_var/",names(rpkm)[q] ,"_edgelist.txt",sep=""),sep="\t",row.names = F,col.names = F,quote = F)

      rm(Cor)
      collectGarbage()
    }

    names(full_edgelist)<-c("gene1","gene2","weight")
    return(full_edgelist[order(full_edgelist[,1],full_edgelist[,2]),])
  } else {
    return(print("Check rpkm to make sure it is a list"))
  }
}

#' Plot MultiSpp co-appearance matrix
#'
#' Plot co-appearance heatmap of co-expression network generated across multiple species
#' @param GeneMeta Data frame that contains the project metadata. See \link{createGeneMeta}
#' @param order The order in which genes should be plotted.
#' @param CA_keep Sparse matrix of louvain community assignments generated from many runs of \code{\link{louvain}}
#' @param sb An integer specifying the increment of the sequence to plot. Example plot every 12th gene in the order.
#' @param remove_0 A logical value indicateing if the unclusterable genes should be plotted.
#' @param text_rotate Angle of module labels.
#' @param lwd Line weight
#' @param cex Font proportion
#' @param my_palette Specify your own color ramp for the heatmap
#' @import methods
#' @importFrom Matrix tcrossprod
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @seealso \code{\link[igraph]{cluster_louvain}}
#' @export
plot_MultiSpp<-function(GeneMeta,order,CA_keep,sb=12, remove_0=TRUE, text_rotate=NULL, lwd=0.5, cex=0.2, my_palette=colorRampPalette(c("white", "lightyellow", "red","black"))(n = 100))
{
  #reorder the data frame
  gene_order<-GeneMeta[order,]

  #remove the unclusterable genes in the _0 or grey group
  if(remove_0==TRUE)
  {
    gene_order<-gene_order[grep("_0",gene_order[,3],invert = T),]
  }

  #subset of the data
  i<-seq(1,length(gene_order$ID),sb)

  IDs<-gene_order$ID[i]
  Wmatrix<-as.matrix(tcrossprod(CA_keep[IDs,]))

  uMods<-unique(gene_order[i,3])
  uMods<-uMods[grep("_0",uMods,invert = T)]

  #generate the coordinates of where the modules located in the Wmatrix
  coord<-data.frame(matrix(NA,0,2))
  for(u in uMods)
  {
    coord<-rbind(coord,range(which(gene_order[i,3]==u)))
  }
  coordN<-coord/length(i)

  #generate the coordinates for each species
  nG<-table(gene_order[i,1])[unique(gene_order[i,1])]
  nm<-cumsum(nG)/sum(nG)
  #nm<-nm[-length(nm)]
  nm<-c(0,nm)
  tmp<-data.frame(matrix(NA,0,4))
  #vertical
  for(j in 1:(length(nm)-1))
  {
    tmp[j,]<-c(nm[j],0,nm[j],nm[j+1])
  }
  #horizontal
  for(j in 1:(length(nm)-1))
  {
    tmp[j+(length(nm)-1),]<-c(nm[j],nm[j+1],1,nm[j+1])
  }

  #set lower tri to NA
  nc<-ncol(Wmatrix)
  ns<-c(0,cumsum(nG))
  for(s in 1:(length(ns)-1))
  {
    Wmatrix[ns[s]:ns[s+1],ns[s+1]:nc]<-NA
  }



  #pdf(file=paste(output_dir,"Co_appearance_louvain_lower.pdf",sep=""),w=10,h=10)
  par(mar=c(5,4,2,9))
  image(Wmatrix, col=my_palette, useRaster=TRUE, axes=FALSE)
  rect(coordN[,1], coordN[,1], coordN[,2], coordN[,2],border = "black",lwd = lwd)
  #Add line segments
  segments(tmp[,1],tmp[,2],tmp[,3],tmp[,4],lwd=lwd)
  abline(h=0,lwd=lwd); abline(v=1,lwd=lwd);
  text(coordN[,1]-0.1,rowMeans(coordN),uMods,srt = text_rotate, cex=cex)
  #dev.off()

}


#calculate kME ie correlation of expression to eigengene
kME<-function(gene_merged,rpkm,MEs)
{
  kME=NULL
  for(sp in unique(gene_merged[,1]))
  {
    tmp_rpkm<-rpkm[[sp]]
    n=ncol(tmp_rpkm)
    tmp_ME<-MEs[[sp]]
    tmp_genes<-gene_merged[which(gene_merged[,1]==sp),]
    tmp<-cbind(tmp_rpkm[tmp_genes[,2],],t(tmp_ME[,tmp_genes[,3]]))
    out<-apply(tmp,1,function(x) cor(x[1:(length(x)/2)],x[((length(x)/2)+1):length(x)]))
    kME<-c(kME,out)

  }
  return(cbind(gene_merged,kME))
}

getEdgelist_from_GeneNames<-function(edges,GeneMeta)
{
  tmp_names<-cbind(GeneMeta,seq(nrow(GeneMeta)))
  row.names(tmp_names)<-GeneMeta[,2]
  names(tmp_names)<-c("spp","gene_names","index")

  if(is.list(edges) & length(edges)>0)
  {
    tmp_edgelist<-data.frame(matrix(nrow = 0,ncol = 3))
    for(q in 1:length(edges))
    {
      tmp<-edges[[q]][[1]]
      edgelist<-data.frame(matrix(nrow = length(tmp),ncol = 3))
      edgelist[,1]<-tmp_names[edges[[q]][[1]],3]
      edgelist[,2]<-tmp_names[edges[[q]][[2]],3]
      edgelist[,3]<-1
      tmp_edgelist<-rbind(tmp_edgelist,edgelist)
    }
  }

  Rtmp_edgelist<-tmp_edgelist[,c(2,1,3)]
  names(Rtmp_edgelist)<-names(tmp_edgelist)
  fullE<-rbind(tmp_edgelist,Rtmp_edgelist)
  names(fullE)<-c("variable","value","weight")
  return(fullE[order(fullE[,1],fullE[,2]),])
}

#'Plot heatmap color scale
#'
#'Plot a user defined color scale
#' @usage color.bar(lut, min, max=-min, nticks=11, title='')
#' @param lut User specified color ramp
#' @param min Minimum value
#' @param max Maximum value
#' @param nticks Number of ticks
#' @param title Title for color scale.
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
#'
color.bar <- function(lut, min, max=-min, nticks=11, title='') {
  ticks=seq(min, max, len=nticks)
  scale = (length(lut)-1)/(max-min)

  plot(c(min,max),c(0,2),  type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(1, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    #rect(0,y,2,y+1/scale, col=lut[i], border=NA)
    rect(ybottom = 0,xleft = y,ytop = 2,xright = y+1/scale, col=lut[i], border=NA)
  }
}

#'Hierarchical clustering of gene co-appearance in Louvain communities
#'
#'For each species calculate the a co-appeance matrix representing how often genes are assigned to the same Louvain communities and perform hierarchical clustering to determine which genes have similar co-appearance.
#' @usage multiSppHclust(occurance, nRuns, GeneMeta)
#' @param occurance Sparse matrix that contain the presence (1) and absence (0) of each gene in each Louvain community.
#' @param nRuns Integer specifying the number of Louvain runs.
#' @param GeneMeta Data frame that contains the project metadata. See \link{createGeneMeta}
#' @import methods
#' @importFrom Matrix tcrossprod
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
multiSppHclust<-function(occurance, nRuns, GeneMeta)
{
  order=NULL
  MultiSpp_trees<-list()
  for(s in unique(GeneMeta[,1]))
  {
    Sgenes<-which(GeneMeta[,1]==s)
    Smatrix<-tcrossprod(occurance[Sgenes,])
    Smatrix_sim<-Smatrix/nRuns
    rm(Smatrix); collectGarbage();
    SNmatrix<-1-Smatrix_sim;
    rm(Smatrix_sim); collectGarbage();
    tree<-fastcluster::hclust(as.dist(SNmatrix),method = "average")
    tree$height<-round(tree$height,12)
    MultiSpp_trees[[s]]<-tree
    order<-c(order,GeneMeta$ID[Sgenes[tree$order]])
    rm(SNmatrix); rm(tree); collectGarbage();
    print(paste("Finished tree: ",s,sep=""))
  }
  return(list(order=order,trees=MultiSpp_trees))
}


#'Identify modules with high co-appearance
#'
#'Use cutreeDynamic to identify gene modules that have high co-appearance.
#' @usage multiSppModules(multiSpp_results, GeneMeta, minModuleSize, cut)
#' @param multiSpp_results A list object containing the gene order and hclust dendrograms for each species. See \link{multiSppHclust}
#' @param GeneMeta Data frame that contains the project metadata. See \link{createGeneMeta}
#' @param minModuleSize A vector containing the minimum module size to be used for each species in the data set.
#' @param cut A vector containing the cut height to be used for species in the data set
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
multiSppModules<-function(multiSpp_results, GeneMeta, minModuleSize, cut)
{
  trees<-multiSpp_results$trees
  order<-multiSpp_results$order
  MultiSpp_Mods<-list()
  if(length(minModuleSize)==length(trees) & length(cut)==length(trees))
     {
  for(p in 1:length(unique(GeneMeta[,1])))
  {
    Sgenes<-which(GeneMeta[,1]==unique(GeneMeta[,1])[p])
    spp<-unique(GeneMeta[,1])[p]
    Mods= cutreeDynamic(dendro = trees[[p]], method="tree",
                        deepSplit = 1, pamRespectsDendro = FALSE,
                        minClusterSize = minModuleSize[p],cutHeight=cut[p]);
    MultiSpp_Mods[[spp]] <-paste(unique(GeneMeta[,1])[p],Mods,sep="_")

  }
  } else {
    return("Please make sure there is a minModuleSize and cut specified for each species")
  }
  return(MultiSpp_Mods)
}

#'Parse inParanoid table
#'
#'Convert an inParanoid table format and output a table of ortholog edges
#' @usage parseInParanoid(MetaDataInParanoid, outDir = ".")
#' @param MetaDataInParanoid Provide a data frame where each row provides the meta data for a single inParanoid run. The resulting data frame must have 3 columns: 1) file path to table.inparanoid_output 2) Species A alias and 3) Species B alias
#' @param outDir Directory to write files
#' @return The function will convert the table.inParanoid format to a two column edgelist and print the results to a text file SppA_SppB_orthologs.txt
#' @examples
#' load("Data/inParanoid_meta.rdata")
#' parseInParanoid(Meta_Data)
#' @author Matthew Zinkgraf, \email{mzinkgraf@gmail.com}
#' @export
parseInParanoid<-function(MetaDataInParanoid, outDir = ".")
{
  options(stringsAsFactors = FALSE);
if(ncol(MetaDataInParanoid)!=3) {stop("meta data does not contain 3 columns")}
  for(f in 1:nrow(MetaDataInParanoid))
  {
    if(!file.exists(MetaDataInParanoid[f,1])) {stop(paste(MetaDataInParanoid[f,1],"File does not exist!"))}

    tb<-read.table(MetaDataInParanoid[f,1],sep="\t",header=T)

    #split each species column

    spp1<-strsplit(as.character(tb$OrtoA), split="\\s\\d\\.\\d+\\s", perl=T)

    spp2<-strsplit(as.character(tb$OrtoB), split="\\s\\d\\.\\d+\\s", perl=T)

    n1<-length(spp1)
    n2<-length(spp2)

    output<-data.frame(matrix(nrow = 0,ncol=2))
    names(output)<-c("spp1","spp2")
    if(n1==n2)
      for(i in 1:n1)
      {
        for(j in spp1[[i]])
        {
          for(k in spp2[[i]])
          {
            tmp<-c(j,k)
            names(tmp)<-c("spp1","spp2")
            output<-rbind(output,tmp)
          }
        }
      }
    names(output)<-c("spp1","spp2")
    file.out<-paste(MetaDataInParanoid[f,2],MetaDataInParanoid[f,3],"orthologs.txt",sep="_")
    write.table(output,file=paste(outDir,file.out,sep="/"),sep="\t",col.names = F,row.names = F,quote = F)

  }
}
