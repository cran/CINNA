## ------------------------------------------------------------------------

library(CINNA)


## ----warning=FALSE,message=FALSE-----------------------------------------

data("zachary")
zachary


## ----warning=FALSE,message=FALSE-----------------------------------------

graph.extract.components(zachary)


## ------------------------------------------------------------------------

data(drugtarget)

drug.comp<-graph.extract.components( drugtarget, directed = TRUE,bipartite.proj=TRUE ,num.proj=2)
head(drug.comp)


## ----warning=FALSE,message=FALSE-----------------------------------------
library(igraph)
zachary.edgelist<-as_edgelist(zachary)

misc.extract.components(zachary.edgelist)


## ----warning=FALSE,message=FALSE-----------------------------------------

giant.component.extract(zachary)


## ----warning=FALSE,message=FALSE-----------------------------------------

proper.centralities(zachary)


## ----warning=FALSE,message=FALSE-----------------------------------------

calculate.centralities(zachary)


## ----warning=FALSE,message=FALSE-----------------------------------------

pr.cent<-proper.centralities(zachary)

calc.cent<-calculate.centralities(zachary, except=pr.cent[5:40])


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "A display of most informative centrality measures based on principal component analysis. The red line indicates the random threshold of contribution. This barplot represents contribution of variable values based on the number of dimensions."----

pca.centralities( calc.cent )


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "A representation of most informative centrality measures based on principal component analysis between unscaled(not normalized) centrality values."----

pca.centralities( calc.cent , scale.unit = FALSE )


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "A display of most informative centrality measures based on t-Distributed Stochastic Neighbor Embedding analysis among scaled(not normalized) centrality values."----

tsne.centralities( calc.cent, dims = 2, perplexity = 1, scale=TRUE)


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = " Graph illustration based on centrality measure. The size of nodes represent the degree centrality values."----

visualize.graph( zachary , centrality.type="Degree Centrality")


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "Observed centrality measure heatmap. The colors from blue to red displays scaled centrality values."----

visualize.heatmap( calc.cent , scale = TRUE  )


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "A display of correlation among computed centrality measures.  The red to blue highlighted circles represent the top to bottom Pearson correlation coefficients[@Benesty2009] which differ from -1 to 1. The higher the value becomes larger, circles' sizes get larger too."----

visualize.correlations(calc.cent,"pearson")


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "Circular dendrogram plot of vertices based on specified centrality measure. Each color represents a cluster."----

visualize.dendogram(zachary,k=4)


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "Association plot between two centrality variables. The red line is an indicator of linear regression line among them."----

subgraph_cent<-calc.cent[[1]]
Topological_coef<-calc.cent[[2]]

visualize.association(  subgraph_cent, Topological_coef)


## ----fig.width=7, fig.height=6,message=FALSE, fig.cap = "Pairwise Pearson correlation between two centrality values."----

visualize.pair.correlation( subgraph_cent ,Topological_coef)


