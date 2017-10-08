% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CINNA.R
\name{graph_extract_components}
\alias{graph_extract_components}
\title{Component extraction of a graph}
\usage{
graph_extract_components(x, directed = TRUE, bipartite_proj = FALSE,
  num_proj = 1)
}
\arguments{
\item{x}{An igraph or a network object}

\item{directed}{Whether to create a directed graph(default=TRUE)}

\item{bipartite_proj}{Whether the bipartite network must be projected or not(default=FALSE)}

\item{num_proj}{Numbers 1 or 2 which shows the number of projects for bipartite graphs.(default=1)}
}
\value{
a list including the componets of the input as igraph objects
}
\description{
This function extracts all connected components of the
input which can be an "igraph" object or a "network" object
and convert them as "igraph" objects.
}
\details{
This function seperates different components of an "igraph" or a "network" object and
illustrates them as a list of independent graphs. If the input graph was bipartite and the
"bipartite_proj" was TRUE, it will project it and you can decide in which project you want
to continue to work with.
}
\examples{

data(zachary)

graph_extract_components(zachary)

}
\seealso{
\code{\link[igraph]{induced.subgraph}}, \code{\link[igraph]{components}}
}
\author{
Minoo Ashtiani, Mohieddin Jafari
}