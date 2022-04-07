#' @title Component extraction of a graph
#'
#' @description This function extracts all connected components of the
#' input which can be an "igraph" object or a "network" object
#' and convert them as "igraph" objects.
#' @param x An igraph or a network object
#' @param directed Whether to create a directed graph(default=TRUE)
#' @param bipartite_proj Whether the bipartite network must be projected or not(default=FALSE)
#' @param num_proj Numbers 1 or 2 which shows the number of projects for bipartite graphs.(default=1)
#' @details
#' This function seperates different components of an "igraph" or a "network" object and
#' illustrates them as a list of independent graphs. If the input graph was bipartite and the
#' "bipartite_proj" was TRUE, it will project it and you can decide in which project you want
#' to continue to work with.
#' @seealso \code{\link[igraph]{induced.subgraph}},\code{\link[igraph]{components}}
#' @return a list including the components of the input as igraph objects
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @examples
#'
#' data(zachary)
#'
#' graph_extract_components(zachary)
#'
#'
#' @export
#' @importFrom igraph is_igraph
#' @importFrom igraph is_bipartite
#' @importFrom igraph bipartite.projection
#' @importFrom igraph simplify
#' @importFrom igraph is_simple
#' @importFrom igraph clusters
#' @importFrom igraph induced.subgraph
#' @importFrom igraph graph_from_edgelist
#' @importFrom network is.network
#' @importFrom network as.edgelist
#' @importFrom network network

graph_extract_components <- function( x, directed = TRUE, bipartite_proj=FALSE, num_proj=1){

  if(!("igraph" %in% class(x)|| "network" %in% class(x))) stop("The input is not an igraph or a
                                                           network object")

  if(is_igraph(x)){

    if( bipartite_proj){

      if(is_bipartite(x)){

        x<-bipartite.projection(x)[[num_proj]]

        if (!is_simple(x))   x<-simplify(x)

        cl <- clusters(x)

        graph.splitting <- function(k, x, cl){
          induced.subgraph(x, cl$membership == k)
        }

        components<-sapply(1:max(cl$membership), graph.splitting, x = x, cl = cl, simplify = FALSE)

      }

    }

    else{

      if (!is_simple(x))   x<-simplify(x)

      cl <- clusters(x)

      graph_splitting <- function(k, x, cl){
        induced.subgraph(x, cl$membership == k)
      }

      components<-sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

    }

  }

  if( is.network(x)){

    edgelist<-as.edgelist(x)

    x<-graph_from_edgelist(edgelist, directed = TRUE)

    if (!is_simple(x))  x<-simplify(x)

    cl <- clusters(x)

    graph_splitting <- function(k, x, cl){
      induced.subgraph(x, cl$membership == k)
    }

    components<-sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

  }

  return(components)

}

#' @title Component extraction of miscellaneous graph formats
#'
#' @description This function extracts all components of the input with various formats
#' and convert them as "igraph" objects.
#' @param x The input could be an edgelist and an adjacency matrix
#' @param directed Whether to create a directed graph.(default=TRUE)
#' @param mode Character scalar, explain how should demonstrate the supplied matrix.
#' Possible values are: directed, undirected, upper, lower, max, min, plus.(default="directed")
#' @param weighted An argument for specifying whether the graph should be weighted or not.
#'  If it is NULL then an unweighted graph is created.(default=NULL)
#' @param unibipartite A boolean parameter describing whether the input edge list is corresponding
#' to a bipartite graph. TRUE value specifies the biprtite graph and vice versa.(default=FALSE)
#' @param diag Logical scalar, whether to consider the diagonal of the matrix or not.
#' If it was FALSE then the diagonal spotted as zeros.(default=TRUE)
#' @details
#' This function assert components of the input object which can be an edgelist,
#' an adjacency matrix and a graphNEL object.
#' The result would be a list including components as seperated graphs.
#' @seealso \code{\link[igraph]{induced.subgraph}},\code{\link[igraph]{components}},
#' \code{\link[igraph]{graph_from_adjacency_matrix}}
#'
#' @return a list including the components of the input graph as igraph objects
#'
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#'
#' @export
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph %>%
#' @importFrom igraph simplify
#' @importFrom igraph clusters
#' @importFrom igraph induced.subgraph
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph as_edgelist
#' @importFrom network is.network
#' @importFrom network as.edgelist
#' @importFrom network network
#' @importFrom qdapTools mtabulate
#' @importFrom igraph graph_from_incidence_matrix
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom plyr .

misc_extract_components <- function( x ,directed = TRUE, mode = "directed",
                                     weighted = NULL, unibipartite = FALSE,
                                     diag = TRUE){
  if (ncol(x)%in%2) {

    if (unibipartite%in%FALSE){

      x <- graph_from_edgelist(x, directed = directed)

      if (!is_simple(x))   x<-simplify(x)

      cl <- clusters(x)

      graph_splitting <- function(k, x, cl){

        induced.subgraph(x, cl$membership == k)

      }

      components <- sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

      return(components)

    }
    else{

      el <- cbind(x,1)

      incidence_mat <- el[rep(seq_len(nrow(el)), el[,'1']), c(colnames(el)[1], colnames(el)[2])] %>%
      {split(.[,colnames(el)[2]], .[,colnames(el)[1]])} %>%
        mtabulate()

      x <- graph_from_incidence_matrix(incidence_mat, directed=directed)


      if (!is_simple(x))   x <- simplify(x)

      cl <- clusters(x)

      graph_splitting <- function(k, x, cl){

        induced.subgraph(x, cl$membership == k)

      }

      components <- sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

      return(components)

    }

  }
  if (ncol(x)>2||class(x)%in%"dgCMatrix") {

    x <- graph_from_adjacency_matrix(x, mode = mode, weighted = weighted, diag = diag,
                                     add.colnames = NULL, add.rownames = NA)

    if (!is_simple(x))   x<-simplify(x)

    cl <- clusters(x)

    graph_splitting <- function(k, x, cl){

      induced.subgraph(x, cl$membership == k)

    }

    components <- sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

    return(components)

  }
}

#' @title Giant component extraction of a graph
#'
#' @description This function extracts the largest connected or the giant component
#' of the input graph which can be an "igraph" object or a "network" object
#' and convert them as "igraph" objects. For the bipartite graph, this will
#' apply projection before extracting the components.
#' @param x An igraph or a network object
#' @param directed Whether to create a directed graph(default=TRUE)
#' @param bipartite_proj Whether the bipartite network must be projected or not(default=FALSE)
#' @param num_proj A number which shows the number of projects especifically for
#' bipartite graphs.(default=1)
#'
#' @details
#' This function distinguishes the largest component of an "igraph" or a "network" object and
#' illustrates them as a list which contains the edgelist of the giant component. If the input
#' graph was bipartite and the "bipartite_proj" was TRUE, it will projet it and you can decide
#' to which project you want to continue to work with that.
#' @seealso \code{\link[igraph]{induced.subgraph}},\code{\link[igraph]{clusters}}
#'
#' @return a list contatining the giant componet of the input graph. The first element is an igraph object and the second is the edgelist of that.
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Newman, M. (2010). Networks. Oxford University Press.
#' @examples
#' # a graph with 4 vertices
#'
#' data(zachary)
#' giant_component_extract(zachary)
#'
#' @export
#' @importFrom igraph is_igraph
#' @importFrom igraph is_bipartite
#' @importFrom igraph bipartite.projection
#' @importFrom igraph is_simple
#' @importFrom igraph simplify
#' @importFrom igraph clusters
#' @importFrom igraph induced.subgraph
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph as_edgelist
#' @importFrom network is.network
#' @importFrom network as.edgelist
#' @importFrom network network


giant_component_extract <- function( x, directed = TRUE, bipartite_proj=FALSE ,num_proj=1){

  if (!("igraph" %in% class(x)|| "network" %in% class(x))) stop("The input is not an igraph or a network object")

  if (is_igraph(x)){

    if ( bipartite_proj){

      if (is_bipartite(x)){

        x <- bipartite.projection(x)[[num_proj]]

        if (!is_simple(x))   x <- simplify(x)

        cl <- clusters(x)

        giant_comp <- induced.subgraph(x, which(cl$membership == which.max(cl$csize)))

        giant_comp_edgelist <- as_edgelist(giant_comp, names = TRUE)

      }

      else stop("The graph is not bipartite")

    }

    else{

      if (!is_simple(x))   x<-simplify(x)

      cl <- clusters(x)

      giant_comp <- induced.subgraph(x, which(cl$membership == which.max(cl$csize)))

      giant_comp_edgelist <- as_edgelist(giant_comp, names = TRUE)

    }

  }

  if ( is.network(x)){

    edgelist<-as.edgelist(x)

    x <- graph_from_edgelist(edgelist, directed = directed)

    if (!is_simple(x))  x<-simplify(x)

    cl <- clusters(x)

    giant_comp <- induced.subgraph(x, which(cl$membership == which.max(cl$csize)))

    giant_comp_edgelist <- as_edgelist(giant_comp, names = TRUE)

  }

  result <- list(giant_comp,giant_comp_edgelist)

  return(result)

}


#' @title Proper centrality measure representation
#'
#' @description This function indicates proper centrality measures of an
#' igraph object based on the network topology
#' @param x an igraph object
#' @details
#' This function represents a list including the names of centrality measures which are applicable
#' for the input graph based on the topology
#' @seealso \code{\link[CINNA]{calculate_centralities}}
#' @return a list including the name of centrality measures which are suitable for the input graph
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @examples
#'
#' data("zachary")
#' proper_centralities(zachary)
#'
#' @export
#' @importFrom igraph is.igraph
#' @importFrom igraph is_directed
#' @importFrom igraph is.weighted

proper_centralities<-function(x){

  if (!is.igraph(x)) stop(" Error: x must be a class of igraph object ")

  if (is_directed(x) && is_weighted(x)){

    proper_centralities <- c("Alpha Centrality",
                             "Burt's Constraint",
                             "Page Rank" ,
                             "Average Distance",
                             "Barycenter Centrality" ,
                             "BottleNeck Centrality",
                             "Centroid value" ,
                             "Closeness Centrality (Freeman)" ,
                             "ClusterRank" ,
                             "Decay Centrality",
                             "Degree Centrality" ,
                             "Diffusion Degree",
                             "DMNC - Density of Maximum Neighborhood Component" ,
                             "Eccentricity Centrality" ,
                             "Harary Centrality",
                             "eigenvector centralities" ,
                             "K-core Decomposition" ,
                             "Geodesic K-Path Centrality" ,
                             "Katz Centrality (Katz Status Index)" ,
                             "Kleinberg's authority centrality scores",
                             "Kleinberg's hub centrality scores" ,
                             "clustering coefficient" ,
                             "Lin Centrality" ,
                             "Lobby Index (Centrality)" ,
                             "Markov Centrality" ,
                             "Radiality Centrality",
                             "Shortest-Paths Betweenness Centrality" ,
                             "Current-Flow Closeness Centrality",
                             "Closeness centrality (Latora)" ,
                             "Communicability Betweenness Centrality",
                             "Community Centrality" ,
                             "Cross-Clique Connectivity" ,
                             "Entropy Centrality" ,
                             "EPC - Edge Percolated Component" ,
                             "Laplacian Centrality" ,
                             "Leverage Centrality" ,
                             "MNC - Maximum Neighborhood Component" ,
                             "Hubbell Index" ,
                             "Semi Local Centrality",
                             "Closeness Vitality" ,
                             "Residual Closeness Centrality",
                             "Stress Centrality",
                             "Load Centrality",
                             "Flow Betweenness Centrality",
                             "Information Centrality",
                             "Dangalchev Closeness Centrality",
                             "Group Centrality",
                             "Harmonic Centrality",
                             "Local Bridging Centrality",
                             "Wiener Index Centrality",
                             "Weighted Vertex Degree"

    )
  }

  if (!is_directed(x) && is_weighted(x)) {

    proper_centralities<-c(    "subgraph centrality scores",
                               "Topological Coefficient",
                               "Average Distance",
                               "Barycenter Centrality" ,
                               "BottleNeck Centrality",
                               "Centroid value" ,
                               "Closeness Centrality (Freeman)" ,
                               "ClusterRank" ,
                               "Decay Centrality",
                               "Degree Centrality" ,
                               "Diffusion Degree",
                               "DMNC - Density of Maximum Neighborhood Component" ,
                               "Eccentricity Centrality" ,
                               "Harary Centrality",
                               "eigenvector centralities" ,
                               "K-core Decomposition" ,
                               "Geodesic K-Path Centrality" ,
                               "Katz Centrality (Katz Status Index)" ,
                               "Kleinberg's authority centrality scores",
                               "Kleinberg's hub centrality scores" ,
                               "clustering coefficient" ,
                               "Lin Centrality" ,
                               "Lobby Index (Centrality)" ,
                               "Markov Centrality" ,
                               "Radiality Centrality",
                               "Shortest-Paths Betweenness Centrality" ,
                               "Current-Flow Closeness Centrality",
                               "Closeness centrality (Latora)" ,
                               "Communicability Betweenness Centrality",
                               "Community Centrality" ,
                               "Cross-Clique Connectivity" ,
                               "Entropy Centrality" ,
                               "EPC - Edge Percolated Component" ,
                               "Laplacian Centrality" ,
                               "Leverage Centrality" ,
                               "MNC - Maximum Neighborhood Component" ,
                               "Hubbell Index" ,
                               "Semi Local Centrality",
                               "Closeness Vitality" ,
                               "Residual Closeness Centrality",
                               "Stress Centrality",
                               "Load Centrality",
                               "Flow Betweenness Centrality",
                               "Information Centrality",
                               "Dangalchev Closeness Centrality",
                               "Group Centrality",
                               "Harmonic Centrality",
                               "Local Bridging Centrality",
                               "Wiener Index Centrality",
                               "Weighted Vertex Degree"
    )

  }

  if (!is_directed(x) && !is_weighted(x)) {
    proper_centralities <- c( "subgraph centrality scores",
                              "Topological Coefficient",
                              "Average Distance",
                              "Barycenter Centrality" ,
                              "BottleNeck Centrality",
                              "Centroid value" ,
                              "Closeness Centrality (Freeman)" ,
                              "ClusterRank" ,
                              "Decay Centrality",
                              "Degree Centrality" ,
                              "Diffusion Degree",
                              "DMNC - Density of Maximum Neighborhood Component" ,
                              "Eccentricity Centrality" ,
                              "Harary Centrality",
                              "eigenvector centralities" ,
                              "K-core Decomposition" ,
                              "Geodesic K-Path Centrality" ,
                              "Katz Centrality (Katz Status Index)" ,
                              "Kleinberg's authority centrality scores",
                              "Kleinberg's hub centrality scores" ,
                              "clustering coefficient" ,
                              "Lin Centrality" ,
                              "Lobby Index (Centrality)" ,
                              "Markov Centrality" ,
                              "Radiality Centrality",
                              "Shortest-Paths Betweenness Centrality" ,
                              "Current-Flow Closeness Centrality",
                              "Closeness centrality (Latora)" ,
                              "Communicability Betweenness Centrality",
                              "Community Centrality" ,
                              "Cross-Clique Connectivity" ,
                              "Entropy Centrality" ,
                              "EPC - Edge Percolated Component" ,
                              "Laplacian Centrality" ,
                              "Leverage Centrality" ,
                              "MNC - Maximum Neighborhood Component" ,
                              "Hubbell Index" ,
                              "Semi Local Centrality",
                              "Closeness Vitality" ,
                              "Residual Closeness Centrality",
                              "Stress Centrality",
                              "Load Centrality",
                              "Flow Betweenness Centrality",
                              "Information Centrality",
                              "Dangalchev Closeness Centrality",
                              "Group Centrality",
                              "Harmonic Centrality",
                              "Local Bridging Centrality",
                              "Wiener Index Centrality"
    )
  }

  if (is_directed(x) && !is_weighted(x) ) {
    proper_centralities <- c("Alpha Centrality" ,
                             "Bonacich power centralities of positions" ,
                             "Page Rank" ,
                             "Average Distance",
                             "Barycenter Centrality" ,
                             "BottleNeck Centrality",
                             "Centroid value" ,
                             "Closeness Centrality (Freeman)" ,
                             "ClusterRank" ,
                             "Decay Centrality",
                             "Degree Centrality" ,
                             "Diffusion Degree",
                             "DMNC - Density of Maximum Neighborhood Component" ,
                             "Eccentricity Centrality" ,
                             "Harary Centrality",
                             "eigenvector centralities" ,
                             "K-core Decomposition" ,
                             "Geodesic K-Path Centrality" ,
                             "Katz Centrality (Katz Status Index)" ,
                             "Kleinberg's authority centrality scores",
                             "Kleinberg's hub centrality scores" ,
                             "clustering coefficient" ,
                             "Lin Centrality" ,
                             "Lobby Index (Centrality)" ,
                             "Markov Centrality" ,
                             "Radiality Centrality",
                             "Shortest-Paths Betweenness Centrality" ,
                             "Current-Flow Closeness Centrality",
                             "Closeness centrality (Latora)" ,
                             "Communicability Betweenness Centrality",
                             "Community Centrality" ,
                             "Cross-Clique Connectivity" ,
                             "Entropy Centrality" ,
                             "EPC - Edge Percolated Component" ,
                             "Laplacian Centrality" ,
                             "Leverage Centrality" ,
                             "MNC - Maximum Neighborhood Component" ,
                             "Hubbell Index" ,
                             "Semi Local Centrality",
                             "Closeness Vitality" ,
                             "Residual Closeness Centrality",
                             "Stress Centrality",
                             "Load Centrality",
                             "Flow Betweenness Centrality",
                             "Information Centrality",
                             "Dangalchev Closeness Centrality",
                             "Group Centrality",
                             "Harmonic Centrality",
                             "Local Bridging Centrality",
                             "Wiener Index Centrality"
    )
  }

  print(proper_centralities)

}


#' @title Centrality measure calculation
#'
#' @description This function computes multitude centrality measures of an igraph object.
#' @param x the component of a network as an igraph object
#' @param except A vector containing names of centrality measures which could be omitted from the calculations.
#' @param include A vector including names of centrality measures which should be computed.
#' @param weights A character scalar specifying the edge attribute to use.(default=NULL)
#' @details
#' This function calculates various types of centrality measures which are applicable to the network topology
#' and returns the results as a list. In "except" argument, you can specify centrality measures which is
#' not necessary to calculate.
#' @seealso \code{ \link[igraph]{alpha.centrality}},\code{ \link[igraph]{bonpow}},\code{ \link[igraph]{constraint}},
#' \code{ \link[igraph]{centr_degree}},\code{ \link[igraph]{eccentricity}},\code{ \link[igraph]{eigen_centrality}},
#' \code{ \link[igraph]{coreness}},\code{ \link[igraph]{authority_score}},\code{ \link[igraph]{hub_score}},
#' \code{ \link[igraph]{transitivity}},\code{ \link[igraph]{page_rank}},\code{ \link[igraph]{betweenness}} ,
#' \code{ \link[igraph]{subgraph.centrality}}, \code{ \link[sna]{flowbet}},\code{ \link[sna]{infocent}},
#' \code{ \link[sna]{loadcent}},\code{ \link[sna]{stresscent}},\code{ \link[sna]{graphcent}}, \code{ \link[centiserve]{topocoefficient}},
#' \code{ \link[centiserve]{closeness.currentflow}},\code{ \link[centiserve]{closeness.latora}},
#' \code{ \link[centiserve]{communibet}}, \code{ \link[centiserve]{communitycent}},
#' \code{ \link[centiserve]{crossclique}},\code{ \link[centiserve]{entropy}},
#' \code{ \link[centiserve]{epc}},\code{ \link[centiserve]{laplacian}},\code{ \link[centiserve]{leverage}},
#' \code{ \link[centiserve]{mnc}},\code{ \link[centiserve]{hubbell}}, \code{ \link[centiserve]{semilocal}},
#' \code{ \link[centiserve]{closeness.vitality}},
#' \code{ \link[centiserve]{closeness.residual}},\code{ \link[centiserve]{lobby}},
#' \code{ \link[centiserve]{markovcent}},\code{ \link[centiserve]{radiality}},\code{ \link[centiserve]{lincent}},
#' \code{ \link[centiserve]{geokpath}}, \code{ \link[centiserve]{katzcent}},\code{ \link[centiserve]{diffusion.degree}},
#' \code{ \link[centiserve]{dmnc}},\code{ \link[centiserve]{centroid}},\code{ \link[centiserve]{closeness.freeman}},
#' \code{ \link[centiserve]{clusterrank}},\code{ \link[centiserve]{decay}},
#' \code{ \link[centiserve]{barycenter}},\code{ \link[centiserve]{bottleneck}},\code{ \link[centiserve]{averagedis}},
#' \code{ \link[CINNA]{local_bridging_centrality}},\code{ \link[CINNA]{wiener_index_centrality}},\code{ \link[CINNA]{group_centrality}},
#' \code{ \link[CINNA]{dangalchev_closeness_centrality}},\code{ \link[CINNA]{harmonic_centrality}},\code{ \link[igraph]{strength}}
#' @return A list concluding centrality measure values in which the columns indicate centralities
#' and the rows show the vertices.
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Bonacich, P., & Lloyd, P. (2001). Eigenvector like measures of centrality for asymmetric relations. Social Networks, 23(3), 191–201.
#'
#' Bonacich, P. (1972). Factoring and weighting approaches to status scores and clique identification. The Journal of Mathematical Sociology, 2(1), 113–120.
#'
#' Bonacich, P. (1987). Power and Centrality: A Family of Measures. American Journal of Sociology, 92(5), 1170–1182.
#'
#' Burt, R. S. (2004). Structural Holes and Good Ideas. American Journal of Sociology, 110(2), 349–399.
#'
#' Batagelj, V., & Zaversnik, M. (2003). An O(m) Algorithm for Cores Decomposition of Networks, 1–9. Retrieved from
#'
#' Seidman, S. B. (1983). Network structure and minimum degree. Social Networks, 5(3), 269–287.
#'
#' Kleinberg, J. M. (1999). Authoritative sources in a hyperlinked environment. Journal of the ACM, 46(5), 604–632.
#'
#' Wasserman, S., & Faust, K. (1994). Social network analysis : methods and applications. American Ethnologist (Vol. 24).
#' Barrat, A., Barthélemy, M., Pastor Satorras, R., & Vespignani, A. (2004). The architecture of complex weighted networks. Proceedings of the National Academy of Sciences of the United States of America , 101(11), 3747–3752.
#'
#' Brin, S., & Page, L. (2010). The Anatomy of a Large Scale Hypertextual Web Search Engine The Anatomy of a Search Engine. Search, 30(June 2000), 1–7.
#'
#' Freeman, L. C. (1978). Centrality in social networks conceptual clarification. Social Networks, 1(3), 215–239.
#'
#' Brandes, U. (2001). A faster algorithm for betweenness centrality*. The Journal of Mathematical Sociology, 25(2), 163–177.
#'
#' Estrada E., Rodriguez-Velazquez J. A.: Subgraph centrality in Complex Networks. Physical Review E 71, 056103.
#'
#' Freeman, L. C., Borgatti, S. P., & White, D. R. (1991). Centrality in valued graphs: A measure of betweenness based on network flow. Social Networks, 13(2), 141–154.
#'
#' Brandes, U., & Erlebach, T. (Eds.). (2005). Network Analysis (Vol. 3418). Berlin, Heidelberg: Springer Berlin Heidelberg.
#'
#' Stephenson, K., & Zelen, M. (1989). Rethinking centrality: Methods and examples. Social Networks, 11(1), 1–37.
#'
#' Wasserman, S., & Faust, K. (1994). Social network analysis : methods and applications. American Ethnologist (Vol. 24).
#'
#' Brandes, U. (2008). On variants of shortest path betweenness centrality and their generic computation. Social Networks, 30(2), 136–145.
#'
#' Goh, K.-I., Kahng, B., & Kim, D. (2001). Universal Behavior of Load Distribution in Scale Free Networks. Physical Review Letters, 87(27), 278701.
#'
#' Shimbel, A. (1953). Structural parameters of communication networks. The Bulletin of Mathematical Biophysics, 15(4), 501–507.
#'
#' Assenov, Y., Ramrez, F., Schelhorn, S.-E., Lengauer, T., & Albrecht, M. (2008). Computing topological parameters of biological networks. Bioinformatics, 24(2), 282–284.
#'
#' Diekert, V., & Durand, B. (Eds.). (2005). STACS 2005 (Vol. 3404). Berlin, Heidelberg: Springer Berlin Heidelberg.
#'
#' Gräßler, J., Koschützki, D., & Schreiber, F. (2012). CentiLib: comprehensive analysis and exploration of network centralities. Bioinformatics (Oxford, England), 28(8), 1178–9.
#'
#' Latora, V., & Marchiori, M. (2001). Efficient Behavior of Small World Networks. Physical Review Letters, 87(19), 198701.
#'
#' Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32(3), 245–251.
#'
#' Estrada, E., Higham, D. J., & Hatano, N. (2009). Communicability betweenness in complex networks. Physica A: Statistical Mechanics and Its Applications, 388(5), 764–774.
#'
#' Hagberg, Aric, Pieter Swart, and Daniel S Chult. Exploring network structure, dynamics, and function using NetworkX. No. LA-UR-08-05495; LA-UR-08-5495. Los Alamos National Laboratory (LANL), 2008.
#'
#' Kalinka, A. T., & Tomancak, P. (2011). linkcomm: an R package for the generation, visualization, and analysis of link communities in networks of arbitrary size and type. Bioinformatics, 27(14), 2011–2012.
#'
#' Faghani, M. R., & Nguyen, U. T. (2013). A Study of XSS Worm Propagation and Detection Mechanisms in Online Social Networks. IEEE Transactions on Information Forensics and Security, 8(11), 1815–1826.
#'
#' Brandes, U., & Erlebach, T. (Eds.). (2005). Network Analysis (Vol. 3418). Berlin, Heidelberg: Springer Berlin Heidelberg.
#'
#' Lin, C.-Y., Chin, C.-H., Wu, H.-H., Chen, S.-H., Ho, C.-W., & Ko, M.-T. (2008). Hubba: hub objects analyzer--a framework of interactome hubs identification for network biology. Nucleic Acids Research, 36(Web Server), W438–W443.
#'
#' Chin, C., Chen, S., & Wu, H. (2009). cyto Hubba: A Cytoscape Plug in for Hub Object Analysis in Network Biology. Genome Informatics …, 5(Java 5), 2–3.
#'
#' Qi, X., Fuller, E., Wu, Q., Wu, Y., & Zhang, C.-Q. (2012). Laplacian centrality: A new centrality measure for weighted networks. Information Sciences, 194, 240–253.
#'
#' Joyce, K. E., Laurienti, P. J., Burdette, J. H., & Hayasaka, S. (2010). A New Measure of Centrality for Brain Networks. PLoS ONE, 5(8), e12200.
#'
#' Lin, C.-Y., Chin, C.-H., Wu, H.-H., Chen, S.-H., Ho, C.-W., & Ko, M.-T. (2008). Hubba: hub objects analyzer--a framework of interactome hubs identification for network biology. Nucleic Acids Research, 36(Web Server), W438–W443.
#'
#' Hubbell, C. H. (1965). An Input Output Approach to Clique Identification. Sociometry, 28(4), 377.
#'
#' Dangalchev, C. (2006). Residual closeness in networks. Physica A: Statistical Mechanics and Its Applications, 365(2), 556–564.
#'
#' Brandes, U. & Erlebach, T. (2005). Network Analysis: Methodological Foundations, U.S. Government Printing Office.
#'
#' Korn, A., Schubert, A., & Telcs, A. (2009). Lobby index in networks. Physica A: Statistical Mechanics and Its Applications, 388(11), 2221–2226.
#'
#' White, S., & Smyth, P. (2003). Algorithms for estimating relative importance in networks. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining   KDD ’03 (p. 266). New York, New York, USA: ACM Press.
#'
#' Cornish, A. J., & Markowetz, F. (2014). SANTA: Quantifying the Functional Content of Molecular Networks. PLoS Computational Biology, 10(9), e1003808.
#'
#' Scardoni, G., Petterlini, M., & Laudanna, C. (2009). Analyzing biological network parameters with CentiScaPe. Bioinformatics, 25(21), 2857–2859.
#'
#' Lin, N. (1976). Foundations of Social Research. Mcgraw Hill.
#'
#' Borgatti, S. P., & Everett, M. G. (2006). A Graph theoretic perspective on centrality. Social Networks, 28(4), 466–484.
#'
#' Newman, M. (2010). Networks. Oxford University Press.
#'
#' Junker, Bjorn H., Dirk Koschutzki, and Falk Schreiber(2006). "Exploration of biological network centralities with CentiBiN." BMC bioinformatics 7.1 : 219.
#'
#' Pal, S. K., Kundu, S., & Murthy, C. A. (2014). Centrality measures, upper bound, and influence maximization in large scale directed social networks. Fundamenta Informaticae, 130(3), 317–342.
#'
#' Lin, C.-Y., Chin, C.-H., Wu, H.-H., Chen, S.-H., Ho, C.-W., & Ko, M.-T. (2008). Hubba: hub objects analyzer--a framework of interactome hubs identification for network biology. Nucleic Acids Research, 36(Web Server), W438–W443.
#'
#' Scardoni, G., Petterlini, M., & Laudanna, C. (2009). Analyzing biological network parameters with CentiScaPe. Bioinformatics, 25(21), 2857–2859.
#'
#' Freeman, L. C. (1978). Centrality in social networks conceptual clarification. Social Networks, 1(3), 215–239.
#'
#' Chen, D.-B., Gao, H., L?, L., & Zhou, T. (2013). Identifying Influential Nodes in Large Scale Directed Networks: The Role of Clustering. PLoS ONE, 8(10), e77455.
#'
#' Jana Hurajova, S. G. and T. M. (2014). Decay Centrality. In 15th Conference of Kosice Mathematicians. Herlany.
#'
#' Viswanath, M. (2009). ONTOLOGY BASED AUTOMATIC TEXT SUMMARIZATION. Vishweshwaraiah Institute of Technology.
#'
#' Przulj, N., Wigle, D. A., & Jurisica, I. (2004). Functional topology in a network of protein interactions. Bioinformatics, 20(3), 340–348.
#'
#' del Rio, G., Koschtzki, D., & Coello, G. (2009). How to identify essential genes from molecular networks BMC Systems Biology, 3(1), 102.
#'
#' Scardoni, G. and Carlo Laudanna, C.B.M.C., 2011. Network centralities for Cytoscape. University of Verona.
#'
#' BOLDI, P. & VIGNA, S. 2014. Axioms for centrality. Internet Mathematics, 00-00.
#'
#' MARCHIORI, M. & LATORA, V. 2000. Harmony in the small-world. Physica A: Statistical Mechanics and its Applications, 285, 539-546.
#'
#' OPSAHL, T., AGNEESSENS, F. & SKVORETZ, J. 2010. Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32, 245-251.
#'
#' OPSAHL, T. 2010. Closeness centrality in networks with disconnected components (http://toreopsahl.com/2010/03/20/closeness-centrality-in-networks-with-disconnected-components/)
#'
#' Michalak, T.P., Aadithya, K.V., Szczepanski, P.L., Ravindran, B. and Jennings, N.R., 2013. Efficient computation of the Shapley value for game-theoretic network centrality. Journal of Artificial Intelligence Research, 46, pp.607-650.
#'
#' Macker, J.P., 2016, November. An improved local bridging centrality model for distributed network analytics. In Military Communications Conference, MILCOM 2016-2016 IEEE (pp. 600-605). IEEE. DOI: 10.1109/MILCOM.2016.7795393
#'
#' DANGALCHEV, C. 2006. Residual closeness in networks. Physica A: Statistical Mechanics and its Applications, 365, 556-564. DOI: 10.1016/j.physa.2005.12.020
#'
#' Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro Vespignani: The architecture of complex weighted networks, Proc. Natl. Acad. Sci. USA 101, 3747 (2004)
#'
#' @examples
#'
#' data("zachary")
#' p <- proper_centralities(zachary)
#' calculate_centralities(zachary, include = "Degree Centrality")
#'
#' @export
#' @importFrom igraph is_directed
#' @importFrom igraph is_connected
#' @importFrom igraph is_weighted
#' @importFrom igraph as_edgelist
#' @importFrom stats setNames
#' @importFrom igraph alpha.centrality
#' @importFrom igraph bonpow
#' @importFrom igraph constraint
#' @importFrom igraph centr_degree
#' @importFrom igraph eccentricity
#' @importFrom igraph eigen_centrality
#' @importFrom igraph coreness
#' @importFrom igraph authority_score
#' @importFrom igraph hub_score
#' @importFrom igraph transitivity
#' @importFrom igraph page_rank
#' @importFrom igraph betweenness
#' @importFrom igraph subgraph.centrality
#' @importFrom igraph strength
#' @importFrom network network
#' @importFrom sna flowbet
#' @importFrom sna infocent
#' @importFrom sna loadcent
#' @importFrom sna stresscent
#' @importFrom sna graphcent
#' @importFrom centiserve topocoefficient
#' @importFrom centiserve closeness.currentflow
#' @importFrom centiserve closeness.latora
#' @importFrom centiserve communibet
#' @importFrom centiserve communitycent
#' @importFrom centiserve crossclique
#' @importFrom centiserve entropy
#' @importFrom centiserve epc
#' @importFrom centiserve laplacian
#' @importFrom centiserve leverage
#' @importFrom centiserve mnc
#' @importFrom centiserve hubbell
#' @importFrom centiserve semilocal
#' @importFrom centiserve closeness.vitality
#' @importFrom centiserve closeness.residual
#' @importFrom centiserve lobby
#' @importFrom centiserve markovcent
#' @importFrom centiserve radiality
#' @importFrom centiserve lincent
#' @importFrom centiserve geokpath
#' @importFrom centiserve katzcent
#' @importFrom centiserve diffusion.degree
#' @importFrom centiserve dmnc
#' @importFrom centiserve centroid
#' @importFrom centiserve closeness.freeman
#' @importFrom centiserve clusterrank
#' @importFrom centiserve decay
#' @importFrom centiserve barycenter
#' @importFrom centiserve bottleneck
#' @importFrom centiserve averagedis

calculate_centralities <- function( x, except = NULL, include = NULL, weights = NULL){

  if (!("igraph" %in% class(x) && is_connected(x))) stop("The input is not an igraph object
                                                       or may not be connected.")
  y <- as_edgelist(x)
  y <- network(y)

  if (is_directed(x) && is_weighted(x) ){

    centrality_funcs <- list(
      "Alpha Centrality" = function(x) alpha.centrality(x, weights = weights),
      "Burt's Constraint" = function(x)constraint(x, weights = weights),
      "Page Rank" = function(x)page_rank(x)$vector,
      "Average Distance" = function(x)averagedis(x, weights = weights),
      "Barycenter Centrality" = function(x)barycenter(x, weights = weights),
      "BottleNeck Centrality" = function(x)bottleneck(x),
      "Centroid value" = function(x)centroid(x, weights = weights),
      "Closeness Centrality (Freeman)" = function(x)closeness.freeman(x, weights = weights),
      "ClusterRank" = function(x)clusterrank(x),
      "Decay Centrality" = function(x)decay(x, weights = weights),
      "Degree Centrality" = function(x)centr_degree(x)$res,
      "Diffusion Degree" = function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component" = function(x)dmnc(x),
      "Eccentricity Centrality" = function(x)eccentricity(x),
      "eigenvector centralities" = function(x)eigen_centrality(x, weights = weights)$vector,
      "K-core Decomposition" = function(x)coreness(x),
      "Geodesic K-Path Centrality" = function(x)geokpath(x, weights = weights),
      "Katz Centrality (Katz Status Index)" = function(x)katzcent(x),
      "Kleinberg's authority centrality scores" = function(x)authority_score(x, weights = weights)$vector,
      "Kleinberg's hub centrality scores" = function(x)hub_score(x, weights = weights)$vector,
      "clustering coefficient" = function(x)transitivity(x, weights = weights,type="local"),
      "Lin Centrality" = function(x)lincent(x, weights = weights),
      "Lobby Index (Centrality)" = function(x)lobby(x),
      "Markov Centrality" = function(x)markovcent(x),
      "Radiality Centrality" = function(x)radiality(x, weights = weights),
      "Shortest-Paths Betweenness Centrality" = function(x)betweenness(x),
      "Current-Flow Closeness Centrality" = function(x)closeness.currentflow(x, weights = weights),
      "Closeness centrality (Latora)" = function(x)closeness.latora(x, weights = weights),
      "Communicability Betweenness Centrality" = function(x)communibet(x),
      "Community Centrality" = function(x)communitycent(x),
      "Cross-Clique Connectivity" = function(x)crossclique(x),
      "Entropy Centrality" = function(x)entropy(x, weights = weights),
      "EPC - Edge Percolated Component" = function(x)epc(x),
      "Laplacian Centrality" = function(x)laplacian(x),
      "Leverage Centrality" = function(x)leverage(x),
      "MNC - Maximum Neighborhood Component" = function(x)mnc(x),
      "Hubbell Index" = function(x)hubbell(x, weights = weights),
      "Semi Local Centrality" = function(x)semilocal(x),
      "Closeness Vitality" = function(x)closeness.vitality(x, weights = weights),
      "Residual Closeness Centrality" = function(x)closeness.residual(x, weights = weights),
      "Stress Centrality" = function(x)stresscent(y),
      "Load Centrality" = function(x)loadcent(y),
      "Flow Betweenness Centrality" = function(x)flowbet(y),
      "Information Centrality" = function(x)infocent(y),
      "Weighted Vertex Degree" = function(x)strength(x, vids = V(x), mode ="all", weights = weights),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = weights),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = weights),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = weights)
    )

    if (!is.null(include)){

      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), include)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
    else{
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), except)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
  }

  if (!is_directed(x) && is_weighted(x)){

    centrality_funcs <- list(
      "subgraph centrality scores" = function(x)subgraph.centrality(x),
      "Topological Coefficient" = function(x)topocoefficient(x),
      "Average Distance" = function(x)averagedis(x, weights = weights),
      "Barycenter Centrality" = function(x)barycenter(x, weights = weights),
      "BottleNeck Centrality" = function(x)bottleneck(x),
      "Centroid value" = function(x)centroid(x, weights = weights),
      "Closeness Centrality (Freeman)" = function(x)closeness.freeman(x, weights = weights),
      "ClusterRank" = function(x)clusterrank(x),
      "Decay Centrality" = function(x)decay(x, weights = weights),
      "Degree Centrality" = function(x)centr_degree(x)$res,
      "Diffusion Degree" = function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component" = function(x)dmnc(x),
      "Eccentricity Centrality" = function(x)eccentricity(x),
      "eigenvector centralities" = function(x)eigen_centrality(x, weights = weights)$vector,
      "K-core Decomposition" = function(x)coreness(x),
      "Geodesic K-Path Centrality" = function(x)geokpath(x, weights = weights),
      "Katz Centrality (Katz Status Index)" = function(x)katzcent(x),
      "Kleinberg's authority centrality scores" = function(x)authority_score(x, weights = weights)$vector,
      "Kleinberg's hub centrality scores" = function(x)hub_score(x, weights = weights)$vector,
      "clustering coefficient" = function(x)transitivity(x, weights = weights,type="local"),
      "Lin Centrality" = function(x)lincent(x, weights = weights),
      "Lobby Index (Centrality)" = function(x)lobby(x),
      "Markov Centrality" = function(x)markovcent(x),
      "Radiality Centrality" = function(x)radiality(x, weights = weights),
      "Shortest-Paths Betweenness Centrality" = function(x)betweenness(x),
      "Current-Flow Closeness Centrality" = function(x)closeness.currentflow(x, weights = weights),
      "Closeness centrality (Latora)" = function(x)closeness.latora(x, weights = weights),
      "Communicability Betweenness Centrality" = function(x)communibet(x),
      "Community Centrality" = function(x)communitycent(x),
      "Cross-Clique Connectivity" = function(x)crossclique(x),
      "Entropy Centrality" = function(x)entropy(x, weights = weights),
      "EPC - Edge Percolated Component" = function(x)epc(x),
      "Laplacian Centrality" = function(x)laplacian(x),
      "Leverage Centrality" = function(x)leverage(x),
      "MNC - Maximum Neighborhood Component" = function(x)mnc(x),
      "Hubbell Index" = function(x)hubbell(x, weights = weights),
      "Semi Local Centrality" = function(x)semilocal(x),
      "Closeness Vitality" = function(x)closeness.vitality(x, weights = weights),
      "Residual Closeness Centrality" = function(x)closeness.residual(x, weights = weights),
      "Stress Centrality" = function(x)stresscent(y),
      "Load Centrality" = function(x)loadcent(y),
      "Flow Betweenness Centrality" = function(x)flowbet(y),
      "Information Centrality" = function(x)infocent(y),
      "Weighted Vertex Degree" = function(x)strength(x, vids = V(x), mode ="all", weights = weights),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="undirected"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = weights),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = weights),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = weights)
    )

    if (!is.null(include)){

      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), include)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
    else{
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), except)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
  }

  if(!is_directed(x) && !is_weighted(x)) {

    centrality_funcs <- list(
      "subgraph centrality scores" = function(x)subgraph.centrality(x),
      "Topological Coefficient" = function(x)topocoefficient(x),
      "Average Distance" = function(x)averagedis(x, weights = NULL),
      "Barycenter Centrality" = function(x)barycenter(x, weights = NULL),
      "BottleNeck Centrality" = function(x)bottleneck(x),
      "Centroid value" = function(x)centroid(x, weights = NULL),
      "Closeness Centrality (Freeman)" = function(x)closeness.freeman(x, weights = NULL),
      "ClusterRank" = function(x)clusterrank(x),
      "Decay Centrality" = function(x)decay(x, weights = NULL),
      "Degree Centrality" = function(x)centr_degree(x)$res,
      "Diffusion Degree" = function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component" = function(x)dmnc(x),
      "Eccentricity Centrality" = function(x)eccentricity(x),
      "eigenvector centralities" = function(x)eigen_centrality(x, weights = NULL)$vector,
      "K-core Decomposition" = function(x)coreness(x),
      "Geodesic K-Path Centrality" = function(x)geokpath(x, weights = NULL),
      "Katz Centrality (Katz Status Index)" = function(x)katzcent(x),
      "Kleinberg's authority centrality scores" = function(x)authority_score(x, weights = NULL)$vector,
      "Kleinberg's hub centrality scores" = function(x)hub_score(x, weights = NULL)$vector,
      "clustering coefficient" = function(x)transitivity(x, weights = NULL,type="local"),
      "Lin Centrality" = function(x)lincent(x, weights = NULL),
      "Lobby Index (Centrality)" = function(x)lobby(x),
      "Markov Centrality" = function(x)markovcent(x),
      "Radiality Centrality" = function(x)radiality(x, weights = NULL),
      "Shortest-Paths Betweenness Centrality" = function(x)betweenness(x),
      "Current-Flow Closeness Centrality" = function(x)closeness.currentflow(x, weights = NULL),
      "Closeness centrality (Latora)" = function(x)closeness.latora(x, weights = NULL),
      "Communicability Betweenness Centrality" = function(x)communibet(x),
      "Community Centrality" = function(x)communitycent(x),
      "Cross-Clique Connectivity" = function(x)crossclique(x),
      "Entropy Centrality" = function(x)entropy(x, weights = NULL),
      "EPC - Edge Percolated Component" = function(x)epc(x),
      "Laplacian Centrality" = function(x)laplacian(x),
      "Leverage Centrality" = function(x)leverage(x),
      "MNC - Maximum Neighborhood Component" = function(x)mnc(x),
      "Hubbell Index" = function(x)hubbell(x, weights = NULL),
      "Semi Local Centrality" = function(x)semilocal(x),
      "Closeness Vitality" = function(x)closeness.vitality(x, weights = NULL),
      "Residual Closeness Centrality" = function(x)closeness.residual(x, weights = NULL),
      "Stress Centrality" = function(x)stresscent(y),
      "Load Centrality" = function(x)loadcent(y),
      "Flow Betweenness Centrality" = function(x)flowbet(y),
      "Information Centrality" = function(x)infocent(y),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = NULL)
    )

    if (!is.null(include)){

      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), include)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
    else{
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), except)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
  }

  if(is_directed(x) && !is_weighted(x) ) {

    centrality_funcs <- list(
      "Alpha Centrality" = function(x) alpha.centrality(x, weights = NULL),
      "Bonacich power centralities of positions" = function(x)bonpow(x),
      "Page Rank" = function(x)page_rank(x)$vector,
      "Average Distance" = function(x)averagedis(x, weights = NULL),
      "Barycenter Centrality" = function(x)barycenter(x, weights = NULL),
      "BottleNeck Centrality" = function(x)bottleneck(x),
      "Centroid value" = function(x)centroid(x, weights = NULL),
      "Closeness Centrality (Freeman)" = function(x)closeness.freeman(x, weights = NULL),
      "ClusterRank" = function(x)clusterrank(x),
      "Decay Centrality" = function(x)decay(x, weights = NULL),
      "Degree Centrality" = function(x)centr_degree(x)$res,
      "Diffusion Degree" = function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component" = function(x)dmnc(x),
      "Eccentricity Centrality" = function(x)eccentricity(x),
      "eigenvector centralities" = function(x)eigen_centrality(x, weights = NULL)$vector,
      "K-core Decomposition" = function(x)coreness(x),
      "Geodesic K-Path Centrality" = function(x)geokpath(x, weights = NULL),
      "Katz Centrality (Katz Status Index)" = function(x)katzcent(x),
      "Kleinberg's authority centrality scores" = function(x)authority_score(x, weights = NULL)$vector,
      "Kleinberg's hub centrality scores" = function(x)hub_score(x, weights = NULL)$vector,
      "clustering coefficient" = function(x)transitivity(x, weights = NULL,type="local"),
      "Lin Centrality" = function(x)lincent(x, weights = NULL),
      "Lobby Index (Centrality)" = function(x)lobby(x),
      "Markov Centrality" = function(x)markovcent(x),
      "Radiality Centrality" = function(x)radiality(x, weights = NULL),
      "Shortest-Paths Betweenness Centrality" = function(x)betweenness(x),
      "Current-Flow Closeness Centrality" = function(x)closeness.currentflow(x, weights = NULL),
      "Closeness centrality (Latora)" = function(x)closeness.latora(x, weights = NULL),
      "Communicability Betweenness Centrality" = function(x)communibet(x),
      "Community Centrality" = function(x)communitycent(x),
      "Cross-Clique Connectivity" = function(x)crossclique(x),
      "Entropy Centrality" = function(x)entropy(x, weights = NULL),
      "EPC - Edge Percolated Component" = function(x)epc(x),
      "Laplacian Centrality" = function(x)laplacian(x),
      "Leverage Centrality" = function(x)leverage(x),
      "MNC - Maximum Neighborhood Component" = function(x)mnc(x),
      "Hubbell Index" = function(x)hubbell(x, weights = NULL),
      "Semi Local Centrality" = function(x)semilocal(x),
      "Closeness Vitality" = function(x)closeness.vitality(x, weights = NULL),
      "Residual Closeness Centrality" = function(x)closeness.residual(x, weights = NULL),
      "Stress Centrality" = function(x)stresscent(y),
      "Load Centrality" = function(x)loadcent(y),
      "Flow Betweenness Centrality" = function(x)flowbet(y),
      "Information Centrality" = function(x)infocent(y),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = NULL)
    )

    if (!is.null(include)){

      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), include)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
    else{
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), except)]

      n <- names(centrality_funcs)

      warningsText <- ""
      result <- lapply(setNames(n, n),
                       function(functionName, x) {
                         f <- centrality_funcs[[functionName]]
                         tryCatch(f(x),
                                  error = function(e) {
                                    warningsText <- paste0(warningsText,
                                                           "\nError in ", functionName, ":\n", e$message)
                                    return(NULL)
                                  })
                       }, x)

      if (nchar(warningsText) > 0)
        warning(warningsText)
      return(result)
    }
  }
}


#' @title Ranking centrality measure based on contributions
#'
#' @description This function demonstrates ranks of centrality measures in order of information levels.
#' @param x a list containg the computed centrality values
#' @param scale.unit	a boolean constant, whether data should be scaled to unit
#' variance(default=TRUE)
#' @param cut.off The intensity that must be exceeded in cumulative percentage of variance of eigen values.(default=80)
#' @param ncp	number of dimensions in final results (default=5)
#' @param graph	a boolean constant, whether the graph shoul be displayed
#' @param axes a length 2 vector describing the number of components to plot(default=c(1,2))
#' @details This function represents centralities in the ranking list based on variable contribution to
#' make principal components.
#' PCA is a method for drawing out important variables from a data set.
#' It helps user to reduced the dimensions in high dimensional data.
#' It is more common to use for more than 3 dimensional datasets.
#' @seealso \code{\link[FactoMineR]{PCA}}
#' @return a plot illustrating significant centralities in the order of contribution
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Husson, F., Lê, S., & Pagès, J. (2010). Exploratory Multivariate Analysis by Example using R. Chapman & Hall/CRC Computer Science & Data Analysis, 40(April), 240.
#'
#' http://www.sthda.com/english/
#'
#' @export
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_contrib
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_grey
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_rect

pca_centralities <- function( x , scale.unit = TRUE, cut.off = 80, ncp = 5,graph = FALSE, axes = c(1,2)){

  x <- x[!sapply(x,is.null)]

  x <- as.data.frame(x)

  x <- na.omit(x)

  res_pca <- PCA(x, scale.unit = scale.unit, ncp = ncp, graph = graph, axes = axes)

  PCs <- table(res_pca$eig[,"cumulative percentage of variance"] > cut.off)["FALSE"] + 1
  if(PCs %in% NA) PCs=1

  fviz_contrib(res_pca, choice = "var", axes = 1:PCs, color = "black", fill = "turquoise")+
    labs(x="\nCentrality measures", y="Contributions\n") +
    theme_grey() +ggtitle("Contribution of variables via PCA")+
    theme(axis.text.x=element_text(angle=45, vjust=0.5),plot.title = element_text(hjust = 0.5)
          ,panel.border = element_rect(colour = "black", fill=NA, size=1.5)
    )
}

#' @title Graph visualization based on a specific centrality measure
#'
#' @description This function demonstrates the input graph in which the size of nodes
#' indicates calculated centrality value.
#' @param x an igraph object
#' @param computed_centrality_value A vector containing the values of calculated centrality measure for each node.
#' @param centrality.type	The type of centrality which should be calculated.
#' @details
#' This function represents the graph in which size of nodes are based on computed centrality value. If the values of wanted centrality
#' measure were computed then by placing them in computed_centrality_value argument to use it for drawing the plot. Otherwise, by only giving the
#' name of favorite centrality measure in centrality.type argument, this function will calculate it and  then demonstrates the corresponding graph.
#' @return a plot illustrating the graph
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @export
#' @importFrom igraph is.directed
#' @importFrom igraph alpha.centrality
#' @importFrom igraph bonpow
#' @importFrom igraph constraint
#' @importFrom igraph centr_degree
#' @importFrom igraph eccentricity
#' @importFrom igraph eigen_centrality
#' @importFrom igraph coreness
#' @importFrom igraph authority_score
#' @importFrom igraph hub_score
#' @importFrom igraph transitivity
#' @importFrom igraph page_rank
#' @importFrom igraph betweenness
#' @importFrom igraph subgraph.centrality
#' @importFrom igraph strength
#' @importFrom sna flowbet
#' @importFrom sna infocent
#' @importFrom sna loadcent
#' @importFrom sna stresscent
#' @importFrom sna graphcent
#' @importFrom centiserve topocoefficient
#' @importFrom centiserve closeness.currentflow
#' @importFrom centiserve closeness.latora
#' @importFrom centiserve communibet
#' @importFrom centiserve communitycent
#' @importFrom centiserve crossclique
#' @importFrom centiserve entropy
#' @importFrom centiserve epc
#' @importFrom centiserve laplacian
#' @importFrom centiserve leverage
#' @importFrom centiserve mnc
#' @importFrom centiserve hubbell
#' @importFrom centiserve semilocal
#' @importFrom centiserve closeness.vitality
#' @importFrom centiserve closeness.residual
#' @importFrom centiserve lobby
#' @importFrom centiserve markovcent
#' @importFrom centiserve radiality
#' @importFrom centiserve lincent
#' @importFrom centiserve geokpath
#' @importFrom centiserve katzcent
#' @importFrom centiserve diffusion.degree
#' @importFrom centiserve dmnc
#' @importFrom centiserve centroid
#' @importFrom centiserve closeness.freeman
#' @importFrom centiserve clusterrank
#' @importFrom centiserve decay
#' @importFrom centiserve barycenter
#' @importFrom centiserve bottleneck
#' @importFrom centiserve averagedis
#' @importFrom igraph layout_in_circle
#' @importFrom graphics plot
#' @importFrom igraph strength

visualize_graph <- function( x , computed_centrality_value=NULL , centrality.type="Degree Centrality"){

  if (is.null(computed_centrality_value)){

    if (!("igraph" %in% class(x) && is_connected(x) )) stop("The input is not an igraph object
                                                          or may not be connected.")

    y <- as_edgelist(x)
    y <- network(y)

    centrality_funcs <- list(
      "subgraph centrality scores"=function(x)subgraph.centrality(x),
      "Topological Coefficient"=function(x)topocoefficient(x),
      "Alpha Centrality"=function(x) alpha.centrality(x),
      "Bonacich power centralities of positions"=function(x)bonpow(x),
      "Burt's Constraint"=function(x)constraint(x),
      "Page Rank"=function(x)page_rank(x)$vector,
      "Average Distance"=function(x)averagedis(x),
      "Barycenter Centrality"=function(x)barycenter(x),
      "BottleNeck Centrality"=function(x)bottleneck(x),
      "Centroid value"=function(x)centroid(x),
      "Closeness Centrality (Freeman)"=function(x)closeness.freeman(x),
      "ClusterRank"=function(x)clusterrank(x),
      "Decay Centrality"=function(x)decay(x),
      "Degree Centrality"=function(x)centr_degree(x)$res,
      "Diffusion Degree"=function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component"=function(x)dmnc(x),
      "Eccentricity Centrality"=function(x)eccentricity(x),
      "eigenvector centralities"=function(x)eigen_centrality(x)$vector,
      "K-core Decomposition"=function(x)coreness(x),
      "Geodesic K-Path Centrality"=function(x)geokpath(x),
      "Katz Centrality (Katz Status Index)"=function(x)katzcent(x),
      "Kleinberg's authority centrality scores"=function(x)authority_score(x)$vector,
      "Kleinberg's hub centrality scores"=function(x)hub_score(x)$vector,
      "clustering coefficient"=function(x)transitivity(x,type="local"),
      "Lin Centrality"=function(x)lincent(x),
      "Lobby Index (Centrality)"=function(x)lobby(x),
      "Markov Centrality"=function(x)markovcent(x),
      "Radiality Centrality"=function(x)radiality(x),
      "Shortest-Paths Betweenness Centrality"=function(x)betweenness(x),
      "Current-Flow Closeness Centrality"=function(x)closeness.currentflow(x),
      "Closeness centrality (Latora)"=function(x)closeness.latora(x),
      "Communicability Betweenness Centrality"=function(x)communibet(x),
      "Community Centrality"=function(x)communitycent(x),
      "Cross-Clique Connectivity"=function(x)crossclique(x),
      "Entropy Centrality"=function(x)entropy(x),
      "EPC - Edge Percolated Component"=function(x)epc(x),
      "Laplacian Centrality"=function(x)laplacian(x),
      "Leverage Centrality"=function(x)leverage(x),
      "MNC - Maximum Neighborhood Component"=function(x)mnc(x),
      "Hubbell Index"=function(x)hubbell(x),
      "Semi Local Centrality"=function(x)semilocal(x),
      "Closeness Vitality"=function(x)closeness.vitality(x),
      "Residual Closeness Centrality"=function(x)closeness.residual(x),
      "Stress Centrality"=function(x)stresscent(y),
      "Load Centrality"=function(x)loadcent(y),
      "Flow Betweenness Centrality"=function(x)flowbet(y),
      "Information Centrality"=function(x)infocent(y),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = NULL),
      "Weighted Vertex Degree" = function(x)strength(x, vids = V(x), mode ="all", weights = NULL)

    )
    centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), centrality.type)]

    n <- names(centrality_funcs)

    warningsText <- ""
    result <- lapply(setNames(n, n),
                     function(functionName, x) {
                       f <- centrality_funcs[[functionName]]
                       tryCatch(f(x),
                                error = function(e) {
                                  warningsText <- paste0(warningsText,
                                                         "\nError in ", functionName, ":\n", e$message)
                                  return(NULL)
                                })
                     }, x)

    if (nchar(warningsText) > 0)
      warning(warningsText)

    result <- result[!sapply(result,is.null)]

    result <- as.data.frame(result)

    result <- scale(result, center = FALSE, scale = TRUE)

    plot(x, vertex.size=result*10, layout= layout_in_circle(x,order(result*10) ) ,
         vertex.color="turquoise",vertex.frame.color="orange2")
  }

  else{

    computed_centrality_value<-scale(computed_centrality_value, center = FALSE, scale = TRUE)

    plot(x, vertex.size=computed_centrality_value*10,
         layout= layout_in_circle(x,order(computed_centrality_value*10) ) ,
         vertex.color="turquoise",vertex.frame.color="orange2")

  }

}

#' @title Pairwise association plot between centrality measures
#'
#' @description This function computes regression between pair of centrality measures
#' to show more details of association among them.
#' @param x a vector containing a centrality measure as independent variable
#' @param y a vector containing a centrality measure as dependent variable
#' @param scale Whether the centrality values should be scaled or not
#' @details
#' This function applies regression analysis on two different centrality values in order to find out the corresponding association
#' between them. Regression analysis is a kind of statitiscal method for approximation the association
#' between variables.It asserts that the value of dependent variable changes when the
#' value of independent variable varies.
#' @return The regression plot, and the values resulted by the regression process.
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' CHAMBERS, & M., J. (1992). Statistical Models in S. Wadsworth. Pacific Grove, California. Retrieved from
#'
#' Wilkinson, G. N., & Rogers, C. E. (1973). Symbolic Description of Factorial Models for Analysis of Variance. Applied Statistics, 22(3), 392.
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 stat_summary
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 mean_cl_normal
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_point

visualize_association <- function( x , y, scale=TRUE){

  xname <- substitute(x)
  yname <- substitute(y)

  df <- as.data.frame(cbind(x,y))

  names(df) <- c(xname, yname)

  if (scale%in%TRUE){

    df <- scale(df, center = TRUE, scale = scale)

    df <- as.data.frame(df)

    visualization <- ggplot(data=df,aes(x,y)) +
      stat_summary(fun.data=mean_cl_normal,geom="point",fill="skyblue3",shape=21,colour="black",  size = 3) +
      geom_smooth(method='lm',colour = 'red') + xlab(xname) +
      ylab(yname)+geom_point(shape=21, fill="blue", color="darkred", size=3)

    linear.regression <- lm(df[,2]~df[,1])
    return (list(linear.regression=linear.regression, visualization=visualization))

  }

  else        {visualization <- ggplot(data=df,aes(x,y)) +
    stat_summary(fun.data=mean_cl_normal,geom="point",fill="skyblue3",shape=21,colour="black",  size = 3) +
    geom_smooth(method='lm',colour = 'red') + xlab(xname) +
    ylab(yname)+geom_point(shape=21, fill="blue", color="darkred", size=3)

  linear.regression <- lm(y~x)
  return (list(linear.regression=linear.regression, visualization=visualization))}

}


#' @title Pairwise correlation plot between two centrality measures
#'
#' @description This function computes and plots correlation between pair of centrality
#'  measures and histogram plots.
#' @param x a vector containing a centrality measure
#' @param y a vector containing another centrality measure
#' @param scale Whether the centrality values should be scaled or not
#' @details
#' This function illustrates the correlation value between two centrality measures and their
#' corresponding scatterplot and histograms.
#' @seealso \code{\link[GGally]{ggpairs}}
#' @return The correlation plot
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Emerson, J. W., Green, W. A., Schloerke, B., Crowley, J., Cook, D., Hofmann, H., & Wickham, H. (2013). The Generalized Pairs Plot. Journal of Computational and Graphical Statistics, 22(1), 79–91.
#' @export
#' @importFrom GGally ggpairs

visualize_pair_correlation <- function( x , y, scale=TRUE){

  xname <- substitute(x)
  yname <- substitute(y)

  df <- as.data.frame(cbind(x,y))

  names(df) <- c(xname, yname)

  if (scale%in%TRUE){

    df <- scale(df, center = TRUE, scale = scale)

    df <- as.data.frame(df)

    ggpairs(df, columns=1:2)

  }

  else  ggpairs(df, columns=1:2)
}

#' @title Heatmap plot between centrality measures
#'
#' @description This function draws heatmap between pair of centrality measures
#' @param x a list indicating calculated centrality measures
#' @param scale Whether the centrality values should be scaled or not
#' @details
#' This function illustrates the heatmap plot of computed centrality measures.
#' @return The correlation plot
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @export
#' @importFrom pheatmap pheatmap

visualize_heatmap <- function( x, scale = TRUE  ){

  if(scale%in%TRUE){

    x <- x[!sapply(x,is.null)]

    d <- as.data.frame(x)

    x <- scale(d, center = TRUE, scale = scale)

    rownames(x)<- rownames(d)

    pheatmap(x, legend = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8)


  }

  else{

    pheatmap(x, legend = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8)

  }

}

#' @title Correlation plot between centrality measures
#'
#' @description This function draw correlation plot between pair of centrality measures
#' @param x a list indicating calculated centrality measures
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @param method a character string describing the type of correlation coefficient (or covariance)
#' to be computed. The proper values are "pearson", "kendall", or "spearman". (default="pearson")
#' @details
#' This function illustrates pairwise correlation plot of computed centrality measures.
#' The names of centralities shown in the result plot is abbreviated and compelete names
#' can be seen in "proper_centralities" function.
#' Colors from red to blue indicate the intensity of correlation value. If two centrality
#' measures have an inverse relationship then their correspnding color in plot have to be red
#' and vice versa.
#'
#' @seealso \code{\link[GGally]{ggpairs}}
#' @return The pairwise correlation plot
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @export
#' @importFrom corrplot corrplot.mixed

visualize_correlations <- function(x, scale=TRUE,method = "pearson"){

  if (scale%in%TRUE){

    x <- x[!sapply(x,is.null)]

    d <- as.data.frame(x)

    x <- scale(d, center = TRUE, scale = scale)

    name <- unlist(lapply(strsplit(colnames(x), '.', fixed = TRUE), '[', 1))

    colnames(x) <- name

    M <- cor(x,method = method)

    corrplot.mixed(M,tl.cex	=0.75)

  }

  else{

    x <- x[!sapply(x,is.null)]

    d <- as.data.frame(x)

    x <- as.matrix(d)

    name <- unlist(lapply(strsplit(colnames(x), '.', fixed = TRUE), '[', 1))

    colnames(x) <- name

    M <- cor(x,method = method)

    corrplot.mixed(M,tl.cex	=0.75)

  }

}

#' @title Dendrogram plot among centrality measures
#'
#' @description This function demonstrates the vertice dendrogram of a graph based
#' on a centrality type.
#' @param x an igraph object
#' @param centrality.type	The type of centrality which should be considered.(default="Degree Centrality")
#' @param computed_centrality_value A vector containing the values of calculated centrality measure for each node.(default=NULL)
#' @param k number of clusters(default=4)
#' @details
#' This function represents node dendrogram of a graph based on a centrality measure.
#' If the favor centrality is not computed yet, by specifying the name of that it will compute it and show the result.
#' @seealso \code{\link[stats]{dendrogram}}
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Galili, T. (2015). dendextend: an R package for visualizing, adjusting and comparing trees of hierarchical clustering. Bioinformatics, 31(22), 3718–3720.
#'
#' @export
#' @importFrom igraph alpha.centrality
#' @importFrom igraph bonpow
#' @importFrom igraph constraint
#' @importFrom igraph centr_degree
#' @importFrom igraph eccentricity
#' @importFrom igraph eigen_centrality
#' @importFrom igraph coreness
#' @importFrom igraph authority_score
#' @importFrom igraph hub_score
#' @importFrom igraph transitivity
#' @importFrom igraph page_rank
#' @importFrom igraph betweenness
#' @importFrom igraph subgraph.centrality
#' @importFrom igraph strength
#' @importFrom sna flowbet
#' @importFrom sna infocent
#' @importFrom sna loadcent
#' @importFrom sna stresscent
#' @importFrom sna graphcent
#' @importFrom centiserve topocoefficient
#' @importFrom centiserve closeness.currentflow
#' @importFrom centiserve closeness.latora
#' @importFrom centiserve communibet
#' @importFrom centiserve communitycent
#' @importFrom centiserve crossclique
#' @importFrom centiserve entropy
#' @importFrom centiserve epc
#' @importFrom centiserve laplacian
#' @importFrom centiserve leverage
#' @importFrom centiserve mnc
#' @importFrom centiserve hubbell
#' @importFrom centiserve semilocal
#' @importFrom centiserve closeness.vitality
#' @importFrom centiserve closeness.residual
#' @importFrom centiserve lobby
#' @importFrom centiserve markovcent
#' @importFrom centiserve radiality
#' @importFrom centiserve lincent
#' @importFrom centiserve geokpath
#' @importFrom centiserve katzcent
#' @importFrom centiserve diffusion.degree
#' @importFrom centiserve dmnc
#' @importFrom centiserve centroid
#' @importFrom centiserve closeness.freeman
#' @importFrom centiserve clusterrank
#' @importFrom centiserve decay
#' @importFrom centiserve barycenter
#' @importFrom centiserve bottleneck
#' @importFrom centiserve averagedis
#' @importFrom dendextend circlize_dendrogram
#' @importFrom dendextend color_branches
#' @importFrom dendextend color_labels
#' @importFrom dendextend highlight_branches_col
#' @importFrom dendextend set
#' @importFrom viridis viridis
#' @importFrom stats as.dendrogram
#' @importFrom stats as.dendrogram
#' @importFrom igraph V
#' @importFrom stats cor
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom stats reorder

visualize_dendrogram <- function( x, centrality.type="Degree Centrality", computed_centrality_value=NULL , k=4){

  if (is.null(computed_centrality_value)){

    if (!("igraph" %in% class(x) && is_connected(x) )) stop("The input is not an igraph object
                                                          or may not be connected.")

    y <- as_edgelist(x)
    y <- network(y)

    centrality_funcs <- list(
      "subgraph centrality scores"=function(x)subgraph.centrality(x),
      "Topological Coefficient"=function(x)topocoefficient(x),
      "Alpha Centrality"=function(x) alpha.centrality(x),
      "Bonacich power centralities of positions"=function(x)bonpow(x),
      "Burt's Constraint"=function(x)constraint(x),
      "Page Rank"=function(x)page_rank(x)$vector,
      "Average Distance"=function(x)averagedis(x),
      "Barycenter Centrality"=function(x)barycenter(x),
      "BottleNeck Centrality"=function(x)bottleneck(x),
      "Centroid value"=function(x)centroid(x),
      "Closeness Centrality (Freeman)"=function(x)closeness.freeman(x),
      "ClusterRank"=function(x)clusterrank(x),
      "Decay Centrality"=function(x)decay(x),
      "Degree Centrality"=function(x)centr_degree(x)$res,
      "Diffusion Degree"=function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component"=function(x)dmnc(x),
      "Eccentricity Centrality"=function(x)eccentricity(x),
      "eigenvector centralities"=function(x)eigen_centrality(x)$vector,
      "K-core Decomposition"=function(x)coreness(x),
      "Geodesic K-Path Centrality"=function(x)geokpath(x),
      "Katz Centrality (Katz Status Index)"=function(x)katzcent(x),
      "Kleinberg's authority centrality scores"=function(x)authority_score(x)$vector,
      "Kleinberg's hub centrality scores"=function(x)hub_score(x)$vector,
      "clustering coefficient"=function(x)transitivity(x,type="local"),
      "Lin Centrality"=function(x)lincent(x),
      "Lobby Index (Centrality)"=function(x)lobby(x),
      "Markov Centrality"=function(x)markovcent(x),
      "Radiality Centrality"=function(x)radiality(x),
      "Shortest-Paths Betweenness Centrality"=function(x)betweenness(x),
      "Current-Flow Closeness Centrality"=function(x)closeness.currentflow(x),
      "Closeness centrality (Latora)"=function(x)closeness.latora(x),
      "Communicability Betweenness Centrality"=function(x)communibet(x),
      "Community Centrality"=function(x)communitycent(x),
      "Cross-Clique Connectivity"=function(x)crossclique(x),
      "Entropy Centrality"=function(x)entropy(x),
      "EPC - Edge Percolated Component"=function(x)epc(x),
      "Laplacian Centrality"=function(x)laplacian(x),
      "Leverage Centrality"=function(x)leverage(x),
      "MNC - Maximum Neighborhood Component"=function(x)mnc(x),
      "Hubbell Index"=function(x)hubbell(x),
      "Semi Local Centrality"=function(x)semilocal(x),
      "Closeness Vitality"=function(x)closeness.vitality(x),
      "Residual Closeness Centrality"=function(x)closeness.residual(x),
      "Stress Centrality"=function(x)stresscent(y),
      "Load Centrality"=function(x)loadcent(y),
      "Flow Betweenness Centrality"=function(x)flowbet(y),
      "Information Centrality"=function(x)infocent(y),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = NULL),
      "Weighted Vertex Degree" = function(x)strength(x, vids = V(x), mode ="all", weights = NULL)

    )
    centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), centrality.type)]

    n <- names(centrality_funcs)

    warningsText <- ""
    result <- lapply(setNames(n, n),
                     function(functionName, x) {
                       f <- centrality_funcs[[functionName]]
                       tryCatch(f(x),
                                error = function(e) {
                                  warningsText <- paste0(warningsText,
                                                         "\nError in ", functionName, ":\n", e$message)
                                  return(NULL)
                                })
                     }, x)

    if (nchar(warningsText) > 0)
      warning(warningsText)

    result <-result[!sapply(result,is.null)]

    result <- as.data.frame(result)

    rownames(result) = V(x)$name

    result <- scale(result, center = TRUE, scale = TRUE)

    dend <- result%>% dist %>% hclust %>% as.dendrogram %>%
      highlight_branches_col(viridis(100))  %>%
      set("branches_k_color", k=3)%>%
      set("labels_colors")%>%set("nodes_pch", 20)

    circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)

  }

  else{

    computed_centrality_value <- scale(computed_centrality_value, center = FALSE, scale = TRUE)

    names(computed_centrality_value) = V(x)$name
    dend <- computed_centrality_value%>% dist %>% hclust %>% as.dendrogram %>%
      highlight_branches_col(viridis(100))  %>%  set("branches_k_color", k=3)%>%
      set("labels_colors")%>%set("nodes_pch", 20)

    circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)


  }

}


#' @title Print pairwise association plot among centrality measures into a file
#'
#' @description This function prints regression plot between pair of centrality measures
#' @param x a vector containing a centrality values as independent variable
#' @param y a vector containing a centrality values as dependent variable
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @param file A character string naming the file to print into. If NULL the result would be printed to the exist directory(default=NULL)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The resulted plot of \code{ \link[CINNA]{visualize_association}}function will be saved in the given directory.
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 stat_summary
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 mean_cl_normal
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggsave

print_visualize_association<-function( x , y, scale = TRUE, file = NULL){

  xname <- substitute(x)
  yname <- substitute(y)
  df <- as.data.frame(cbind(x,y))
  names(df) <- c(xname, yname)

  if (scale%in%TRUE){

    df <- scale(df, center = TRUE, scale = scale)

    df <- as.data.frame(df)

    reg.plot <- ggplot(data=df,aes(x,y))+stat_summary(fun.data=mean_cl_normal,geom="point",
                                                      fill="skyblue3",shape=21,colour="black",  size = 3)+
      geom_smooth(method='lm',colour = 'red')+xlab(xname) +
      ylab(yname)

    if (is.null(file)){

      ggsave(filename=paste(getwd(),"/Linear regression plot between two centrality measure.pdf", sep = ""), plot=reg.plot)
    }

    else ggsave(filename=file, plot=reg.plot)

  }

  else {

    reg.plot <- ggplot(data=df,aes(x,y))+stat_summary(fun.data=mean_cl_normal,
                                                      geom = "point",fill="skyblue3",shape=21,colour="black",  size = 3)+
      geom_smooth(method='lm',colour = 'red')+xlab(xname) + ylab(yname)

    if (is.null(file)){

      ggsave(filename=paste(getwd(),"/Linear regression plot between two centrality measure.pdf", sep = ""), plot=reg.plot)

    }

    else ggsave(filename=file, plot=reg.plot)

  }
}


#' @title Print pairwise correlation and histogram plots between two centrality measures
#'
#' @description This function prints pairwise correlation of centrality measures and histogram plot.
#' @param x a vector containing a centrality measure
#' @param y a vector containing another centrality measure
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @param file A character string naming the .pdf file to print into. If NULL the result would be printed to the exist directory.
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return The resulted plot of \code{ \link[CINNA]{visualize_pair_correlation}}function will be saved in the given directory.
#' @importFrom GGally ggpairs
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf

print_visualize_pair_correlation <- function( x , y, scale=TRUE, file=NULL){

  xname <- substitute(x)
  yname <- substitute(y)

  df <- as.data.frame(cbind(x,y))

  names(df) <- c(xname, yname)

  if (scale%in%TRUE){

    df<-scale(df, center = TRUE, scale = scale)

    df<-as.data.frame(df)

    if (is.null(file)){

      variables <- paste( as.character(xname), "and", as.character(yname), sep = " ")
      pdf(file=paste(getwd(),"/Pairwise correlation plot between ", variables,".pdf", sep = ""))
      g <-ggpairs(df, columns=1:2)
      print(g)
      dev.off()
    }
    else
      pdf(file=file)
    g <- ggpairs(df, columns=1:2)
    print(g)
    dev.off()
  }

  else{

    if (is.null(file)){
      variables<-paste( as.character(xname), "and", as.character(yname), sep = " ")
      pdf(file=paste(getwd(),"/Pairwise correlation plot between ", variables,".pdf", sep = ""))
      g <- ggpairs(df, columns=1:2)
      print(g)
      dev.off()
    }
    else
      pdf(file=file)
    g <- ggpairs(df, columns=1:2)
    print(g)
    dev.off()
  }
}


#' @title Print the heatmap plot of centrality measures
#'
#' @description This function prints the heatmap plot
#' @param x a list indicating calculated centrality measures
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @param file A character string naming the .pdf file to print into. If NULL
#'  the result would be printed to the exist directory.(default=NULL)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The resulted plot of \code{ \link[CINNA]{visualize_heatmap}}function will be saved in the given directory.
#' @importFrom pheatmap pheatmap

print_visualize_heatmap <- function( x, scale = TRUE, file=NULL ){

  if (scale%in%TRUE){

    x <- x[!sapply(x,is.null)]

    d <- as.data.frame(x)

    x <- scale(d, center = TRUE, scale = scale)

    rownames(x) <- rownames(d)

    if(is.null(file)){

      pdf(file=paste(getwd(),"/Heatmap plot among centrality values.pdf", sep = ""), onefile=FALSE)
      pheatmap(x, legend = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8)
      dev.off()
    }
    else
     pdf(file=file, onefile=FALSE)
     pheatmap(x, legend = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8)
     dev.off()

  }

  else{

    if(is.null(file)){

      pdf(file=paste(getwd(),"/Heatmap plot among centrality values.pdf", sep = ""), onefile=FALSE)
      pheatmap(x, legend = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8)
      dev.off()
    }
    else
    pdf(file=file, onefile=FALSE)
    pheatmap(x, legend = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8)
    dev.off()
  }

}


#' @title Print dendrogram plot of a graph
#'
#' @description This function prints dendrogram plot of a graph based on predefined
#' centrality measures.
#' @param x an igraph object
#' @param centrality.type	The type of centrality which should be calculated(default="Degree Centrality")
#' @param computed_centrality_value A vector containing the values of calculated
#' centrality measure for each node(default=NULL)
#' @param k number of clusters
#' @param file A character string naming the .pdf file to print into.
#' If NULL the result would be printed to the exist directory.(default=NULL)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The resulted plot of \code{ \link[CINNA]{visualize_dendrogram}}function will be saved in the given directory.#' @importFrom igraph alpha.centrality
#' @importFrom igraph bonpow
#' @importFrom igraph constraint
#' @importFrom igraph centr_degree
#' @importFrom igraph eccentricity
#' @importFrom igraph eigen_centrality
#' @importFrom igraph coreness
#' @importFrom igraph authority_score
#' @importFrom igraph hub_score
#' @importFrom igraph transitivity
#' @importFrom igraph page_rank
#' @importFrom igraph betweenness
#' @importFrom igraph subgraph.centrality
#' @importFrom igraph strength
#' @importFrom sna flowbet
#' @importFrom sna infocent
#' @importFrom sna loadcent
#' @importFrom sna stresscent
#' @importFrom sna graphcent
#' @importFrom centiserve topocoefficient
#' @importFrom centiserve closeness.currentflow
#' @importFrom centiserve closeness.latora
#' @importFrom centiserve communibet
#' @importFrom centiserve communitycent
#' @importFrom centiserve crossclique
#' @importFrom centiserve entropy
#' @importFrom centiserve epc
#' @importFrom centiserve laplacian
#' @importFrom centiserve leverage
#' @importFrom centiserve mnc
#' @importFrom centiserve hubbell
#' @importFrom centiserve semilocal
#' @importFrom centiserve closeness.vitality
#' @importFrom centiserve closeness.residual
#' @importFrom centiserve lobby
#' @importFrom centiserve markovcent
#' @importFrom centiserve radiality
#' @importFrom centiserve lincent
#' @importFrom centiserve geokpath
#' @importFrom centiserve katzcent
#' @importFrom centiserve diffusion.degree
#' @importFrom centiserve dmnc
#' @importFrom centiserve centroid
#' @importFrom centiserve closeness.freeman
#' @importFrom centiserve clusterrank
#' @importFrom centiserve decay
#' @importFrom centiserve barycenter
#' @importFrom centiserve bottleneck
#' @importFrom centiserve averagedis
#' @importFrom dendextend circlize_dendrogram
#' @importFrom dendextend color_branches
#' @importFrom dendextend color_labels
#' @importFrom dendextend highlight_branches_col
#' @importFrom viridis viridis

print_visualize_dendrogram <- function( x, centrality.type="Degree Centrality",
                                     computed_centrality_value=NULL , k=4, file=NULL){

  if (is.null(computed_centrality_value)){

    if (!("igraph" %in% class(x) && is_connected(x) )) stop("The input is not an igraph object
                                                       or may not be connected.")
    y <- as_edgelist(x)
    y <- network(y)

    centrality_funcs <- list(
      "subgraph centrality scores"=function(x)subgraph.centrality(x),
      "Topological Coefficient"=function(x)topocoefficient(x),
      "Alpha Centrality"=function(x) alpha.centrality(x),
      "Bonacich power centralities of positions"=function(x)bonpow(x),
      "Burt's Constraint"=function(x)constraint(x),
      "Page Rank"=function(x)page_rank(x)$vector,
      "Average Distance"=function(x)averagedis(x),
      "Barycenter Centrality"=function(x)barycenter(x),
      "BottleNeck Centrality"=function(x)bottleneck(x),
      "Centroid value"=function(x)centroid(x),
      "Closeness Centrality (Freeman)"=function(x)closeness.freeman(x),
      "ClusterRank"=function(x)clusterrank(x),
      "Decay Centrality"=function(x)decay(x),
      "Degree Centrality"=function(x)centr_degree(x)$res,
      "Diffusion Degree"=function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component"=function(x)dmnc(x),
      "Eccentricity Centrality"=function(x)eccentricity(x),
      "eigenvector centralities"=function(x)eigen_centrality(x)$vector,
      "K-core Decomposition"=function(x)coreness(x),
      "Geodesic K-Path Centrality"=function(x)geokpath(x),
      "Katz Centrality (Katz Status Index)"=function(x)katzcent(x),
      "Kleinberg's authority centrality scores"=function(x)authority_score(x)$vector,
      "Kleinberg's hub centrality scores"=function(x)hub_score(x)$vector,
      "clustering coefficient"=function(x)transitivity(x,type="local"),
      "Lin Centrality"=function(x)lincent(x),
      "Lobby Index (Centrality)"=function(x)lobby(x),
      "Markov Centrality"=function(x)markovcent(x),
      "Radiality Centrality"=function(x)radiality(x),
      "Shortest-Paths Betweenness Centrality"=function(x)betweenness(x),
      "Current-Flow Closeness Centrality"=function(x)closeness.currentflow(x),
      "Closeness centrality (Latora)"=function(x)closeness.latora(x),
      "Communicability Betweenness Centrality"=function(x)communibet(x),
      "Community Centrality"=function(x)communitycent(x),
      "Cross-Clique Connectivity"=function(x)crossclique(x),
      "Entropy Centrality"=function(x)entropy(x),
      "EPC - Edge Percolated Component"=function(x)epc(x),
      "Laplacian Centrality"=function(x)laplacian(x),
      "Leverage Centrality"=function(x)leverage(x),
      "MNC - Maximum Neighborhood Component"=function(x)mnc(x),
      "Hubbell Index"=function(x)hubbell(x),
      "Semi Local Centrality"=function(x)semilocal(x),
      "Closeness Vitality"=function(x)closeness.vitality(x),
      "Residual Closeness Centrality"=function(x)closeness.residual(x),
      "Stress Centrality"=function(x)stresscent(y),
      "Load Centrality"=function(x)loadcent(y),
      "Flow Betweenness Centrality"=function(x)flowbet(y),
      "Information Centrality"=function(x)infocent(y),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = NULL),
      "Weighted Vertex Degree" = function(x)strength(x, vids = V(x), mode ="all", weights = NULL)


    )
    centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), centrality.type)]

    n <- names(centrality_funcs)

    warningsText <- ""
    result <- lapply(setNames(n, n),
                     function(functionName, x) {
                       f <- centrality_funcs[[functionName]]
                       tryCatch(f(x),
                                error = function(e) {
                                warningsText <- paste0(warningsText,
                                                         "\nError in ", functionName, ":\n", e$message)
                                  return(NULL)
                                })
                                }, x)

    if (nchar(warningsText) > 0)
      warning(warningsText)

    result <- result[!sapply(result,is.null)]

    result <- as.data.frame(result)

    rownames(result) = V(x)$name

    result <- scale(result, center = TRUE, scale = TRUE)

    dend <- result%>% dist %>% hclust %>% as.dendrogram %>%
      highlight_branches_col(viridis(100))  %>%  set("branches_k_color", k=3)%>%
      set("labels_colors")%>%set("nodes_pch", 20)

    if(is.null(file)){
      pdf(file=paste(getwd(),"/Centrality dendrogram.pdf", sep = ""), onefile=FALSE)
      circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
      dev.off()
    }
    else
      pdf(file=file, onefile=FALSE)
    circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
    dev.off()

  }

  else{

    computed_centrality_value <- scale(computed_centrality_value, center = FALSE, scale = TRUE)

    names(computed_centrality_value) = V(x)$name
    dend <- computed_centrality_value%>% dist %>% hclust %>% as.dendrogram %>%
      highlight_branches_col(viridis(100))  %>%  set("branches_k_color", k=3)%>%
      set("labels_colors")%>%set("nodes_pch", 20)

    if (is.null(file)){
      pdf(file=paste(getwd(),"/Centrality dendrogram.pdf", sep = ""), onefile=FALSE)
      circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
      dev.off()
    }
    else
      pdf(file=file, onefile=FALSE)
    circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
    dev.off()

  }
}

#' @title Print centrality correlation plot
#'
#' @description This function prints a plot including all pairwise correlation between
#' centrality measures
#' @param x a list indicating calculated centrality measures which is the output of
#' "calculate_centralities" function
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @param method   character string describing the type of correlation coefficient (or covariance)
#' to be computed. The proper values are "pearson", "kendall", or "spearman". (default="pearson")
#' @param file A character string naming the .pdf file to print into. If NULL the result would
#' be printed to the exist directory.(default=NULL)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The resulted plot of \code{ \link[CINNA]{visualize_correlations}}function will be saved in the given directory.#' @importFrom igraph alpha.centrality
#' @importFrom corrplot corrplot.mixed

print_visualize_correlations <- function(x, scale=TRUE,method = c("pearson", "kendall", "spearman"),file=NULL){

  if (scale%in%TRUE){

    x <- x[!sapply(x,is.null)]

    d <- as.data.frame(x)

    x <- scale(d, center = TRUE, scale = TRUE)

    name <- unlist(lapply(strsplit(colnames(x), '.', fixed = TRUE), '[', 1))

    colnames(x) <- name

    M <- cor(x,method = "pearson")

    if (is.null(file)){
      pdf(file=paste(getwd(),"/Centrality correlation plot.pdf", sep = ""), onefile=FALSE)
      corrplot.mixed(M,tl.cex	=0.75)
      dev.off()
    }
    else
      pdf(file=file, onefile=FALSE)
    corrplot.mixed(M,tl.cex	=0.75)
    dev.off()
  }
  else{

    x <- x[!sapply(x,is.null)]

    d <- as.data.frame(x)

    x <- as.matrix(d)

    name <- unlist(lapply(strsplit(colnames(x), '.', fixed = TRUE), '[', 1))

    colnames(x) <- name

    M <- cor(x,method = "pearson")

    if (is.null(file)){
      pdf(file=paste(getwd(),"/Centrality correlation plot.pdf", sep = ""), onefile=FALSE)
      corrplot.mixed(M,tl.cex	=0.75)
      dev.off()
    }
    else
      pdf(file=file, onefile=FALSE)
    corrplot.mixed(M,tl.cex	=0.75)
    dev.off()

  }
}

#' @title Print visualized based on centrality values graph
#'
#' @description This function prints visualized based on centrality values graph.
#' @param x an igraph object
#' @param computed_centrality_value A vector containing the values of calculated centrality measure for each node(default=NULL)
#' @param centrality.type	The type of centrality which should be calculated(default="Degree Centrality")
#' @param file A character string naming the .pdf file to print into. If NULL the result would be printed to the exist directory.(default=NULL)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The resulted plot of \code{ \link[CINNA]{visualize_graph}}function will be saved in the given directory.#' @importFrom igraph alpha.centrality
#' @importFrom igraph is.directed
#' @importFrom igraph alpha.centrality
#' @importFrom igraph bonpow
#' @importFrom igraph constraint
#' @importFrom igraph centr_degree
#' @importFrom igraph eccentricity
#' @importFrom igraph eigen_centrality
#' @importFrom igraph coreness
#' @importFrom igraph authority_score
#' @importFrom igraph hub_score
#' @importFrom igraph transitivity
#' @importFrom igraph page_rank
#' @importFrom igraph betweenness
#' @importFrom igraph subgraph.centrality
#' @importFrom igraph strength
#' @importFrom sna flowbet
#' @importFrom sna infocent
#' @importFrom sna loadcent
#' @importFrom sna stresscent
#' @importFrom sna graphcent
#' @importFrom centiserve topocoefficient
#' @importFrom centiserve closeness.currentflow
#' @importFrom centiserve closeness.latora
#' @importFrom centiserve communibet
#' @importFrom centiserve communitycent
#' @importFrom centiserve crossclique
#' @importFrom centiserve entropy
#' @importFrom centiserve epc
#' @importFrom centiserve laplacian
#' @importFrom centiserve leverage
#' @importFrom centiserve mnc
#' @importFrom centiserve hubbell
#' @importFrom centiserve semilocal
#' @importFrom centiserve closeness.vitality
#' @importFrom centiserve closeness.residual
#' @importFrom centiserve lobby
#' @importFrom centiserve markovcent
#' @importFrom centiserve radiality
#' @importFrom centiserve lincent
#' @importFrom centiserve geokpath
#' @importFrom centiserve katzcent
#' @importFrom centiserve diffusion.degree
#' @importFrom centiserve dmnc
#' @importFrom centiserve centroid
#' @importFrom centiserve closeness.freeman
#' @importFrom centiserve clusterrank
#' @importFrom centiserve decay
#' @importFrom centiserve barycenter
#' @importFrom centiserve bottleneck
#' @importFrom centiserve averagedis


print_visualize_graph <- function( x , computed_centrality_value=NULL , centrality.type="Degree Centrality", file=NULL){

  if (is.null(computed_centrality_value)){

    if (!("igraph" %in% class(x) && is_connected(x) )) stop("The input is not an igraph object
                                                       or may not be connected.")

    y <- as_edgelist(x)
    y <- network(y)

    centrality_funcs <- list(
      "subgraph centrality scores"=function(x)subgraph.centrality(x),
      "Topological Coefficient"=function(x)topocoefficient(x),
      "Alpha Centrality"=function(x) alpha.centrality(x),
      "Bonacich power centralities of positions"=function(x)bonpow(x),
      "Burt's Constraint"=function(x)constraint(x),
      "Page Rank"=function(x)page_rank(x)$vector,
      "Average Distance"=function(x)averagedis(x),
      "Barycenter Centrality"=function(x)barycenter(x),
      "BottleNeck Centrality"=function(x)bottleneck(x),
      "Centroid value"=function(x)centroid(x),
      "Closeness Centrality (Freeman)"=function(x)closeness.freeman(x),
      "ClusterRank"=function(x)clusterrank(x),
      "Decay Centrality"=function(x)decay(x),
      "Degree Centrality"=function(x)centr_degree(x)$res,
      "Diffusion Degree"=function(x)diffusion.degree(x),
      "DMNC - Density of Maximum Neighborhood Component"=function(x)dmnc(x),
      "Eccentricity Centrality"=function(x)eccentricity(x),
      "eigenvector centralities"=function(x)eigen_centrality(x)$vector,
      "K-core Decomposition"=function(x)coreness(x),
      "Geodesic K-Path Centrality"=function(x)geokpath(x),
      "Katz Centrality (Katz Status Index)"=function(x)katzcent(x),
      "Kleinberg's authority centrality scores"=function(x)authority_score(x)$vector,
      "Kleinberg's hub centrality scores"=function(x)hub_score(x)$vector,
      "clustering coefficient"=function(x)transitivity(x,type="local"),
      "Lin Centrality"=function(x)lincent(x),
      "Lobby Index (Centrality)"=function(x)lobby(x),
      "Markov Centrality"=function(x)markovcent(x),
      "Radiality Centrality"=function(x)radiality(x),
      "Shortest-Paths Betweenness Centrality"=function(x)betweenness(x),
      "Current-Flow Closeness Centrality"=function(x)closeness.currentflow(x),
      "Closeness centrality (Latora)"=function(x)closeness.latora(x),
      "Communicability Betweenness Centrality"=function(x)communibet(x),
      "Community Centrality"=function(x)communitycent(x),
      "Cross-Clique Connectivity"=function(x)crossclique(x),
      "Entropy Centrality"=function(x)entropy(x),
      "EPC - Edge Percolated Component"=function(x)epc(x),
      "Laplacian Centrality"=function(x)laplacian(x),
      "Leverage Centrality"=function(x)leverage(x),
      "MNC - Maximum Neighborhood Component"=function(x)mnc(x),
      "Hubbell Index"=function(x)hubbell(x),
      "Semi Local Centrality"=function(x)semilocal(x),
      "Closeness Vitality"=function(x)closeness.vitality(x),
      "Residual Closeness Centrality"=function(x)closeness.residual(x),
      "Stress Centrality"=function(x)stresscent(y),
      "Load Centrality"=function(x)loadcent(y),
      "Flow Betweenness Centrality"=function(x)flowbet(y),
      "Information Centrality"=function(x)infocent(y),
      "Harary Centrality" = function(x)graphcent(y, gmode="graph", diag=T, cmode="directed"),
      "Dangalchev Closeness Centrality"= function(x)dangalchev_closeness_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Group Centrality"= function(x)group_centrality(x, vids = V(x)),
      "Harmonic Centrality"= function(x)harmonic_centrality(x, vids = V(x), mode = "all", weights = NULL),
      "Local Bridging Centrality"= function(x)local_bridging_centrality(x, vids = V(x)),
      "Wiener Index Centrality"= function(x)wiener_index_centrality(x, vids = V(x), mode ="all", weights = NULL),
      "Weighted Vertex Degree" = function(x)strength(x, vids = V(x), mode ="all", weights = NULL)


    )
    centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), centrality.type)]

    n <- names(centrality_funcs)

    warningsText <- ""
    result <- lapply(setNames(n, n),
                     function(functionName, x) {
                       f <- centrality_funcs[[functionName]]
                       tryCatch(f(x),
                                error = function(e) {
                                  warningsText <- paste0(warningsText,
                                                         "\nError in ", functionName, ":\n", e$message)
                                  return(NULL)
                                })
                                }, x)

    if (nchar(warningsText) > 0)
      warning(warningsText)

    result<-result[!sapply(result,is.null)]

    result<-as.data.frame(result)


    result<-scale(result, center = FALSE, scale = TRUE)

    pdf(file=paste(getwd(),"/Graph centrality visualization.pdf", sep = ""), onefile=FALSE)
    plot(x, vertex.size=result*10, layout= layout_in_circle(x,order(result*10) ) ,
         vertex.color="turquoise",vertex.frame.color="orange2")
         dev.off()

  }

  else{

    computed_centrality_value<-scale(computed_centrality_value, center = FALSE, scale = TRUE)

    pdf(file=paste(getwd(),"/Graph centrality visualization.pdf", sep = ""), onefile=FALSE)
    plot(x, vertex.size=computed_centrality_value*10, layout= layout_in_circle(x,order(computed_centrality_value*10) ) ,
         vertex.color="turquoise",vertex.frame.color="orange2")
         dev.off()
  }

}


#' @title Summarize PCA result related to centrality measures
#'
#' @description This function summarizes PCA result related to centrality measures.
#' @param x a list containg the computed centrality values
#' @param scale.unit	a boolean constant, whether data should be scaled to unit
#' variance(default=TRUE)
#' @param ncp	number of dimensions in final results (default=5)
#' @return  The result values of \code{ \link[CINNA]{pca_centralities}}function will be saved in the given directory.#' @importFrom igraph alpha.centrality
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @importFrom FactoMineR PCA

summary_pca_centralities <- function( x , scale.unit = TRUE,ncp = 5){

  x <- x[!sapply(x,is.null)]

  x <- as.data.frame(x)

  x <- na.omit(x)

  res_pca <- PCA(x, scale.unit = scale.unit, ncp = ncp, graph = FALSE)

  l <- list(eigen= res_pca$eig, contribution= res_pca$var$contrib)
  print(l)

}

#' @title Summarize component extraction of a graph
#'
#' @description This function summarizes all components of the input which can be an "igraph" object or a "network" object
#' @param x An igraph or a network object
#' @param directed a boolean constant, Whether to create a directed graph(default=TRUE)
#' @param bipartite_proj Whether the bipartite network must be projected or not(default=FALSE)
#' @param num_proj A number which shows the number of projects especifically for bipartite graphs.(default=1)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The result values of \code{ \link[CINNA]{graph_extract_components}}function will be saved in the given directory.#' @importFrom igraph alpha.centrality
#' @importFrom igraph is_igraph
#' @importFrom igraph is_bipartite
#' @importFrom igraph bipartite.projection
#' @importFrom igraph is_simple
#' @importFrom igraph simplify
#' @importFrom igraph clusters
#' @importFrom igraph induced.subgraph
#' @importFrom igraph graph_from_edgelist
#' @importFrom network is.network
#' @importFrom network as.edgelist
#' @importFrom network is.network
#' @importFrom network network

summary_graph_extract_components <- function( x, directed = TRUE, bipartite_proj = FALSE , num_proj = 1){

  if (!("igraph" %in% class(x) || ("network" %in% class(x) ))) stop("The input is not an igraph or a network object")

  if (is_igraph(x)) {

    if (bipartite_proj){

      if (is_bipartite(x)){

        x <- bipartite.projection(x)[[num_proj]]

        if (!is_simple(x))   x<-simplify(x)

        cl <- clusters(x)

        graph_splitting <- function(k, x, cl){
          induced.subgraph(x, cl$membership == k)
        }

        components <- sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

      }

    }

    else{

      if (!is_simple(x))   x<-simplify(x)

      cl <- clusters(x)

      graph_splitting <- function(k, x, cl){
        induced.subgraph(x, cl$membership == k)
      }

      components <- sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

    }

  }

  if( is.network(x)){

    edgelist <- as.edgelist(x)

    x <- graph_from_edgelist(edgelist, directed = directed)

    if (!is_simple(x))  gr <- simplify(x)

    cl <- clusters(x)

    graph_splitting <- function(k, x, cl){
      induced.subgraph(x, cl$membership == k)
    }

    components <- sapply(1:max(cl$membership), graph_splitting, x = x, cl = cl, simplify = FALSE)

  }

  summary.list <- lapply(components, function(x) summary(x))

}

#' @title Summarize centrality measure calculation results
#' @description This function computes  minimum, first quater, median,
#' mean, third qarter and maximum values of computed centrality measures.
#' @param x centrality measure calculation results
#' @return a list concluding summary results for each centrality measure value
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  The result values of \code{ \link[CINNA]{calculate_centralities}}function will be saved in the given directory.#' @importFrom igraph alpha.centrality

summary_calculate_centralities <- function(x){

  x <- x[!sapply(x, is.null)]

  x <- na.omit(x)

  summary_list <- lapply(x, function(x) summary(x))
  return(summary_list)
}

#' @title Summarize t-Distributed Stochastic Neighbor Embedding (t-SNE) on centrality measures
#'
#' @description This function summarize tsne analysis results on centrality measures
#' @param x a list containg the computed cetrality values
#' @param dims	integer; number of the outpu dimensions(default=2)
#' @param perplexity numeric; A flexible measure of the efficient number of neighbors.
#' The performance of SNE is fairly robust to changes in the perplexity, and typical
#' values are between 5 and 50.(default=5)
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' @return It resturns a list containing below values:
#'
#' Y	Matrix containing the new representations for the objects
#'
#' costs	The cost for every object after the final iteration
#'
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @importFrom Rtsne Rtsne


summary_tsne_centralities<-function( x , dims = 2, perplexity = 5, scale = TRUE){

  x <- x[!sapply(x, is.null)]

  x <- as.data.frame(x)

  x <- na.omit(x)

  if (scale%in%TRUE){

    x <- scale(x, center = TRUE, scale = scale)

    x <- x[!duplicated(x), ]

    tsne.Y <- Rtsne(t(x), dims = dims, perplexity = perplexity, check_duplicates = FALSE)$Y

    rownames(tsne.Y)<-colnames(x)

    cost <- Rtsne(t(x), dims = dims, perplexity = perplexity ,check_duplicates = FALSE)$cost

    names(cost) <- colnames(x)

    tsne_cost <- sort(cost)

    res_tsne <- list(tsne.Y = tsne.Y, tsne_cost = tsne_cost)
    return(res_tsne)

  }

  else{

    x <- x[!duplicated(x), ]

    tsne.Y <- Rtsne(t(x), dims = dims, perplexity = perplexity,check_duplicates = FALSE)$Y

    rownames(tsne.Y) <- colnames(x)

    cost <- Rtsne(t(x), dims = dims, perplexity = perplexity,check_duplicates = FALSE)$cost

    names(cost) <- colnames(x)

    tsne_cost <- sort(cost)

    res_tsne <- list(tsne.Y = tsne.Y,tsne_cost = tsne_cost)

    return(res_tsne)

  }
}


#' @title t-Distributed Stochastic Neighbor Embedding (t-SNE) on centrality measures
#'
#' @description This function applies t-SNE, dimensionality reduction algorithm, on
#' centrality measures.
#' @param x a list containg the computed cetrality values
#' @param dims	integer; number of the outpu dimensions(default=2)
#' @param perplexity numeric; A flexible measure of the efficient number of neighbors.
#' The performance of SNE is fairly robust to changes in the perplexity, and typical
#' values are between 5 and 50.(default=5)
#' @param scale Whether the centrality values should be scaled or not(default=TRUE)
#' @details t-SNE is a non-linear dimensionality reduction algorithm used for exploring high-dimensional data. Here, It maps multi-dimensional centrality measure data to less dimensions suitable to work with it.
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' @return It resturns cost plot of tsne results which displays centralities in order of their corresponding costs.
#'
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' van der Maaten, L. (2014). Accelerating t SNE using Tree Based Algorithms. Journal of Machine Learning Research, 15, 3221–3245.
#' Van Der Maaten, L. J. P., & Hinton, G. E. (2008). Visualizing high dimensional data using t sne. Journal of Machine Learning Research, 9, 2579–2605.
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme

tsne_centralities <- function( x , dims = 2, perplexity = 5, scale = TRUE){

  x <- x[!sapply(x,is.null)]

  x <- as.data.frame(x)

  x <- na.omit(x)

  if (scale%in%TRUE){

    x <- scale(x, center = TRUE, scale = scale)

    x <- x[!duplicated(x), ]

    cost <- Rtsne(t(x), dims = dims, perplexity = perplexity,check_duplicates = FALSE)$cost

    names(cost) <- colnames(x)

    df <- data.frame(sort(cost))

    ggplot(data=df, aes(x=reorder(rownames(df), -cost), y=cost, fill=rownames(df))) +
      geom_bar(colour="black", fill="turquoise", width=.8, stat="identity") +
      guides(fill=FALSE) +
      xlab("Centrality measures") + ylab("tsne costs") +
      ggtitle("tsne cost results")+  theme(axis.text.x=element_text(angle=45, vjust=0.5),plot.title = element_text(hjust = 0.5)
                                            ,panel.border = element_rect(colour = "black", fill=NA, size=1.5))

  }

  else{

    x <- x[!duplicated(x), ]

    cost <- Rtsne(t(x), dims = dims, perplexity = perplexity,check_duplicates = FALSE)$cost

    names(cost) <- colnames(x)

    df <- data.frame(sort(cost))

    ggplot(data=df, aes(x=reorder(rownames(df), -cost), y=cost, fill=rownames(df))) +
      geom_bar(colour="black", fill="turquoise", width=.8, stat="identity") +
      guides(fill=FALSE) +
      xlab("Centrality measures") + ylab("tsne costs") +
      ggtitle("tsne cost results")+  theme(axis.text.x=element_text(angle=45, vjust=0.5),
                                            plot.title = element_text(hjust = 0.5)
                                            ,panel.border = element_rect(colour = "black", fill=NA, size=1.5))

  }
}

#' @title Print computed centrality measures results into a file
#'
#' @description This function prints all centrality measure results into a file
#' @param x a list containing the centrality measure values
#' @param file A character string naming the .pdf file to print into.
#' If NULL the result would be printed to the exist directory.(default=NULL)
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @return  Print out  \code{ \link[CINNA]{calculate_centralities}}function will be saved in the given directory.
#' @importFrom utils write.csv

print_calculate_centralities<- function(x , file = NULL){

  if (is.null(file)){

    c = do.call("cbind",x)
    write.csv(c,"Computed.centrality.values.txt")
  }
  else{
    c = do.call("cbind",x)
    write.csv(c,"Computed.centrality.values.txt", file = file)
  }
}


#' @title Dangalchev Closeness Centrality
#'
#' @description This function computes Dangalchev Closeness Centrality.
#' This can be access by computing a network resistance. More specifically, it measures the closeness by removing nodes and edges.
#' The evalauation of this measure of closeness will be easier
#' and this can be useful for unconnected graphs too.
#'
#' @param x An igraph or a network object
#' @param vids Nodes to be considered in the calculation
#' @param mode A Character value, indicating whether the shortest paths "in" or "out" of the nodes in the directed graphs should be considered. For undirected graphs we use "all".
#' @param weights Numeric vector indicating weights of the edges
#'
#' @seealso \code{\link[centiserve]{closeness.residual}}
#' @return
#' a vector including centrality values for each node
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' DANGALCHEV, C. 2006. Residual closeness in networks. Physica A: Statistical Mechanics and its Applications, 365, 556-564. DOI: 10.1016/j.physa.2005.12.020
#'
#' @examples
#'
#' data(zachary)
#'
#' dangalchev_closeness_centrality(zachary)
#'
#' @export
#' @importFrom igraph distances
#' @importFrom igraph V
#' @importFrom igraph is_igraph
#' @importFrom igraph is_named
#' @importFrom igraph getIgraphOpt
#' @importFrom intergraph asIgraph

dangalchev_closeness_centrality<-function (x, vids = V(x), mode = c("all", "out", "in"), weights = NULL){

  if (!("igraph" %in% class(x) || "network" %in% class(x))) stop("The input is not an igraph or a
                                                            network object")

  if (is_igraph(x)){

    distMat<-1/2^distances(x)

    res <- rowSums(distMat) - diag(distMat)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(res) <- V(x)$name[vids]

    }
  }

  else{

    x<-asIgraph(x)

    distMat<-1/2^distances(x)

    res <- rowSums(distMat) - diag(distMat)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(res) <- V(x)$name[vids]

    }
  }
  return(res)
}


#' @title Local Bridging Centrality
#'
#' @description This function computes Local Bridging Centrality. This classifies nodes regarding their structural links among the dense components.
#'
#' @param x An igraph or a network object
#' @param vids Nodes to be considered in the calculation
#'
#' @seealso \code{\link[igraph]{betweenness}}
#' @return
#' a vector including centrality values for each node
#'
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Macker, J.P., 2016, November. An improved local bridging centrality model for distributed network analytics. In Military Communications Conference, MILCOM 2016-2016 IEEE (pp. 600-605). IEEE. DOI: 10.1109/MILCOM.2016.7795393
#'
#' @examples
#'
#' data(zachary)
#'
#' local_bridging_centrality(zachary)
#'
#' @export
#' @importFrom igraph degree
#' @importFrom igraph neighbors
#' @importFrom igraph V
#' @importFrom igraph is_igraph
#' @importFrom igraph is_named
#' @importFrom igraph getIgraphOpt
#' @importFrom intergraph asIgraph

local_bridging_centrality<-function (x, vids = V(x)){

  if (!("igraph" %in% class(x)|| "network" %in% class(x))) stop("The input is not an igraph or a
                                                            network object")

  f <- function(v){
    (1/degree(x,v))/sum(1/degree(x,neighbors(x,v,'all')))}

  if (is_igraph(x)){

    results <- sapply(V(x),f)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(results) <- V(x)$name[vids]
    }
  }

  else{

    x<-asIgraph(x)

    results <- sapply(V(x),f)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(results) <- V(x)$name[vids]
    }
  }
  return(results)
}

#' @title Wiener Index Centrality
#'
#' @description This function computes Wiener Index Centrality.
#' The Wiener index computes the sum of the all shortest paths
#' between a node v and all other related nodes in the graph.
#' Fundementally, it's like to the closeness but here since the reciprocal is not computed, the value has the opposite meaning.
#'
#' @param x An igraph or a network object
#' @param vids Nodes to be considered in the calculation
#' @param mode A Character value, indicating whether the shortest paths "in" or "out" of the nodes in the directed graphs should be considered. For undirected graphs we use "all".
#' @param weights Numeric vector indicating weights of the edges
#'
#' @return
#' a vector including centrality values for each node
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Scardoni, G. and Carlo Laudanna, C.B.M.C., 2011. Network centralities for Cytoscape. University of Verona.
#'
#' @examples
#'
#' data(zachary)
#'
#' wiener_index_centrality(zachary)
#'
#' @export
#' @importFrom igraph shortest.paths
#' @importFrom igraph V
#' @importFrom igraph is_igraph
#' @importFrom igraph is_named
#' @importFrom igraph getIgraphOpt
#' @importFrom intergraph asIgraph

wiener_index_centrality<-function (x, vids = V(x), mode = c("all", "out", "in"), weights = NULL){

  if (!("igraph" %in% class(x)||"network" %in% class(x))) stop("The input is not an igraph or a
                                                            network object")
  if (is_igraph(x)){

    distMat<- shortest.paths(x , mode = mode[1], weights = weights)
    diag(distMat)<-0

    res <- rowSums(distMat) - diag(distMat)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(res) <- V(x)$name[vids]

    }

  }

  else{

    x<-asIgraph(x)

    distMat<- shortest.paths(x , mode = mode[1], weights = weights)
    diag(distMat)<-0

    res <- rowSums(distMat) - diag(distMat)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(res) <- V(x)$name[vids]

    }

  }
  return(res=res)
}


#' @title Harmonic Centrality
#'
#' @description This function computes Harmonic Centrality.
#' The harmonic metric defines as the denormalized reciprocal of the harmonic mean of all distances.
#'
#' @param x An igraph or a network object
#' @param vids Nodes to be considered in the calculation
#' @param mode a character value, “out” for out-degree, “in” for in-degree or “total” for the sum of the two. For undirected graphs this argument is ignored. “all” is a synonym of “total”.
#' @param weights  Numeric vector indicating weights of the edges
#'
#' @return
#' a vector including centrality values for each node
#'
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#'
#' @references
#' BOLDI, P. & VIGNA, S. 2014. Axioms for centrality. Internet Mathematics, 00-00.
#'
#' MARCHIORI, M. & LATORA, V. 2000. Harmony in the small-world. Physica A: Statistical Mechanics and its Applications, 285, 539-546.
#'
#' OPSAHL, T., AGNEESSENS, F. & SKVORETZ, J. 2010. Node centrality in weighted networks: Generalizing degree and shortest paths. Social Networks, 32, 245-251.
#'
#' OPSAHL, T. 2010. Closeness centrality in networks with disconnected components (http://toreopsahl.com/2010/03/20/closeness-centrality-in-networks-with-disconnected-components/)
#'
#' @examples
#'
#' data(zachary)
#'
#' harmonic_centrality(zachary)
#'
#' @export
#' @importFrom igraph distances
#' @importFrom igraph vcount
#' @importFrom intergraph asIgraph

harmonic_centrality <- function (x, vids = V(x), mode = c("all", "out", "in"), weights = NULL){

  if (!("igraph" %in% class(x) || "network" %in% class(x))) stop("The input is not an igraph or a
                                                            network object")
  if (is_igraph(x)){

    distMat<-1/distances(x, mode = mode[1], weights = weights)
    diag(distMat)<-0

    res <- rowSums(distMat) - diag(distMat)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(res) <- V(x)$name[vids]

    }

  }
  else{

    x<-asIgraph(x)

    distMat<-1/distances(x, mode = mode[1], weights = weights)
    diag(distMat)<-0

    res <- rowSums(distMat) - diag(distMat)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(res) <- V(x)$name[vids]

    }
  }

  return(res)
}

#' @title Group Centrality
#'
#' @description This function computes group Centrality. So, it considers a consistent ranking of
#' each node to be calculated such that scores diverse possible synergies among possible groups of vertices.
#'
#' @param x An igraph or a network object
#' @param vids Nodes to be considered in the calculation
#'
#' @return
#' a vector including centrality values for each node
#'
#' @author Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
#' @references
#' Michalak, T.P., Aadithya, K.V., Szczepanski, P.L., Ravindran, B. and Jennings, N.R., 2013. Efficient computation of the Shapley value for game-theoretic network centrality. Journal of Artificial Intelligence Research, 46, pp.607-650.
#'
#' https://www.civilica.com/Paper-IBIS07-IBIS07_127.html
#'
#' @examples
#'
#' data(zachary)
#'
#' group_centrality(zachary)
#'
#' @export
#' @importFrom igraph neighbors
#' @importFrom igraph degree
#' @importFrom igraph V
#' @importFrom igraph is_igraph
#' @importFrom igraph is_named
#' @importFrom igraph getIgraphOpt
#' @importFrom intergraph asIgraph

group_centrality<-function (x, vids = V(x)){

  if (!("igraph" %in% class(x)|| "network" %in% class(x))) stop("The input is not an igraph or a
                                                            network object")

  f <- function(v){ 1/(1 + degree(x,v)) + sum(1/(1+degree(x,neighbors(x,v,'all'))))}

  if (is_igraph(x)){

    results <- sapply(V(x),f)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(results) <- V(x)$name[vids]
    }
  }

  else{

    x<-asIgraph(x)

    results <- sapply(V(x),f)

    if (getIgraphOpt("add.vertex.names") && is_named(x)) {
      names(results) <- V(x)$name[vids]
    }
  }
  return(results)
}


