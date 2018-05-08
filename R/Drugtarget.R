#' @title Drug Target Network
#'
#' @description A bipartite graph extracted from DrugBank 1.0 database. The network includes two set of nodes
#' including Food and Drug Administration (FDA)-approved drugs and their corresponding protein targets designated by their Uniprot ID.
#' The 1080 drugs and their 519 target proteins nodes are connected via 3766 interactions.
#' Please note that it is a shrunken network in which metabolizing enzymes,
#' carriers and transporters associated with drug metabolism are filtered and solely
#' targets directly related to their pharmacological effects are included.
#' It is also an example of unconnected graphs.
#' @name drugTarget
#' @docType data
#' @usage drugTarget
#' @format an igraph object with "gml" format
#' @references Barneh, F., Jafari, M., & Mirzaie, M. (2015). Updates on drug–target network; facilitating polypharmacology and data integration by growth of DrugBank database. Briefings in Bioinformatics, bbv094. https://doi.org/10.1093/bib/bbv094
#' @keywords datasets
#' @examples
#' data("drugTarget")
#' print(drugTarget)
NULL
