#' @title Moreno Rhesus Network
#'
#' @description A directed graph including observed grooming episodes between
#' free ranging rhesus macaques (Macaca mulatta) in Cayo Santiago
#' during a two month period in 1963. Cayo Santiago is an island off the coast of Puerto Rico,
#' which also is named as Isla de los monos (Island of the monkeys).
#' A node indicates a monkey and a directed edge in which
#' a rhesus macaque groomed another rhesus macaque.
#' The weights of edges demonstrates how often this behaviour
#' was seen.
#' @name rhesus
#' @docType data
#' @usage rhesus
#' @format an igraph object with "gml" format
#' @references Rhesus network dataset -- KONECT, October 2016.
#'
#' DS Sade. Sociometrics of macaca mulatta I. linkages and cliques in grooming matrices. Folia Primatologica, 18(3-4):196--223, 1972.
#' @keywords datasets
#' @examples
#' data("rhesus")
#' print(rhesus)
NULL
