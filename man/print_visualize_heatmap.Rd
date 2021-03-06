% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CINNA.R
\name{print_visualize_heatmap}
\alias{print_visualize_heatmap}
\title{Print the heatmap plot of centrality measures}
\usage{
print_visualize_heatmap(x, scale = TRUE, file = NULL)
}
\arguments{
\item{x}{a list indicating calculated centrality measures}

\item{scale}{Whether the centrality values should be scaled or not(default=TRUE)}

\item{file}{A character string naming the .pdf file to print into. If NULL
the result would be printed to the exist directory.(default=NULL)}
}
\value{
The resulted plot of \code{ \link[CINNA]{visualize_heatmap}}function will be saved in the given directory.
}
\description{
This function prints the heatmap plot
}
\author{
Minoo Ashtiani, Mehdi Mirzaie, Mohieddin Jafari
}
