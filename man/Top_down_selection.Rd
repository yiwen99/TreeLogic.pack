\name{Top_down_selection}
\alias{Top_down_selection}
\title{Top_down_selection}
\usage{
Top_down_selection(y, X, cvrt = NULL, tree, type='default')
}
\arguments{
  \item{y}{the response vector}
  \item{X}{the matrix including the infered P/A status (through "OR" aggregation) of all nodes in the used tree structure, note that the Top_down_selection function requires X to have colnames() attributes}
  \item{cvrt}{a non-microbiome covariate matrix}
  \item{tree}{the tree structure that we use, need to contain all nodes appearing in X matrix}
  \item{type}{'default' or 'global min'}
}
\description{
This is the function for Top-Down Feature Selection. This function include two types: "default" (early stop) and "global min".


The "global min" type spans features from the root until the leaf level, and returns the "spanned" tree at the point where the spanning process yields globally minimal BIC value.


The "default" type, also being referred as "early stop" type, stops the aggregation process at the point right before BIC starts to increase.


The function returns

(1) final.tree: a data frame of the final spanned tree structure, represented with node labels


(2) final.X: a matrix of the (aggregated) P/A status of the leaf nodes in the spanned tree. Feature names selected by the Top-down process can be accessed by colnames(result_list$final.X)

(3) BIC.seq: the BIC trajectory throughout the feature aggregation process
}
\author{
Yiwen Chen
}
\examples{
#please see vignettes
}
