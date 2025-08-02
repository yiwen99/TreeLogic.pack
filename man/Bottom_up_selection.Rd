\name{Bottom_up_selection}
\alias{Bottom_up_selection}
\title{Bottom_up_selection}
\usage{
Bottom_up_selection(y, X, cvrt = NULL, tree, type='default')
}
\arguments{
  \item{y}{the response vector}
  \item{X}{leaf-level microbiome P/A status matrix }
  \item{cvrt}{a non-microbiome covariate matrix }
  \item{tree}{the tree structure that we use, dimension needs to match with number of columns of X}
  \item{type}{'default' or 'global min'}
}
\description{
This is the function for Bottom-Up Feature Aggregation. This function include two types: "default" (early stop) and "global min".

The "global min" type combines features from the leaf up until the root, and returns the trimmed tree at the point where the aggregation process yields globally minimal BIC value.

The "default" type, stops the aggregation process at the point right before BIC starts to increase.


The function returns

(1)final.tree: a data frame of the final trimmed tree structure, represented with node labels


(2) final.X: a matrix of the (aggregated) P/A status of the leaf nodes in the trimmed tree. Feature names selected by the Bottom-up process can be accessed by colnames(result_list$final.X)


(3) BIC.seq: the BIC trajectory throughout the feature aggregation process
}
\author{
Yiwen Chen
}
\examples{
#please see vignettes
}
