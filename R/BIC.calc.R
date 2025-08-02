
#This function is a helper function that calculates the BIC
#The function takes in
                   #y: a response vector
                   #pred: the design matrix
#The function returns
                   #BIC value of the fitted linear model
#' @export
BIC.calc<-function(y, pred){
  # fit a linear regression of y~pred and calc BIC
  out=lm(y~as.matrix(pred))
  value=BIC(out)
  return(value)
}
