#This is the function for Top-Down Feature Selection. This function include two types: "default" (early stop) and "global min".
#The "global min" type spans features from the root until the leaf level, and returns the "spanned" tree at the point where the spanning process yields globally minimal BIC value.
#The "default" type, also being referred as "early stop" type, stops the aggregation process at the point right before BIC starts to increase.
#The input parameters of the function are
                                          #y: the response vector
                                          #X: the matrix including the infered P/A status (through "OR" aggregation) of all nodes in the used tree structure, note that the Top_down_selection function requires X to have colnames() attributes
                                          #cvrt: a non-microbiome covariate matrix
                                          #tree: the tree structure that we use, need to contain all nodes appearing in X matrix
                                          #type: 'default' or 'global min'
#The function returns
                                          #final.tree: a data frame of the final spanned tree structure, represented with node labels
                                          #final.X: a matrix of the (aggregated) P/A status of the leaf nodes in the spanned tree. Feature names selected by the Top-down process can be accessed by colnames(result_list$final.X)
                                          #BIC.seq: the BIC trajectory throughout the feature aggregation process
#Additionally, note that the function uses a function BIC.calc(), also included in the package, to calculate BIC of a linear model
#' @export
Top_down_selection <- function(y, X, cvrt = NULL, tree, type='default'){
  if(type=='default'){
    n1 = length(y)
    n2 = dim(X)[1]
    p1 = dim(X)[2]
    p2 = dim(tree)[1]


    if(is.vector(cvrt)){q=1
    } else if(is.null(cvrt)){q=0
    }  else {q=dim(cvrt)[2]}


    if(n1!=n2) stop('sample size does not match')

    # some basic values
    n=n1
    p=p1
    q=q
    L=dim(tree)[2] # number of tree levels

    #assign node labels (desired: (the immediate candidate parent to span)true leaf(cannot be spanned)->-1, eligible parent->1, spanned parent and unexplored node ->3, root -> 0)
    #stop if there is no eligible parent to span
    node.label=matrix(3,nrow=nrow(tree),ncol=ncol(tree))
    node.label[,L] <- 0
    node.label[,(L-1)] <- 1

    eligible.parent <- unique(tree[node.label==1])

    leaves <- unique(tree[node.label==-1])
    X.curr=matrix(NA,nrow=nrow(X),length(leaves)+length(eligible.parent))

    X.curr <- X[,colnames(X)%in%c(eligible.parent, leaves)]
    returned.X <- X.curr

    #compute current BIC
    BIC.prev=-Inf
    BIC.curr=BIC.calc(y,cbind(cvrt,X.curr))
    BIC.seq=BIC.curr


    while(length(eligible.parent)>0){
      BIC.prev <- BIC.curr
      X.prev <- X.curr

      BIC.temp=c()
      for(i in 1:length(eligible.parent)){
        #span each of the eligible parents
        parent = eligible.parent[i]
        loc=which(tree==parent,arr.ind = TRUE)
        X.span.list <- unique(tree[cbind(loc[,1],loc[,2]-1)])
        X.span <- matrix(X[,match(X.span.list,colnames(X))],nrow=n1)
        colnames(X.span) <- X.span.list
        X.temp <- cbind(X.span, X.curr[,-which(colnames(X.curr)==parent)])

        #calculate temp BIC
        BIC.temp[i] <- BIC.calc(y,cbind(cvrt,X.temp))
      }

      #choose minimum of BIC.temp, get the index, which is the index of the parent we want to span (first combine redundant nodes (who do not cause any change in bic, for example, one phylum-one class))
      min.BIC <- min(BIC.temp)
      if(min.BIC>BIC.prev){
        break
      }

      #which parent to span
      span.parent.idx <- which.min(BIC.temp)
      span.parent <- eligible.parent[span.parent.idx]
      print(paste0(span.parent," is spanned"))
      loc=which(tree==span.parent,arr.ind = TRUE)
      X.span.list <- unique(tree[cbind(loc[,1],loc[,2]-1)])
      X.span <- matrix(X[,match(X.span.list,colnames(X))],nrow=n1)
      colnames(X.span) <- X.span.list
      if(is.null(dim(X.curr[,-which(colnames(X.curr)==span.parent)]))){
        X.curr_new <- cbind(X.span, X.curr[,-which(colnames(X.curr)==span.parent)])
        colnames(X.curr_new) <- c(colnames(X.span), colnames(X.curr)[-which(colnames(X.curr)==span.parent)])
        X.curr <- X.curr_new
      }else{
        X.curr <- cbind(X.span, X.curr[,-which(colnames(X.curr)==span.parent)])
      }

      BIC.curr <- min.BIC
      BIC.seq <- c(BIC.seq, BIC.curr)


      #reassign node labels
      #now X.curr: check if they are true leaves, if not they are eligible parents.
      node.label=matrix(3,nrow=nrow(tree),ncol=ncol(tree))
      node.label[,L] <- 0
      for(j in 1:length(colnames(X.curr))){
        curr.loc <- which(tree==colnames(X.curr)[j],arr.ind=TRUE)
        if(curr.loc[1,2]==1){
          node.label[curr.loc] <- -1
        }else{
          node.label[curr.loc] <- 1
        }
      }

      #reassign eligible parent and leaves (leaves: ones that do not have any children in the full tree)
      eligible.parent <- unique(tree[node.label==1])
      leaves <- unique(tree[node.label==-1])

    }

    return(list(final.tree=node.label, final.X = X.curr, BIC.seq = BIC.seq))
  }
  if(type=='global min'){

    n1 = length(y)
    n2 = dim(X)[1]
    p1 = dim(X)[2]
    p2 = dim(tree)[1]


    if(is.vector(cvrt)){q=1
    } else if(is.null(cvrt)){q=0
    }  else {q=dim(cvrt)[2]}


    if(n1!=n2) stop('sample size does not match')

    # some basic values
    n=n1
    p=p1
    q=q
    L=dim(tree)[2] # number of tree levels

    #assign node labels (desired: (the immediate candidate parent to span)true leaf(cannot be spanned)->-1, eligible parent->1, spanned parent and unexplored node ->3, root -> 0)
    #stop if there is no eligible parent to span
    node.label=matrix(3,nrow=nrow(tree),ncol=ncol(tree))
    node.label[,L] <- 0
    node.label[,(L-1)] <- 1

    #redundant nodes label -99

    eligible.parent <- unique(tree[node.label==1])
    leaves <- unique(tree[node.label==-1])
    X.curr=matrix(NA,nrow=nrow(X),length(leaves)+length(eligible.parent))
    X.curr <- X[,colnames(X)%in%c(eligible.parent, leaves)]
    returned.X <- X.curr

    #compute current BIC
    BIC.prev=-Inf
    BIC.curr=BIC.calc(y,cbind(cvrt,X.curr))
    BIC.seq=BIC.curr


    while(length(eligible.parent)>0){
      BIC.prev <- BIC.curr
      X.prev <- X.curr

      BIC.temp=c()
      for(i in 1:length(eligible.parent)){
        #span each of the eligible parents
        parent = eligible.parent[i]
        loc=which(tree==parent,arr.ind = TRUE)
        X.span.list <- unique(tree[cbind(loc[,1],loc[,2]-1)])
        X.span <- matrix(X[,match(X.span.list,colnames(X))],nrow=n1)
        colnames(X.span) <- X.span.list
        X.temp <- cbind(X.span, X.curr[,-which(colnames(X.curr)==parent)])

        #calculate temp BIC
        BIC.temp[i] <- BIC.calc(y,cbind(cvrt,X.temp))
      }

      #choose minimum of BIC.temp, get the index, which is the index of the parent we want to span (first combine redundant nodes (who do not cause any change in bic, for example, one phylum-one class))
      min.BIC <- min(BIC.temp)
      #which parent to span
      span.parent.idx <- which.min(BIC.temp)
      span.parent <- eligible.parent[span.parent.idx]
      print(paste0(span.parent," is spanned"))
      loc=which(tree==span.parent,arr.ind = TRUE)
      X.span.list <- unique(tree[cbind(loc[,1],loc[,2]-1)])
      X.span <- matrix(X[,match(X.span.list,colnames(X))],nrow=n1)
      colnames(X.span) <- X.span.list
      #X.curr <- cbind(X.span, X.curr[,-which(colnames(X.curr)==span.parent)])
      if(is.null(dim(X.curr[,-which(colnames(X.curr)==span.parent)]))){
        X.curr_new <- cbind(X.span, X.curr[,-which(colnames(X.curr)==span.parent)])
        colnames(X.curr_new) <- c(colnames(X.span), colnames(X.curr)[-which(colnames(X.curr)==span.parent)])
        X.curr <- X.curr_new
      }else{
        X.curr <- cbind(X.span, X.curr[,-which(colnames(X.curr)==span.parent)])
      }

      BIC.curr <- min.BIC
      BIC.seq <- c(BIC.seq, BIC.curr)


      #reassign node labels
      #now X.curr: check if they are true leaves, if not they are eligible parents.
      node.label=matrix(3,nrow=nrow(tree),ncol=ncol(tree))
      node.label[,L] <- 0
      #use a for loop
      for(j in 1:length(colnames(X.curr))){
        curr.loc <- which(tree==colnames(X.curr)[j],arr.ind=TRUE)
        if(curr.loc[1,2]==1){
          node.label[curr.loc] <- -1
        }else{
          node.label[curr.loc] <- 1
        }
      }

      #reassign eligible parent and leaves (leaves: ones that do not have any children in the full tree)
      eligible.parent <- unique(tree[node.label==1])
      leaves <- unique(tree[node.label==-1])

      #if BIC.curr[length(BIC.curr)] <= min(BIC.curr[length(BIC.curr)-1]), update X.curr
      if(dim(X.curr)[2]>=dim(X.curr)[1]){
        return(list(final.tree=node.label,trimmed.X.end=X.prev, bic.seq = BIC.seq[1:(length(BIC.seq)-1)], final.X = returned.X))
      }
      if(BIC.seq[length(BIC.seq)]<= min(BIC.seq[1:length(BIC.seq)-1])){
        returned.X <- X.curr
      }
    }

    return(list(final.tree=node.label, final.X = returned.X, BIC.seq = BIC.seq))
  }
  else{
    print("wrong type")
  }

}

# BIC.calc<-function(y, pred){
#   # fit a linear regression of y~pred and calc BIC
#   out=lm(y~as.matrix(pred))
#   value=BIC(out)
#   return(value)
# }
