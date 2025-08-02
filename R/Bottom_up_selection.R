#This is the function for Bottom-Up Feature Aggregation. This function include two types: "default" (early stop) and "global min".
#The "global min" type combines features from the leaf up until the root, and returns the trimmed tree at the point where the aggregation process yields globally minimal BIC value.
#The "default" type, also being referred as "early stop" type, stops the aggregation process at the point right before BIC starts to increase.
#The input parameters of the function are  
                                          #y: the response vector 
                                          #X: leaf-level microbiome P/A status matrix 
                                          #cvrt: a non-microbiome covariate matrix 
                                          #tree: the tree structure that we use, dimension needs to match with number of columns of X
                                          #type: 'default' or 'global min'
#The function returns 
                      #final.tree: a data frame of the final trimmed tree structure, represented with node labels
                      #final.X: a matrix of the (aggregated) P/A status of the leaf nodes in the trimmed tree. Feature names selected by the Bottom-up process can be accessed by colnames(result_list$final.X)
                      #BIC.seq: the BIC trajectory throughout the feature aggregation process
#Additionally, note that the function uses a function BIC.calc(), also included in the package, to calculate BIC of a linear model 


Bottom_up_selection <- function(y, X, cvrt = NULL, tree, type='default'){
  if(type=='default'){
    n1 = length(y)
    n2 = dim(X)[1]
    p1 = dim(X)[2]
    p2 = dim(tree)[1]
    
    
    if(is.vector(cvrt)){q=1
    } else if(is.null(cvrt)){q=0
    }  else {q=dim(cvrt)[2]}
    
    
    if(n1!=n2) stop('sample size does not match')
    if(p1!=p2) stop('dimension does not match')
    if(p1+q>n1) stop('sample size too small for the full model')
    
    # some basic values
    n=n1
    p=p1
    q=q
    L=dim(tree)[2] # number of tree levels
    BIC.prev=Inf
    BIC.curr=BIC.calc(y,cbind(cvrt,X))
    BIC.seq=BIC.curr
    tree.curr=tree
    # assign node labels/Simplify tree (as initial setup of node labels)
    node.label=matrix(3,nrow=nrow(tree),ncol=ncol(tree)) # true leaf=1, immediate parent=2, other internal=3, redundant leaf=-1
    node.label[,1]=1
    for(level in 2:L){
      for(ind in 1:length(unique(tree[,level]))){
        curr.node=unique(tree[,level])[ind]
        loc=which(tree[,level]==curr.node)
        if(length(loc)==1){
          node.label[loc,level]=-1
          if(node.label[loc,level-1]!=-1){node.label[loc,level-1] = 1}
        }else if(prod(node.label[loc,level-1])==1){
          node.label[loc,level]=2
        }else{node.label[loc,level]=3}
      }
    }
    
    # get current X (not matched with tree)
    X.curr=matrix(NA,nrow=nrow(X),ncol=ncol(X))
    leafs=unique(tree[node.label==1])
    for(j in 1:length(leafs)){
      loc=which(tree==leafs[j],arr.ind = TRUE)
      if(length(loc[1])>1)stop("leaf should only be one node")
      X.curr[,j]=X[,loc[1]]
    }
    colnames(X.curr)=leafs
    
    BIC.rec=BIC.curr
    while(BIC.curr < BIC.prev){
      BIC.prev = BIC.curr
      
      # combine leaf nodes for a parent node
      BIC.temp=c()
      eligible.parent=unique(tree[node.label==2])
      if(length(eligible.parent)>0){
        for(j in 1:length(eligible.parent)){
          parent=eligible.parent[j]
          loc=which(tree==parent,arr.ind = TRUE) # find all its leafs in the original tree
          # first, create new pred from original data
          X.comb=(rowSums(X[,loc[,1]])>0)+0
          # second, remove all children of parent from X.curr
          blklist=unlist(tree[loc[,1],(1:loc[1,2]-1)]) # should be a vector of taxa names
          mch=match(colnames(X.curr),blklist)
          X.clean=X.curr[,is.na(mch)]
          # third, add the new pred (and its higher level taxa name)
          X.temp=cbind(X.clean,X.comb)
          colnames(X.temp)[ncol(X.temp)]=parent
          
          BIC.temp[j]=BIC.calc(y,cbind(cvrt,X.temp))
        }
        sel=which.min(BIC.temp) # selected parent node to be combined
        min.BIC=BIC.temp[sel] # smallest BIC
      }
      
      if(min.BIC >= BIC.prev){break}
      
      
      parent.cb=eligible.parent[sel]
      print(paste("All children of",parent.cb,"are combined."))
      # 1. modify X.curr
      loc=which(tree==parent.cb,arr.ind = TRUE) 
      X.comb=(rowSums(X[,loc[,1]])>0)+0
      blklist=unlist(tree[loc[,1],(1:loc[1,2]-1)]) # should be a vector of taxa names
      mch=match(colnames(X.curr),blklist)
      X.clean=X.curr[,is.na(mch)]
      X.curr=cbind(X.clean,X.comb)
      colnames(X.curr)[ncol(X.curr)]=parent.cb
      # 2. immediate change to node.label
      node.label[loc]=1
      node.label[loc[,1],(1:loc[,2][1]-1)] = ifelse(node.label[loc[,1],(1:loc[,2][1]-1)]==1,-1,node.label[loc[,1],(1:loc[,2][1]-1)])
      
      BIC.curr=min.BIC
      
      for(level in 2:(L-1)){
        checkset=unique(tree[which(node.label[,level]>=1),level])
        if(length(checkset)==0){next}
        for(ind in 1:length(checkset)){
          curr.node=checkset[ind]
          loc=which(tree[,level]==curr.node) # row indices of which row of tree matrix is under this curr node
          if(node.label[loc[1],level]==1 & sum(tree[,level+1]==tree[loc[1],level+1])==length(loc)){
            node.label[loc,level]=-1
            node.label[loc,level+1]=1} # 1->-1
          if(node.label[loc[1],level]==2 & all(node.label[loc, 1:level]!=1)){node.label[loc,level]=1} # 2->1
          if(node.label[loc[1],level]==3 & all(node.label[loc, 1:level]==-1)){ #3->1
            node.label[loc,level]=1
          } else if(node.label[loc[1],level]==3 & all(node.label[loc,1:(level-1)]<=1)){
            node.label[loc,level]=2
          } # 3->2
        }
      }
      
      BIC.rec=c(BIC.rec,BIC.curr)
    }
    dim(X.curr)
    return(list(final.tree=node.label,final.X=X.curr, BIC.seq=BIC.rec))
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
    if(p1!=p2) stop('dimension does not match')
    if(p1+q>n1) stop('sample size too small for the full model')
    
    # some basic values
    n=n1
    p=p1
    q=q
    L=dim(tree)[2] # number of tree levels
    BIC.prev=Inf
    BIC.curr=BIC.calc(y,cbind(cvrt,X))
    BIC.seq=BIC.curr
    tree.curr=tree
    # assign node labels/Simplify tree (OK, as initial setup of node labels)
    node.label=matrix(3,nrow=nrow(tree),ncol=ncol(tree)) # true leaf=1, immediate parent=2, other internal=3,     redundant leaf=-1
    node.label[,1]=1
    for(level in 2:L){
      for(ind in 1:length(unique(tree[,level]))){
        curr.node=unique(tree[,level])[ind]
        loc=which(tree[,level]==curr.node)
        if(length(loc)==1){
          node.label[loc,level]=-1
          if(node.label[loc,level-1]!=-1){node.label[loc,level-1] = 1}
        }else if(prod(node.label[loc,level-1])==1){
          node.label[loc,level]=2
        }else{node.label[loc,level]=3}
      }
    }
    
    # get current X (not matched with tree)
    X.curr=matrix(NA,nrow=nrow(X),ncol=ncol(X))
    leafs=unique(tree[node.label==1])
    for(j in 1:length(leafs)){
      loc=which(tree==leafs[j],arr.ind = TRUE)
      if(length(loc[1])>1)stop("leaf should only be one node")
      X.curr[,j]=X[,loc[1]]
    }
    colnames(X.curr)=leafs
    combined.X <- X.curr
    BIC.rec=BIC.curr
    
    eligible.parent <- -99
    while(length(eligible.parent)>0){
      BIC.prev = BIC.curr
      
      # combine leaf nodes for a parent node
      BIC.temp=c()
      eligible.parent=unique(tree[node.label==2])
      if(length(eligible.parent)>0){
        for(j in 1:length(eligible.parent)){
          parent=eligible.parent[j]
          loc=which(tree==parent,arr.ind = TRUE) # find all its leafs in the original tree
          # first, create new pred from original data
          X.comb=(rowSums(X[,loc[,1]])>0)+0
          # second, remove all children of parent from X.curr
          blklist=unlist(tree[loc[,1],(1:loc[1,2]-1)]) # should be a vector of taxa names
          mch=match(colnames(X.curr),blklist)
          X.clean=X.curr[,is.na(mch)]
          # third, add the new pred (and its higher level taxa name)
          X.temp=cbind(X.clean,X.comb)
          colnames(X.temp)[ncol(X.temp)]=parent
          
          BIC.temp[j]=BIC.calc(y,cbind(cvrt,X.temp))
        }
        sel=which.min(BIC.temp) # selected parent node to be combined
        min.BIC=BIC.temp[sel] # smallest BIC
      }
      
      if(length(eligible.parent)==0){break}
      
      
      parent.cb=eligible.parent[sel]
      
      print(paste("All children of",parent.cb,"are combined."))
      # 1. modify X.curr
      loc=which(tree==parent.cb,arr.ind = TRUE) 
      X.comb=(rowSums(X[,loc[,1]])>0)+0
      blklist=unlist(tree[loc[,1],(1:loc[1,2]-1)]) # should be a vector of taxa names
      mch=match(colnames(X.curr),blklist)
      X.clean=X.curr[,is.na(mch)]
      X.curr=cbind(X.clean,X.comb)
      colnames(X.curr)[ncol(X.curr)]=parent.cb
      # 2. immediate change to node.label

      node.label[loc]=1
      node.label[loc[,1],(1:loc[,2][1]-1)] = ifelse(node.label[loc[,1],(1:loc[,2][1]-1)]==1,-1,node.label[loc[,1],(1:loc[,2][1]-1)])
      
      BIC.curr=min.BIC
      
      for(level in 2:(L-1)){
        checkset=unique(tree[which(node.label[,level]>=1),level])
        if(length(checkset)==0){next}
        for(ind in 1:length(checkset)){
          curr.node=checkset[ind]
          loc=which(tree[,level]==curr.node) # row indices of which row of tree matrix is under this curr node
          if(node.label[loc[1],level]==1 & sum(tree[,level+1]==tree[loc[1],level+1])==length(loc)){
            node.label[loc,level]=-1
            node.label[loc,level+1]=1} # 1->-1
          if(node.label[loc[1],level]==2 & all(node.label[loc, 1:level]!=1)){node.label[loc,level]=1} # 2->1
          if(node.label[loc[1],level]==3 & all(node.label[loc, 1:level]==-1)){ #3->1
            node.label[loc,level]=1
          } else if(node.label[loc[1],level]==3 & all(node.label[loc,1:(level-1)]<=1)){
            node.label[loc,level]=2
          } # 3->2
        }
      }
      
      BIC.rec=c(BIC.rec,BIC.curr)
      if(BIC.rec[length(BIC.rec)]<= min(BIC.rec[1:length(BIC.rec)-1])){
        combined.X <- X.curr
      }
    }

    return(list(final.tree=node.label, final.X = combined.X, BIC.seq=BIC.rec)) #combined.X is the point where it reaches global min for BIC
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