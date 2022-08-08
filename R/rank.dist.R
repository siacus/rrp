rank.dist <- function (X, msplit = 2, cut.in = 0, 
   thr=0.75, weights, asdist=FALSE, verbose = 0) 
{
 n <- dim(X)[1]
 nv <- dim(X)[2]
 if(missing(weights))
  weights<- rep(1/nv,nv)
 if(any(weights<0))
   stop("weights must be positive")
 if(length(weights) != nv){
  stop("weights length must equal number of variables")
 }
 if(sum(weights)>1){
   weights <- weights/sum(weights)
   warning("weights normalized to 1")
   }  

  fv <- NULL

  for(i in 1:nv){
   if(cut.in>0){
    if (is.numeric(X[[i]]) | is.integer(X[[i]]))
     X[[i]] <- cut(X[[i]], seq(min(X[[i]], 
                    na.rm = TRUE), max(X[[i]], na.rm = TRUE), 
                    length = cut.in), include.lowest = TRUE,labels=FALSE)
   }
	if(is.logical(X[[i]]))
	 X[[i]] <- factor(X[[i]])
	if(is.factor(X[[i]]))
	  fv <- c(fv,i)	 		
  }

 px <- vector(n, mode="list")

 for(i in 1:n){
  if(verbose>0) cat(".")
  px[[i]] <- numeric(n)
  for(v in fv){
    idx <- which(X[[v]] == X[i,v])
	tmp <- rep(weights[v], n)
	tmp[-idx] <- 0
	px[[i]] <- px[[i]] + tmp
  }

  for(v in (1:nv)[-fv]){
   tmp <- abs(X[i,v]-X[[v]])/(msplit-1)
   idx <- which(tmp > 1)
   tmp[idx] <- 1
   px[[i]] <- px[[i]] + (1- tmp)*weights[v]
  }

  px[[i]] <- px[[i]]/sum(weights)
  idx <- which(px[[i]]>thr)
  px[[i]] <- (px[[i]])[idx]
  names(px[[i]]) <- idx
  px[[i]] <- sort(px[[i]], decreasing=TRUE)
 }
 if(asdist){
  mat <- matrix(0,n,n)
  for(i in 1:n){
   idx <- as.integer(names(px[[i]]))
   mat[i,idx] <- as.numeric(px[[i]])
  }
  return(invisible(as.dist(mat)))
 } else 
   return(invisible(px))
}
