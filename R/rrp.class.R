"rrp.class" <-
function(x, cl, train, test, k=1){

if( !any(class(x) == "XPtr") )  
 stop("`x' must be a of class `XPtr'")

 if(attr(x, "Size") != (length(train)+length(test)))
  stop("Wrong dimensions")

 if(!is.factor(cl))
  cl <- factor(cl)

 pred <- factor(rep(NA,length(test)), levels=levels(cl))
 
 f <- function(x) {x[which(x==1)] <- NA; tmp <- order(x,na.last=NA); tmp[1:min(k,length(tmp))]}
 nn <- applyXPtr(x, test, train, f)
 
 for(i in 1:length(test)){
  tmp <- nn[[i]]
 
  if(length(tmp)>1){
   votes <- table(cl[tmp])
   idx <- which(votes == max(votes))   
   if(length(idx)>1)
     idx <- sample(idx,1)
	 
    pred.cl <- names(votes)[idx]
    pred[i] <- pred.cl
  } else {
    if(length(tmp)>0)
     pred[i] <- cl[tmp]
  }  
 }
 return(pred)
}

