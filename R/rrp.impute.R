"rrp.impute" <-
function(data, D = NULL, k=1,  msplit=10, Rep=250, cut.in=15){

 n <- dim(data)[1]
 all <- 1:n
 if(is.null(D) | (!any(class(D) == "XPtr")))
	d <- rrp.dist(data,  msplit=msplit, Rep=Rep,cut.in=cut.in, check.bal=FALSE, 
		plot=FALSE)
 else
  d <- D

 miss.obs <- which(apply(data, 1, function(x) length(which(is.na(x)))>0 ) == TRUE)

 y <- data
 comp.obs <- all[-miss.obs]

 f <- function(x) {x[which(x==1)] <- NA; tmp <- order(x,na.last=NA); tmp[1:min(k,length(tmp))]}
 nn <- applyXPtr(d, miss.obs, comp.obs, f)

 g <- function(x) {x[which(x==1)] <- NA; tmp <- order(x,na.last=NA); 
   tmp <- tmp[1:min(k,length(tmp))]; x[tmp]}
 ww <- applyXPtr(d, miss.obs, comp.obs, g)

 for(i in 1:length(miss.obs)){
  tmp <- nn[[i]] 
  wt <- ww[[i]]
  if(length(tmp>0)){
   v.idx <- which(is.na(data[miss.obs[i],]))	
   for(v in v.idx){
    if(!is.factor(data[[v]])){
	  y[miss.obs[i],v] <- weighted.mean(data[comp.obs[tmp],v], w=1-wt, na.rm=TRUE)
	} else {	
	  if(length(tmp)>1){
       votes <- table(data[comp.obs[tmp],v])
       idx <- which(votes == max(votes))   
       if(length(idx)>1)
        idx <- sample(idx,1)
	     y[miss.obs[i],v] <- names(votes)[idx]
      } else {
       y[miss.obs[i],v] <- data[comp.obs[tmp],v]
      }  
   
	}
   } # for
  } # if(length(tmp>0))
 } # for(i in 1:length(miss.obs))
  return(list(new.data = y, dist = d))
}

