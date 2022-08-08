"rrp.predict" <-
function(x, y, train, test, k=1){

if( !any(class(x) == "XPtr") )  
 stop("`x' must be a of class `XPtr'")

 if(attr(x, "Size") != (length(train)+length(test)))
  stop("Wrong dimensions")

 if(!is.numeric(y))
  stop("`y' must be numeric")

 pred <- as.numeric(rep(NA,length(test)))

 f <- function(x) {x[which(x==1)] <- NA; tmp <- order(x,na.last=NA); tmp[1:min(k,length(tmp))]}
 nn <- applyXPtr(x, test, train, f)
 
 g <- function(x) {x[which(x==1)] <- NA; tmp <- order(x,na.last=NA); 
   tmp <- tmp[1:min(k,length(tmp))]; x[tmp]}
 ww <- applyXPtr(x, test, train, g)
 
 for(i in 1:length(test)){
  tmp <- nn[[i]]
  wt <- ww[[i]]
  if(length(tmp)>0)
   pred[i] <- weighted.mean(y[tmp], w=1-wt, na.rm=TRUE)
 }
 return(pred)
}

