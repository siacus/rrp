"setXPtr" <-
function(d, x, k)
{
if(!is.list(x))
  x <- as.list(x)
if( !any(class(d) == "XPtr") )  
 stop("`d' must be a of class `XPtr'")

if(any(unlist(x)>attr(d,"Size")))
 stop("indexes in `x' behond the size of `d'")

 n <- length(x)
 k <- rep(k,n)[1:n]

 x <- lapply(x, sort)

 return(invisible(.Call("setXPtr", d, x, as.numeric(k), PACKAGE="rrp")))
}
