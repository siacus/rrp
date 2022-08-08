"addDist" <-
function(d, x, k)
{
 if(!is.list(x))
  x <- as.list(x)
 if(!(class(d) == "dist"))  
  stop("`d' must be a of class `dist'")
 if(any(unlist(x)>attr(x,"Size")))
  stop("indexes in `x' behond the size of `d'")
 x <- lapply(x, sort)  
 n <- length(x)
 k <- rep(k,n)[1:n]
 return(invisible(.Call("addDist", d, x, as.numeric(k), PACKAGE = "rrp")))
}

