"applyXPtr" <-
function(d, idx, sub, f)
{
if(!is.integer(idx))
  x <- as.integer(idx)
if(!is.integer(sub))
  sub <- as.integer(sub)

if( !any(class(d) == "XPtr") )  
 stop("`d' must be a of class `XPtr'")

if(any(idx>attr(d,"Size")))
 stop("indexes in `idx' behond the size of `d'")

if(any(sub>attr(d,"Size")))
 stop("indexes in `sub' behond the size of `d'")

 return(invisible(.Call("applyXPtr", d, idx, sub, f, .GlobalEnv, PACKAGE="rrp")))
}
