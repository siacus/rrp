"XPtrToDist" <-
function(d)
{
if( !any(class(d) == "XPtr") )  
 stop("`d' must be a of class `XPtr'")
 tmp <- .Call("XPtrToNumeric", d, PACKAGE="rrp")
 attributes(tmp) <- list(Size=attr(d,"Size"), class="dist")
 return(invisible(tmp))
}
