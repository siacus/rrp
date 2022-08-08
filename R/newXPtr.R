"newXPtr" <-
function(n, k=0)
{
 tmp <- .Call("newXPtr", as.integer(n), as.numeric(k), PACKAGE = "rrp")
 attributes(tmp) <- list(class=c(class(tmp), "XPtr"), Size=n)
 return(invisible(tmp))
}

