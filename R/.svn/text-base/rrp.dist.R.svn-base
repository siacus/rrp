
	
RRP.rpart <- function (formula,  data, weights, subset, na.action = na.rpart, 
    method, model = FALSE, x = FALSE, y = TRUE, parms, control, 
    cost, m=NULL, RRPY=NULL, RRPX=NULL, ...) 
{
    call <- match.call()

    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1)) 
        stop("Trees cannot handle interaction terms")
    Y <- RRPY
    wt <- model.extract(m, "weights")
    if (length(wt) == 0) 
        wt <- rep(1, nrow(m))
    offset <- attr(Terms, "offset")

    nobs <- nrow(RRPX)
    nvar <- ncol(RRPX)
  
	   
	
        method.int <- 1 # "anova"
		
        if (missing(parms)) 
            init <- (get(paste("RRP.rpart", method, sep = ".")))(Y, 
                offset, , wt)
        else init <- (get(paste("RRP.rpart", method, sep = ".")))(Y, 
            offset, parms, wt)
    
	Y <- init$y
    xlevels <- attr(RRPX, "column.levels")
    cats <- rep(0, ncol(RRPX))
    if (!is.null(xlevels)) {
        cats[match(names(xlevels), dimnames(RRPX)[[2]])] <- unlist(lapply(xlevels, 
            length))
    }
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(rpart.control))
        indx <- match(names(extraArgs), controlargs, nomatch = 0)
        if (any(indx == 0)) 
            stop("Argument ", names(extraArgs)[indx == 0], "not matched")
    }
    controls <- rpart.control(...)
    if (!missing(control)) 
        controls[names(control)] <- control
    xval <- controls$xval
    if (is.null(xval) || (length(xval) == 1 && xval == 0) || 
        method == "user") {
        xgroups <- 0
        xval <- 0
    }
    if (missing(cost)) 
        cost <- rep(1, nvar)
    else {
        if (length(cost) != nvar) 
            stop("Cost vector is the wrong length")
        if (any(cost <= 0)) 
            stop("Cost vector must be positive")
    }
    tfun <- function(x) {
        if (is.matrix(x)) 
            rep(is.ordered(x), ncol(x))
        else is.ordered(x)
    }
    isord <- unlist(lapply(m[attr(Terms, "term.labels")], tfun))
    rpfit <- .C(rpart:::C_s_to_rp, n = as.integer(nobs), nvarx = as.integer(nvar),
    ncat = as.integer(cats * (!isord)), method = as.integer(method.int),
    as.double(unlist(controls)), parms = as.double(unlist(init$parms)),
    as.integer(xval), as.integer(xgroups), as.double(t(init$y)),
    as.double(RRPX), as.integer(!is.finite(RRPX)), error = character(1),
    wt = as.double(wt), as.integer(init$numy), as.double(cost),
    NAOK = TRUE)

    #rpfit <- .C("s_to_rp", n = as.integer(nobs), nvarx = as.integer(nvar),
    #   ncat = as.integer(cats * (!isord)), method = as.integer(method.int),
    #   as.double(unlist(controls)), parms = as.double(unlist(init$parms)),
    #   as.integer(xval), as.integer(xgroups), as.double(t(init$y)),
    #   as.double(RRPX), as.integer(!is.finite(RRPX)), error = character(1),
    #   wt = as.double(wt), as.integer(init$numy), as.double(cost),
    #   NAOK = TRUE, PACKAGE="rpart")
    if (rpfit$n == -1) 
        stop(rpfit$error)
    nodes <- rpfit$n
    nsplit <- rpfit$nvarx
    numcp <- rpfit$method
    ncat <- rpfit$ncat[1]
    numresp <- init$numresp
    if (nsplit == 0) 
        xval <- 0
    cpcol <- if (xval > 0 && nsplit > 0) 
        5
    else 3
    if (ncat == 0) 
        catmat <- 0
    else catmat <- matrix(integer(1), ncat, max(cats))
    rp <- .C("s_to_rp2", as.integer(nobs), as.integer(nsplit), 
        as.integer(nodes), as.integer(ncat), as.integer(cats * 
            (!isord)), as.integer(max(cats)), as.integer(xval), 
        which = integer(nobs), cptable = matrix(double(numcp * 
            cpcol), nrow = cpcol), dsplit = matrix(double(1), 
            nsplit, 3), isplit = matrix(integer(1), nsplit, 3), 
        csplit = catmat, dnode = matrix(double(1), nodes, 3 + 
            numresp), inode = matrix(integer(1), nodes, 6),
			PACKAGE="rpart")
	return(rp$which)
}


rrp.dist <- function (X, treated = NULL, msplit = 10, Rep = 250, cut.in = 15, 
    check.bal = FALSE, plot = FALSE, asdist = FALSE, verbose=0) 
{
    n <- dim(X)[1]
    n.var <- dim(X)[2]
    RRP <- newXPtr(n, Rep)
    t.att <- NULL
    Y <- X
    if (cut.in > 0) {
     if(verbose>1)
 	  cat("splitting the support of covariates...\n")
        for (i in 1:n.var) {
            if (is.numeric(X[, i]) | is.integer(X[, i])) {
                if (length(unique(X[, i])) > cut.in) 
                  Y[, i] <- ordered(cut(X[, i], seq(min(X[, i], 
                    na.rm = TRUE), max(X[, i], na.rm = TRUE), 
                    length = cut.in), include.lowest = TRUE))
                else Y[, i] <- ordered(X[, i])
            }
        }
    }
    if(verbose>1)
	 cat("preprocessing data...\n")
    for (i in 1:dim(X)[2]) X[, i] <- as.numeric(X[, i])
    x <- cbind(z = numeric(n), Y)
    rm(Y)
	RRPX <- NULL
    if(verbose>0)
	 cat("RRP running now\n")


	RRP.rpart.setup <- function (formula,
		data, weights, subset, na.action = na.rpart, 
		method, model = FALSE, x = FALSE, y = TRUE, parms, control, 
		cost, ...) 
	{
		call <- match.call()
		if (is.data.frame(model)) {
			m <- model
			model <- FALSE
		} else {
			m <- match.call(expand.dots = FALSE)
			m$model <- m$method <- m$control <- NULL
			m$x <- m$y <- m$parms <- m$... <- NULL
			m$cost <- NULL
			m$na.action <- na.action
			m[[1]] <- as.name("model.frame")
			m <- eval(m, parent.frame())
		}
		RRPX <<- RRP.rpart.matrix(m)
		m
	}	

	mRRP <- RRP.rpart.setup(z ~ ., data = x, method = "anova", minsplit = msplit, 
            xval = 0, cp = 0, maxcompete=0, maxsurrogate=0)


    for (K in 1:Rep) {
        if (plot) {
            if (is.null(treated)) 
                plot(X[, 1:2])
            else plot(X[, 1:2], col = treated + 2, pch = ifelse(treated, 
                20, 17), cex = 1)
        }
        mod <- RRP.rpart(z ~ ., data = x, method = "anova", minsplit = msplit, 
            xval = 0, cp = 0, m=mRRP, RRPY=runif(n), RRPX=RRPX)
		group <- sapply(unique(mod), function(x) which(mod == x))

        n.leaves <- length(group)
        if (!is.null(treated) & check.bal) {
            idx.t <- lapply(group, function(x) x[which(treated[x] == 
                1)])
            idx.c <- lapply(group, function(x) x[which(treated[x] == 
                0)])
        }
        new.g <- vector(2 * length(group), mode = "list")
        if (check.bal) {
         check.for.bal <- function(g) {
            similar <- unlist(group[[g]])
            residual <- numeric(0)
            if (length(idx.t[[g]]) > 0 & length(idx.c[[g]]) > 
                0) {
                mins <- apply(X[idx.t[[g]], ], 2, function(x) min(x, 
                  na.rm = TRUE))
                maxs <- apply(X[idx.t[[g]], ], 2, function(x) max(x, 
                  na.rm = TRUE))
                test <- sapply(idx.c[[g]], function(x) prod(mins <= 
                  X[x, ]) * prod(maxs >= X[x, ]))
                test <- as.logical(test)
                similar <- c(idx.t[[g]], (idx.c[[g]])[test])
                residual <- (idx.c[[g]])[!test]
            }
            new.g[[2 * g - 1]] <<- similar
            new.g[[2 * g]] <<- residual
            if (plot) {
                minX <- min(X[similar, 1], na.rm = TRUE)
                maxX <- max(X[similar, 1], na.rm = TRUE)
                minY <- min(X[similar, 2], na.rm = TRUE)
                maxY <- max(X[similar, 2], na.rm = TRUE)
                rect(minX, minY, maxX, maxY, lty = 3)
                if (length(residual) > 0) {
                  minX <- min(X[residual, 1], na.rm = TRUE)
                  maxX <- max(X[residual, 1], na.rm = TRUE)
                  minY <- min(X[residual, 2], na.rm = TRUE)
                  maxY <- max(X[residual, 2], na.rm = TRUE)
                  rect(minX, minY, maxX, maxY, lty = 3)
                }
             }
            }
            sapply(1:n.leaves, check.for.bal)
            group <- new.g
        }
        addXPtr(RRP, group, -1)
        rm(mod)
		if(verbose>0){
         if (K%/%100 == K/100) {
            cat("0\n")
         } else {
            if (K%/%10 == K/10) 
                cat("+")
            else cat(".")
           }
		}
    }
	rm(RRPX)
	tmp <- mulXPtr(RRP, list(1:n), 1/Rep)
    t.att <- NULL
    if (!is.null(treated)) 
        t.att <- as.logical(treated)
    RRPcl <- class(RRP)
    if (asdist) 
        tmp <- XPtrToDist(RRP) # XPtrToDist apparently removes class info
    attributes(tmp) <- list(Size = n, Diag = FALSE, Upper = FALSE, 
        method = "RRP", minsplit = msplit, replications = Rep, 
        cov.cut = cut.in, balanced = check.bal, treated = t.att, 
        call = match.call())
    if (asdist) 
        class(tmp) <- "dist"
    else class(tmp) <- RRPcl
    return(invisible(tmp))
}
	

# The following code is not exported by rpart but we need it
# Original copyright follows for the rpart:::rpart.matrix
#SCCS  @(#)rpart.matrix.s	1.6 04/02/01
#
# This differs from tree.matrix in xlevels -- we don't keep NULLS in
#   the list for all of the non-categoricals
#
RRP.rpart.matrix <- function(frame)
    {
    if(!inherits(frame, "data.frame"))
	    return(as.matrix(frame))
    frame$"(weights)" <- NULL
    terms <- attr(frame, "terms")
    if(is.null(terms)) predictors <- names(frame)
    else {
	a <- attributes(terms)
	predictors <- as.character(a$variables)[-1] # R change
	removals <- NULL
	if((TT <- a$response) > 0) {
	    removals <- TT
	    frame[[predictors[TT]]] <- NULL
	    }
	if(!is.null(TT <- a$offset)) {
	    removals <- c(removals, TT)
	    frame[[predictors[TT]]] <- NULL
	    }
	if(!is.null(removals)) predictors <- predictors[ - removals]
        labels <- a$term.labels
	if(abs(length(labels)-length(predictors))>0)
	  predictors <- predictors[match(labels,predictors)]
	}

    factors <- sapply(frame, function(x) !is.null(levels(x)))
    characters <- sapply(frame, is.character)
    if(any(factors | characters)) {
	# change characters to factors
	for (preds in predictors[characters])
		frame[[preds]] <- as.factor(frame[[preds]])
        factors <- factors | characters
        column.levels <- lapply(frame[factors], levels)

	# Now make them numeric
	for (preds in predictors[factors])
	     frame[[preds]] <- as.numeric(frame[[preds]])
	x <- as.matrix(frame)
	attr(x, "column.levels") <- column.levels
	}
    else x <- as.matrix(frame[predictors])
    class(x) <- "rpart.matrix"
    x
    }


# The following code is not exported by rpart but we need it
# SCCS  @(#)formatg.s	1.3 06/06/01
# format a set of numbers using C's "g" format
#  It is applied on an element by element basis, which is more
#  appropriate for rpart output than the standard Splus format()
#  command.
# For instance if x=(123, 1.23, .00123)
#	  format(x) = "123.00000", "1.23000", "0.00123"
#  but formatg does not add all of those zeros to the first two numbers
#

RRP.formatg <- function (x, digits = unlist(options("digits")), 
    format = paste("%.", digits, "g", sep = "")) 
{
    if (!is.numeric(x)) 
        stop("x must be a numeric vector")
    n <- length(x)
    temp <- sprintf(format, x)
    if (is.matrix(x)) 
        matrix(temp, nrow = nrow(x))
    else temp
}

# The following code is not exported by rpart but we need it
# Original copyright follows for the rpart:::rpart.anova
#SCCS @(#)rpart.anova.s	1.4 05/02/01

RRP.rpart.anova <- function(y, offset, parms, wt) {
    if (!is.null(offset)) y <- y-offset
    list(y=y, parms=0, numresp=1, numy=1,
	 summary= function(yval, dev, wt, ylevel, digits ) {
	     paste("  mean=", RRP.formatg(yval, digits),
		   ", MSE=" , RRP.formatg(dev/wt, digits),
		   sep='')
	     },
	 text= function(yval, dev, wt, ylevel, digits, n, use.n ) {
	     if(use.n) {paste(RRP.formatg(yval,digits),"\nn=", n,sep="")} else
	               {paste(RRP.formatg(yval,digits))}}

	 )
    }
