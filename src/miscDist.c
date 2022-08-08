/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2002	S. M. Iacus
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Exports
 *	addDist(...)
 *  setDist(...)
 *  newXPtr(...)
 *
 * to be called as  .Call(.)  in ../R/miscDist.R
 */



#include <R.h>
#include <Rmath.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Complex.h>


/* dist object manipulation */      
SEXP addDist(SEXP d, SEXP x, SEXP k);
SEXP setDist(SEXP d, SEXP x, SEXP k);

SEXP newXPtr(SEXP n, SEXP k);
SEXP setXPtr(SEXP XPtr, SEXP x, SEXP k);
SEXP addXPtr(SEXP XPtr, SEXP x, SEXP k);
SEXP mulXPtr(SEXP XPtr, SEXP x, SEXP k);
SEXP XPtrToNumeric(SEXP XPtr);

SEXP applyXPtr(SEXP XPtr, SEXP x, SEXP sub, SEXP f, SEXP rho);


static SEXP RRP_dist_tag;

#define CHECK_RRP_OBJECT(s) do { \
    if (TYPEOF(s) != EXTPTRSXP || \
        R_ExternalPtrTag(s) !=  RRP_dist_tag) \
        error("bad RRP dist object"); \
} while (0)

/* addDist/setDist:
   ------
   
   parameters:
   -----------
   
   d   : the pointer to the vector containing elements of the dist object
   x   : a list of indexes to change in the original matrix
   k   : a vector of constants to add vector-wise
*/



SEXP addDist(SEXP d, SEXP x, SEXP k)
{
	int i, j, h, pos, len, N;
	double *dObj, *kappa, *myd;
	int *xIdx, n2, n1, i1, i2;
	
	if(!isNumeric(d)) error("`d' must be numeric");
	if(!isNumeric(k)) error("`k' must be numeric");
	
	PROTECT(d = AS_NUMERIC(d));
	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	
	n2 = LENGTH(x);
	if (n2 == NA_INTEGER)
     error("error");
	
	dObj = NUMERIC_POINTER(d);

	kappa = NUMERIC_POINTER(k);
	len =  LENGTH(d);		  
    N = (int)(0.5*(1+sqrt(1+8*len)));
	myd = REAL(d);

	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   i1 = xIdx[i];
	   i2 = xIdx[j];
	   if(i1 > i2){ 
	    i2 = xIdx[i];
	    i1 = xIdx[j];
	   }
	   pos = (i1-1)*(2*N-i1)/2 + i2-i1;
	   myd[pos-1] = dObj[pos-1] + kappa[h];
	  }
	 }
	}
	UNPROTECT(3);
    return(d);
}

SEXP setDist(SEXP d, SEXP x, SEXP k)
{
	int h, i, j, pos, len, N;
	double  *kappa,  *myd;
	int  *xIdx, n2, n1, i1, i2;
	
	if(!isNumeric(d)) error("`d' must be numeric");
	if(!isNumeric(k)) error("`k' must be numeric");
	
	PROTECT(d = AS_NUMERIC(d));
	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	
	n2 = LENGTH(x);
	if (n2 == NA_INTEGER)
		error("error");
	
	len =  LENGTH(d);		  
    N = (int)(0.5*(1+sqrt(1+8*len)));
	kappa = NUMERIC_POINTER(k);

    myd = REAL(d);			  

	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   i1 = xIdx[i];
	   i2 = xIdx[j];
	   if(i1 > i2){ 
	    i2 = xIdx[i];
	    i1 = xIdx[j];
	   }
	   pos = (i1-1)*(2*N-i1)/2 + i2-i1;
	   myd[pos-1] = kappa[h];
	  }
	 }
	}

	UNPROTECT(3);
	return(d);
}

/* 
   DocInfo for the following functions: 
    newXPtr, addXPtr, mulXPtr, setXPtr, XPtrToNumeric
   R interfaces have the names but XPtrToNumeric which is called 
   by XPtrToDist	
	
   Preamble: XPtr are external pointer to dist objects. Not really objects of
   class 'dist' in the R sense, but something you should think about as.
   XPtr should be tought as the lower triangular part of a symmetrix matrix, 
   say M. The diagonal is not included and the user must assume it to be 
   identically zero. The equivalent for a true dist object 'd' obtained
   from a symmetric matrix 'M' would be  

   > M <- matrix(0, 10, 10)
   > d <- as.dist(M) 
   > str(d)
   Class 'dist'  atomic [1:45] 0 0 0 0 0 0 0 0 0 0 ...
     ..- attr(*, "Size")= int 10
     ..- attr(*, "call")= language as.dist.default(m = M)
     ..- attr(*, "Diag")= logi FALSE
     ..- attr(*, "Upper")= logi FALSE
	
   All the functions access the XPtr as if it was a matrix (see example below)

   parameters:
   -----------
   n    : the length of the numeric vector to be created as external pointer
   XPtr : an already existing XPtr, i.e. external pointer with tag RRP_dist_tag
   x   : a list of vectors of indexes to change in the original matrix
   k   : a numeric vector of constants to add vector-wise
   the length of k and x must match. It is up to the R code to fix this.
   
   value: an XPtr

   newXPtr(n, k)       : creates a new XPtr object of Size n, initializes it 
                         with constant k
   addXPtr(XPtr, x, k) : adds a vector fo constants
   mulXPtr(XPtr, x, k) : multiply by a vector fo constants
   setXPtr(XPtr, x, k) : puts contants in a XPtr object				   
   XPtrToNumeric(XPtr) : return a SEXP containg a numeric vector. This not 
                         copied
example:
   Consider the above matrix M. This can be initialized by the newXPtr as 
   follows
   a <- newXPtr(10, 0)

   if we now set 
   > x <- list(c(1,7,4), 2:3) 
   > k <- c(1,-1)
   
   and we call 
   
   > addXPtr(a, x, k)

   this is equivalent to
   > M[ c(1,7,4), c(1,7,4) ] <- M[ c(1,7,4), c(1,7,4) ] + 1
   > M[ 2:3, 2:3 ] <- M[ 2:3, 2:3 ] + (-1)
   
   then "as.dist(M)" and "(XPtrToDist(a))" will return the same output (of 
   course the diagonal of M has been changed and this could not be recovered
   in a dist object without diagonal. It will be up to the user to fix things)
*/

SEXP newXPtr(SEXP n, SEXP k){
	int i, N, *np;
	SEXP distC, distR, dim;
    double K,*kp, *d;
		
	if(!isNumeric(k)) error("`k' must be numeric");
	if(!isInteger(n)) error("`n' must be an integer");
	
	np = INTEGER_POINTER(n);
	kp = NUMERIC_POINTER(k);
    N = *np;
	K = *kp;

	PROTECT(distC = allocVector(REALSXP, N*(N-1)/2)); 
    d = REAL(distC);
	for(i=0; i<N*(N-1)/2; i++)
	 d[i] = K;

	distR = R_MakeExternalPtr(d, RRP_dist_tag,  distC);
    UNPROTECT(1);
	return(distR);
}


SEXP setXPtr(SEXP XPtr, SEXP x, SEXP k)
{
	int h, i, j, pos;
	double  *kappa, *d;
	int  len, N;
	int  *xIdx, n2, n1, i1, i2;
		
	CHECK_RRP_OBJECT(XPtr); 
    d = R_ExternalPtrAddr(XPtr);

	if(!isNumeric(k)) error("`k' must be numeric");

	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	
	n2 = LENGTH(x); 
    if (n2 == NA_INTEGER)
     error("`x' is zero-length");
	 
	kappa = NUMERIC_POINTER(k);
	len =  LENGTH(R_ExternalPtrProtected(XPtr));		  
    N = (int)(0.5*(1+sqrt(1+8*len)));

	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   i1 = xIdx[i];
	   i2 = xIdx[j];
	   if(i1 > i2){ 
	    i2 = xIdx[i];
	    i1 = xIdx[j];
	   }
	   pos = (i1-1)*(2*N-i1)/2 + i2-i1;
	   d[pos-1] = kappa[h];
	  }
	 }
	}

	UNPROTECT(2);
	return(XPtr);
}



SEXP addXPtr(SEXP XPtr, SEXP x, SEXP k)
{
	int i, j, h, pos;
	double *kappa, *d;
	int  *xIdx, n2, n1;
	int  len, N,i1,i2;

	CHECK_RRP_OBJECT(XPtr); 

    d = R_ExternalPtrAddr(XPtr);
	if(!isNumeric(k)) error("`k' must be numeric");
	
	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	
	n2 = LENGTH(x); 
    if (n2 == NA_INTEGER)
     error("`x' is zero-length");
	 	
	kappa = NUMERIC_POINTER(k);
	len =  LENGTH(R_ExternalPtrProtected(XPtr));		  
    N = (int)(0.5*(1+sqrt(1+8*len)));

	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   i1 = xIdx[i];
	   i2 = xIdx[j];
	   if(i1 > i2){ 
	    i2 = xIdx[i];
	    i1 = xIdx[j];
	   }
	   pos = (i1-1)*(2*N-i1)/2 + i2-i1;
	   d[pos-1] = d[pos-1] + kappa[h];
	  }
	 }
	}
	UNPROTECT(2);
    return(XPtr);
}

SEXP mulXPtr(SEXP XPtr, SEXP x, SEXP k)
{
	int i, j, h, pos;
	double *kappa, *d;
	int  *xIdx, n2, n1;
	int  len, N, i1, i2;

	CHECK_RRP_OBJECT(XPtr); 

    d = R_ExternalPtrAddr(XPtr);
	if(!isNumeric(k)) error("`k' must be numeric");
	
	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	
	n2 = LENGTH(x); 
    if (n2 == NA_INTEGER)
     error("`x' is zero-length");
	 	
	kappa = NUMERIC_POINTER(k);
	len =  LENGTH(R_ExternalPtrProtected(XPtr));		  
    N = (int)(0.5*(1+sqrt(1+8*len)));

	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   i1 = xIdx[i];
	   i2 = xIdx[j];
	   if(i1 > i2){ 
	    i2 = xIdx[i];
	    i1 = xIdx[j];
	   }
	   pos = (i1-1)*(2*N-i1)/2 + i2-i1;
	   d[pos-1] = d[pos-1] * kappa[h];
	  }
	 }
	}
	UNPROTECT(2);
    return(XPtr);
}


/* 
   DocInfo for the following functions: 
   applyXPtr
   the R counterpart has the same name 
   applies function 'f' to the elements M[i, sub] where M is the symmetric
   matrix imagined to be associated with XPtr with 0's on the diagonal and
   'i' varies in 'idx' (a vector of indexes).	   
*/

SEXP applyXPtr(SEXP XPtr, SEXP idx, SEXP sub, SEXP f, SEXP rho)
{
	int i, j, h, pos;
	double *d, tmpval;
	int  IDX, SUB, nidx, nsub, n2, n1;
	int  len, N, np=0;
    SEXP  val, tmp, subval;
	int *mysub, *myidx;
	double *mysubval;
	SEXP R_fcall; 
	
	CHECK_RRP_OBJECT(XPtr); 

	PROTECT(R_fcall = allocList(2));
	np++;
	SETCAR(R_fcall, f);
	SET_TYPEOF(R_fcall, LANGSXP);

    d = R_ExternalPtrAddr(XPtr);
	if(!isInteger(idx)) error("`idx' must be integer");
	if(!isInteger(sub)) error("`sub' must be integer");
	
	PROTECT(idx = AS_INTEGER(idx));
	myidx = INTEGER(idx);
	np++;
	PROTECT(sub = AS_INTEGER(sub));
	mysub = INTEGER(sub);
	np++;
	
	nidx = LENGTH(idx); 
    if (nidx == NA_INTEGER)
     error("`idx' is zero-length");

	nsub = LENGTH(sub); 
    if (nsub == NA_INTEGER)
     error("`sub' is zero-length");
	 	
	len =  LENGTH(R_ExternalPtrProtected(XPtr));		  
    N = (int)(0.5*(1+sqrt(1+8*len)));

    PROTECT(val = allocVector(VECSXP, nidx));
    np++;
	for(i=0; i<nidx; i++){	 
	 IDX = myidx[i];
	 PROTECT(subval = NEW_NUMERIC(nsub));
	 mysubval = REAL(subval);
	 np++;
	 for(j = 0; j<nsub; j++){
	  tmpval = 0.0;
	  SUB = mysub[j];
	  if(IDX<SUB){
		pos = (IDX-1)*(2*N-IDX)/2 + SUB-IDX;
	    tmpval = d[pos-1];
	   }
	   if(IDX>SUB){
		pos = (SUB-1)*(2*N-SUB)/2 + IDX-SUB;
	    tmpval = d[pos-1];
	   }
	  mysubval[j] = tmpval;
     }	  
 	 SETCADR(R_fcall, subval);
	 SET_VECTOR_ELT(val, i, eval(R_fcall, rho) );
	}

	UNPROTECT(np);
    return(val);
}
/* This function does not copy the object, it just returns
   the SEXP cached inside external pointer XPtr
*/
SEXP XPtrToNumeric(SEXP XPtr)
{
	SEXP val;
	CHECK_RRP_OBJECT(XPtr); 
	PROTECT(val =  R_ExternalPtrProtected(XPtr));
	UNPROTECT(1);
	return(val);
}


static R_CMethodDef R_CDef[] = {
   {"addDist", (DL_FUNC)&addDist, 3},
   {"setDist", (DL_FUNC)&setDist, 3},
   {"newXPtr", (DL_FUNC)&newXPtr, 2},   
   {"addXPtr", (DL_FUNC)&addXPtr, 3},
   {"mulXPtr", (DL_FUNC)&mulXPtr, 3},
   {"setXPtr", (DL_FUNC)&setXPtr, 3},
   {"applyXPtr", (DL_FUNC)&applyXPtr, 5},
   {"XPtrToNumeric", (DL_FUNC)&XPtrToNumeric, 1},
   {NULL, NULL, 0},
};


void
R_init_rrp(DllInfo *info)
{
    R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
    RRP_dist_tag = install("RRP_DIST_TAG");
}


