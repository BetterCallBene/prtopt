#ifndef _CS_TYPES_H
#define _CS_TYPES_H

#include <sys_types.h>

#define CS_VER 3                    /* CSparse Version */
#define CS_SUBVER 1
#define CS_SUBSUB 7
#define CS_DATE "Feb 1, 2016"       /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2016"


/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    csi nzmax ;     /* maximum number of entries */
    csi m ;         /* number of rows */
    csi n ;         /* number of columns */
    csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *i ;        /* row indices, size nzmax */
    realtype *x ;     /* numerical values, size nzmax */
    csi nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    csi *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    csi *q ;        /* fill-reducing column permutation for LU and QR */
    csi *parent ;   /* elimination tree for Cholesky and QR */
    csi *cp ;       /* column pointers for Cholesky, row counts for QR */
    csi *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    csi m2 ;        /* # of rows for QR, after adding fictitious rows */
    realtype lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    realtype unz ;    /* # entries in U for LU; in R for QR */
} css ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs *L ;         /* L for LU and Cholesky, V for QR */
    cs *U ;         /* U for LU, R for QR, not used for Cholesky */
    csi *pinv ;     /* partial pivoting for LU */
    realtype *B ;     /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    csi *p ;        /* size m, row permutation */
    csi *q ;        /* size n, column permutation */
    csi *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    csi *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    csi nb ;        /* # of blocks in fine dmperm decomposition */
    csi rr [5] ;    /* coarse row decomposition */
    csi cc [5] ;    /* coarse column decomposition */
} csd ;


#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))
#endif

