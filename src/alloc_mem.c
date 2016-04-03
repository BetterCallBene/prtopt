#include <util.h>
/* wrapper for malloc */
void *wrp_malloc (rqci n, size_t size)
{
    return (malloc (MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *wrp_calloc (rqci n, size_t size)
{
    return (calloc (MAX (n,1), size)) ;
}

/* wrapper for free */
void *wrp_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    p = NULL;
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}

/* wrapper for realloc */
void *wrp_realloc (void *p, rqci n, size_t size, rqci *ok)
{
    void *pnew ;
    pnew = realloc (p, MAX (n,1) * size) ; /* realloc the block */
    *ok = (pnew != NULL) ;                  /* realloc fails if pnew is NULL */
    return ((*ok) ? pnew : p) ;             /* return original p if failure */
}