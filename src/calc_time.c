#include <sys_types.h>
#include <time.h>
realtype tic (void) { return (clock () / (realtype) CLOCKS_PER_SEC) ; }
realtype toc (realtype t) { realtype s = tic () ; return (MAX (0, s-t)) ; }