#include <time.h>
  #include <stdio.h>
#include <string.h>

  int main(void)
  {
#ifdef WITH_MPI
  	printf("YEAH\n");
#endif
    return 0;
  }