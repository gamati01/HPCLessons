#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h> 

#define nn (4096)

#include "inc_precision.h"

//time
struct timeval start, stop;
int64_t sec    ;
int64_t usec   ;
double elapsed ;
  
REAL a[nn][nn];       /** matrixes**/
REAL b[nn][nn];
REAL c[nn][nn];

int main()
{
  int k, i, j;
  float time1, time2, dub_time;

  printf("===============================\n");
  if(sizeof(a[1][1]) == sizeof(float))
    printf(" single precision\n");
  else
    printf("double precision\n");
  printf("size %d \n",nn );
  printf("===============================\n");
  printf("Initialization\n");

  /* initialize matrix */
  time1 = clock();
  for (j = 0; j < nn; j++)
    {
    for (i = 0; i < nn; i++)
      {
      a[j][i] = ((REAL)rand())/((REAL)RAND_MAX);
      b[j][i] = ((REAL)rand())/((REAL)RAND_MAX);
      c[j][i] = 0.0L;		
      }
    }
  time2 = clock();
  dub_time = (time2 - time1)/(double) CLOCKS_PER_SEC;
  printf("Elapsed time for initialization \n");
  printf("Total time -----------------> %f \n", dub_time);

  time1 = clock();
  gettimeofday(&start, 0);
  #pragma acc kernels
  for (i = 0; i < nn; i++) {
       for (k = 0; k < nn; k++) {
            for (j = 0; j < nn; j++) {
                 c[i][j] = c[i][j] + a[i][k]*b[k][j];
	    }
        }
   }

  time2 = clock();
  gettimeofday(&stop, 0);

  dub_time = (time2 - time1)/(double) CLOCKS_PER_SEC;

  // RESULTS
  sec    = stop.tv_sec  - start.tv_sec;
  usec   = stop.tv_usec - start.tv_usec;
  elapsed = sec + usec*1E-6;
  //
  printf("===============================\n");
  printf("Time   (1)-------------> %f \n", dub_time);
  printf("Mflops (1)-------------> %f \n", 
          2.0*nn*nn*nn/(1000*1000*dub_time));
  printf("Time   (2)-------------> %f \n", elapsed);
  printf("Mflops (2)-------------> %f \n", 
          2.0*nn*nn*nn/(1000*1000*elapsed));

  /* simple check */
  printf("Check -------------> %f \n", c[nn/2-1][nn/2-1]);

   return 0;  
}
