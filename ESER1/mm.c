#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define nn (1024)

#include "inc_precision.h"

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
  /*                              */ 
  /* write here 3 nested loop     */ 
  /*                              */ 
                 c[i][j] = c[i][j] + a[i][k]*b[k][j];


  time2 = clock();
  dub_time = (time2 - time1)/(double) CLOCKS_PER_SEC;
  printf("===============================\n");
  printf("Tme -----------------> %f \n", dub_time);
  printf("Mflops ----------------> %f \n", 
          2.0*nn*nn*nn/(1000*1000*dub_time));

  /* simple check */
  printf("Check -------------> %f \n", c[nn/2-1][nn/2-1]);

   return 0;  
}
