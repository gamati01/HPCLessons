#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>

//#define N_BLOCK 16
#define nn 512

#include "inc_precision.h"

struct timeval start, stop;
int64_t sec    ;
int64_t usec   ;
double elapsed ;

REAL a[nn][nn];       /** matrixes**/
REAL b[nn][nn];
REAL c[nn][nn];
REAL check[nn][nn];

/*---------------------------------------------------------*/


__global__ void gpu_mm(REAL* d_a, REAL* d_b, REAL* d_c, int n) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    REAL tmp = 0;

    for (int k = 0; k < n; ++k) {
        REAL d_a_tile = d_a[row*n+k];
        REAL d_b_tile = d_b[k*n+col];
        tmp += d_a_tile * d_b_tile;
    }
    d_c[row*nn+col] = tmp;
}

/*---------------------------------------------------------*/

int main()
{
  int i, j, k;
  int N_BLOCK;
  float time1, time2, dub_time;
  float gpu_elapsed_time_ms;

  // Allocate memory space on the device
  REAL *d_a, *d_b, *d_c;
  cudaMalloc((void **) &d_a, sizeof(REAL)*nn*nn);
  cudaMalloc((void **) &d_b, sizeof(REAL)*nn*nn);
  cudaMalloc((void **) &d_c, sizeof(REAL)*nn*nn);

  N_BLOCK=nn/32;

  dim3 dimGrid(N_BLOCK, N_BLOCK);
  dim3 dimBlock(nn/N_BLOCK, nn/N_BLOCK);

  printf("===============================\n");
  if(sizeof(a[1][1]) == sizeof(float))
    printf(" single precision\n");
  else
    printf("double precision\n");
  printf("size %d \n",nn );
  printf("#   block %d \n", N_BLOCK*N_BLOCK);
  printf("dim block %d \n", (nn/N_BLOCK)*(nn/N_BLOCK));
  printf("===============================\n");
  printf("Initialization\n");

  /* initialize matrix */
  time1 = clock();
  for (j = 0; j < nn; j++) {
    for (i = 0; i < nn; i++) {
      a[j][i] = ((REAL)rand())/((REAL)RAND_MAX);
      b[j][i] = ((REAL)rand())/((REAL)RAND_MAX);
      c[j][i] = 0.0L;		
      check[j][i] = 0.0L;		
    }
  }

  time2 = clock();
  dub_time = (time2 - time1)/(double) CLOCKS_PER_SEC;
  printf("Elapsed time for initialization \n");
  printf("Total time -----------------> %f \n", dub_time);

  gettimeofday(&start, 0);
  time1 = clock();

  // some events to count the execution time
  cudaEvent_t custart, custop;
  cudaEventCreate(&custart);
  cudaEventCreate(&custop);
  cudaEventRecord(custart, 0);
  cudaEventSynchronize(custart);

  // copy matrix A and B from host to device memory
  cudaMemcpy(d_a, a, sizeof(REAL)*nn*nn, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, sizeof(REAL)*nn*nn, cudaMemcpyHostToDevice);

  gpu_mm <<< dimGrid, dimBlock >>> (d_a, d_b, d_c, nn);

  // Transfer results from device to host
  cudaMemcpy(c, d_c, sizeof(REAL)*nn*nn, cudaMemcpyDeviceToHost);

  // time counting terminate
  cudaEventRecord(custop, 0);
  cudaEventSynchronize(custop);

  time2 = clock();
  gettimeofday(&stop, 0);

  // compute time elapsed on GPU computing
  cudaEventElapsedTime(&gpu_elapsed_time_ms, custart, custop);

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
  printf("Time   (3)-------------> %f \n", gpu_elapsed_time_ms/1000);
  printf("Mflops (3)-------------> %f \n", 
          2.0*nn*nn*nn/(1000*1000*gpu_elapsed_time_ms/1000));

  /* simple check */
  printf("Check -----------------> %f \n", c[nn/2][nn/2]);

#ifdef VALIDATION
  for (j = 0; j < nn; j++) { 
      for (k = 0; k < nn; k++) { 
          for (i = 0; i < nn; i++) { 
              check[j][i] = check[j][i]+a[j][k]*b[k][i];		
          } 
      } 
  } 
//
  for (j = 0; j < nn; j++) { 
      for (i = 0; i < nn; i++) { 
          printf("Error ---------> %lf, \n", check[j][i]-c[j][i]); 
      } 
  } 
#else
// do nothing
#endif

   return 0;  
}

