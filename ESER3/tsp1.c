/* tsp1.c -- Experiments on TSP algorithms
    Input lines from standard input:
        <algnum> <n>
    Run algnum on first n points in input file
    FNAME (lines in x,y format)
 */

/* Codice originale di J. Bentley DDJ 1999      */
/* modifiche (piccole) di G. Amati, CASPUR 2005 */

#include <stdio.h>
#include <time.h>
#include <math.h>

/* Main Defines and Globals*/

#define COUNTDISTS 0
#define MAXN 30
#define FNAME "lincoln.txt"
#define INF 1e38
typedef double Dist;
typedef double Dist2;

int n;

typedef struct point {
    float x;
    float y;
} Point;
Point c[MAXN];

int p[MAXN];
Dist minsum;
Dist distarr[MAXN][MAXN];
int minp[MAXN];

/* Edge Operations */

int distcnt;

void getgeom()
{	FILE *fp;
	fp = fopen(FNAME, "r");
	n = 0;
	while (fscanf(fp, "%f %f", &c[n].x, &c[n].y) != EOF) {
		n++;
	}
}

float sqr(float x)
{   return x*x;
}

Dist d(int i, int j)
{   
	return (Dist) (sqrt(sqr(c[i].x-c[j].x)+sqr(c[i].y-c[j].y)));
}

Dist2 d2(int i, int j)
{
        return (Dist2) (sqr(c[i].x-c[j].x)+sqr(c[i].y-c[j].y));
}

void initdistarr()
{
        int i, j;
        for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                        distarr[i][j] = d(i, j);
}


/* Support */

void swap(int i, int j)
{	int t = p[i];
	p[i] = p[j];
	p[j] = t;
}

/* ------------------------------------------------------------*/
/* VERSION 1 -- Base: n*n! */

void save(Dist sum)
{	int i;
	if (sum < minsum) {
		minsum = sum;
		for (i = 0; i < n; i++) {
			minp[i] = p[i];
		}
                printf(" new minimum found = %g \n", minsum);
	}
}

void check1()
{	int i;
	Dist sum = d(p[0], p[n-1]);
	for (i = 1; i < n; i++)
		sum += d(p[i-1], p[i]);
	save(sum);
}

void search1(int m)
{	int i;
	if (m == 1) {
		check1();
	} else {
		for (i = m-1; i >= 0; i--) {
			swap(i, m-1);
			search1(m-1);
			swap(i, m-1);
		}
	}
}

void solve1()
{	search1(n);
}

/* ------------------------------------------------------------*/
/* VERSION 2 -- Fix City n-1: (n-1)!*n */

void solve2()
{	search1(n-1);
}

/* ------------------------------------------------------------*/
/* VERSION 3 -- Keep sum of distances so far: (1+e)*(n-1)! */

void check3(Dist sum)
{	sum += d(p[0], p[n-1]);
	save(sum);
}

void search3(int m, Dist sum)
{	int i;
	if (m == 1) {
		check3(sum + d(p[0], p[1])); /* bug if n==1? */
	} else {
		for (i = m-1; i >= 0; i--) {
			swap(i, m-1);
			search3(m-1, sum + d(p[m-1], p[m]));
			swap(i, m-1);
		}
	}
}

void solve3()
{	search3(n-1, 0);
}

/* ------------------------------------------------------------*/
/* VERSION 4 -- Prune search */

void search4(int m, Dist sum)
{	int i;
	if (sum > minsum) return;
	if (m == 1) {
		check3(sum + d(p[0], p[1]));
	} else {
		for (i = m-1; i >= 0; i--) {
			swap(i, m-1);
			search4(m-1, sum + d(p[m-1], p[m]));
			swap(i, m-1);
		}
	}
}

void solve4()
{	search4(n-1, 0);
}

/* ------------------------------------------------------------*/
/* VERSION 5 -- Prune search + precompute sqrt */

void check5(Dist sum)
{	sum += distarr[p[0]][ p[n-1]];
	save(sum);
}

void search5(int m, Dist sum)
{	int i;
	if (sum > minsum) return;
	if (m == 1) {
		check5(sum + distarr[p[0]][p[1]]);
	} else {
		for (i = m-1; i >= 0; i--) {
			swap(i, m-1);
			search5(m-1, sum + distarr[p[m-1]][p[m]]);
			swap(i, m-1);
		}
	}
}

void solve5()
{	search5(n-1, 0);
}

/* ------------------------------------------------------------*/
/* DRIVERS */

int main()
{	int i, alg, start;
	float secs;
        FILE *fout1,*fout2;
	getgeom();
	while (scanf("%d %d", &alg, &n) != EOF) {
		distcnt = 0;
                fout1 = fopen("original.dat","w");
                fout2 = fopen("minimum.dat","w");
		printf("-------------------------------------\n");
		printf(" Travel Salesman problem             \n");
		printf(" algorithm   = %d                    \n", alg);
		printf(" cities      = %d                    \n", n);
		printf("-------------------------------------\n");
		printf("starting point                       \n");
		for (i = 0; i < n; i++) {
			p[i] = i;
               		printf(" %d %g %g \n", i, c[i].x,  c[i].y);  
               		fprintf(fout1," %d %f %f \n", i, c[i].x,  c[i].y);  
		}
		printf("-------------------------------------\n");
		minsum = INF;
		start = clock();
		switch (alg) {
		case 1: solve1(); break;
		case 2: solve2(); break;
		case 3: solve3(); break;
		case 4: solve4(); break;
		case 5: {
			initdistarr(); 
			solve5(); 
			break;
			}
		}
		secs = ((float) clock() - (float) start)
			/ (float) CLOCKS_PER_SEC;
		printf("------------------------------------\n");
		printf(" time (sec) = %f                    \n", secs);
		printf(" distance   = %f                    \n", (float) minsum);
		printf("------------------------------------\n");
		printf(" minimal path                       \n");
                for (i = 0; i < n; i++) {
                     printf(" %d %lf %lf \n", minp[i], c[minp[i]].x,  c[minp[i]].y);  
               	     fprintf(fout2," %d %f %f \n", minp[i], c[minp[i]].x,  c[minp[i]].y);  
		}
		printf("------------------------------------\n");
                fclose(fout1);
                fclose(fout2);
		break;
	}
}


