#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// v is number of nodes
// A is the adjacent matrix

int main(int argc, char *argv[] ) {

	if( argc == 2 ) {
      		printf("The argument supplied is %s\n", argv[1]);
   	}
   	else if( argc > 2 ) {
      		printf("Too many arguments supplied.\n");
   	}
   	else {
      		printf("One argument expected.\n");
   	}

	struct timespec ts_start, ts_end;
	long delta_sec;

	// int v
	int v = atoi(argv[1]);
	printf("v = %d\n", v);
	// boolean **A;
	int **A = (int **)malloc(v * sizeof(int *));
  for (int i=0; i<v; i++)
      A[i] = (int *)malloc(v * sizeof(int));

	// random graph filling
	srand((unsigned int)time(NULL));
	for(int i = 0; i<v; i++){
		for (int j=i; j<v; j++){
      A[i][j] = rand() % 2;
			A[j][i] = A[i][j];
      printf("A[%d][%d] = %d ", i, j ,A[i][j]);
    }
		printf("\n");
    for(int k = -1; k<i; k++)
			printf("            ");
	}
	printf("\n");

	// vector with number of triangle that correspond to a node
	int *c3 = (int *)malloc(v * sizeof(int *));
	// initialize c3
	for(int i = 0; i<v; i++) {
		c3[i] = 0;
	}

	clock_gettime(CLOCK_MONOTONIC, &ts_start); // get the start time

	for(int i = 0; i<v; i++){
	for(int j = 0; j<v; j++){
	for(int k = 0; k<v; k++){
		if ((i != j)&&(j != k)&&(k != i)){
		if ((A[i][j] == 1)&&(A[j][k] == 1)&&(A[k][i] == 1)){
			c3[i]++;
			c3[j]++;
			c3[k]++;
		}}
	}}}

	clock_gettime(CLOCK_MONOTONIC, &ts_end); // get the finishing time
  delta_sec = ts_end.tv_nsec - ts_start.tv_nsec;
  printf("Execution time: %ld ns\n", delta_sec);

	printf("c3 =");
	for(int i = 0; i<v; i++) {
		printf(" %d", c3[i]);
	}
	printf("\n");

	clock_gettime(CLOCK_MONOTONIC, &ts_start); // get the start time

	// initialize c3
	for(int i = 0; i<v; i++) {
		c3[i] = 0;
	}

	for(int i = 0; i<v-2; i++){
	for(int j = i+1; j<v-1; j++){
	for(int k = j+1; k<v; k++){
		//if ((i != j)&&(j != k)&&(k != i)){
		if ((A[i][j] == 1)&&(A[j][k] == 1)&&(A[k][i] == 1)){
			c3[i]++;
			c3[j]++;
			c3[k]++;
		}//}
	}}}

	clock_gettime(CLOCK_MONOTONIC, &ts_end); // get the finishing time
  delta_sec = ts_end.tv_nsec - ts_start.tv_nsec;
  printf("Execution time: %ld ns\n", delta_sec);

	printf("c3 =");
	for(int i = 0; i<v; i++) {
		printf(" %d", c3[i]);
	}
	printf("\n");


}
