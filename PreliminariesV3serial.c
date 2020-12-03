#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "mmio.h"
#include <inttypes.h>

#define SERIAL_ELEMENT_LIMIT 11

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

int binarySearch(int *vector, int start, int end, int value)
{
  int numElem = end - start + 1;
  int i = start;

  if(numElem < 1)
    return 1;
  else if(numElem<SERIAL_ELEMENT_LIMIT)
  {
    for(int tmp=0; tmp<numElem; tmp++)
    {
      if(value == vector[start+tmp])
        return 0;
    }
    return 1;
  }
  else if(numElem%2 == 0)
    i = i + numElem/2;
  else
    i = i + (numElem+1)/2;


  if(value == vector[i])
    return 0;
  else if(value < vector[i])
  {
    binarySearch(vector, start, i-1, value);
  }
  else
  {
    binarySearch(vector, i+1, end, value);
  }
}

// n is number of nodes
// A is the adjacent matrix

int main(int argc, char *argv[] ) {

	// COO variables

	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz, NNZ;
	int *I, *J;
	int *val;

/*
*	2 arguments for now
*	2nd argument: is the MatrixMarket FILE
*
*/

	if( argc == 2 ) {
		if ((f = fopen(argv[1], "r")) == NULL)
			exit(1);
  }
  else if( argc > 2 )
	{
  		printf("\nToo many arguments supplied.\n");
			exit(1);
  }
  else
	{
  		printf("\nArgument supply error.\n");
			exit(1);
  }

	// get matcode

	if (mm_read_banner(f, &matcode) != 0)
	{
			printf("Could not process Matrix Market banner.\n");
			exit(1);
	}

	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */

	if (!(mm_is_matrix(matcode) && mm_is_coordinate(matcode) &&
					mm_is_symmetric(matcode) ) || mm_is_complex(matcode))
	{
			printf("Sorry, this application does not support ");
			printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
			exit(1);
	}

	// screen more TO DO (done?)

	/* find out size of sparse matrix .... */

	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
			exit(1);

	// COO
	// triagonal to full
	NNZ = 0;

	/* reseve memory for matrices */

	I = (int *) malloc(nz * sizeof(int));
	J = (int *) malloc(nz * sizeof(int));
	val = (int *) malloc(nz * sizeof(int));

	/* Replace missing or double val column with 1s and change the fscanf to match pattern matrices*/

	int flag = 0;
	for (int i=0; i<nz; i++)
	{
		fscanf(f, "%d %d\n", &I[flag], &J[flag]);
		if(I[flag] != J[flag])
		{
      val[flag]=1;
			I[flag]--;  /* adjust from 1-based to 0-based */
			J[flag]--;
			flag++;
			NNZ++;
		}
	}

	if (f !=stdin) fclose(f);

	/************************/
	/* now write out matrix */
	/************************/

	mm_write_banner(stdout, matcode);
	mm_write_mtx_crd_size(stdout, M, N, NNZ);
	/*
	for (int i=0; i<NNZ; i++)
		fprintf(stdout, "%d %d %d\n", I[i]+1, J[i]+1, val[i]);
	*/

	// CSR format equal to CSC

	int *V = (int *)malloc(NNZ * sizeof(int));
	int *ROW_INDEX = (int *)malloc(NNZ * sizeof(int));
	int *COL_INDEX = (int *)malloc((N+1) * sizeof(int));
	COL_INDEX[0] = 0;
	COL_INDEX[N] = NNZ;

	coo2csc(ROW_INDEX, COL_INDEX, I, J, NNZ, N, 0);

	printf("V =");
	for(int i = 0; i<NNZ; i++)
  {
		V[i] = 1;
		printf(" %d", V[i]);
	}
	printf("\n");

	printf("ROW_INDEX =");
	for(int i = 0; i<NNZ; i++)
  {
		printf(" %d", ROW_INDEX[i]);
	}
	printf("\n");

	printf("COL_INDEX =");
	for(int i = 0; i<N+1; i++)
  {
		printf(" %d", COL_INDEX[i]);
	}
	printf("\n");

	// vector with number of triangles that corresponds to a node
	int *c3 = (int *)malloc(N * sizeof(int *));

	// initialize c3
	for(int i=0; i<N; i++)
		c3[i] = 0;

  //int *edgeNum = (int *)malloc(N * sizeof(int *));

  for(int i=0; i<N-1; i++)
  {
    printf("i = %d\n", i);
    int ci_next = COL_INDEX[i+1];
    int ci = COL_INDEX[i];
    int iEdgeNum = ci_next - ci;
    if(iEdgeNum>0)
      int *cineigh = (int *)malloc(iEdgeNum * sizeof(int *));
    // find i's neighbours
    printf("cineigh =");
    for(int flag = 0; flag<iEdgeNum; flag++)
    {
      cineigh[flag] = ROW_INDEX[ci+flag];
      printf(" %d", cineigh[flag]);
    }
    printf("\n");

    // search by edge
    for(int j=0; j<iEdgeNum; j++)
    {
      int jneigh = cineigh[j];
      printf("j = %d\n", jneigh);
      int cj_next = COL_INDEX[jneigh+1];
      int cj = COL_INDEX[jneigh];
      int jEdgeNum = cj_next - cj;
      if(jEdgeNum>0)
        int *cjneigh = (int *)malloc(jEdgeNum * sizeof(int *));
      // find j's neighbours
      printf("cjneigh =");
      for(int flag = 0; flag<jEdgeNum; flag++)
      {
        cjneigh[flag] = ROW_INDEX[cj+flag];
        printf(" %d", cjneigh[flag]);
      }
      printf("\n");

      for(int k=0; k<jEdgeNum; k++)
      {
        int kneigh = cjneigh[k];
        printf("k = %d\n", kneigh);

        int tmp = binarySearch(cineigh, 0, iEdgeNum-1, kneigh);
        printf("tmp = %d\n", tmp);
        if(tmp == 0)
        {
          c3[i]++;
          c3[cineigh[j]]++;
          c3[cjneigh[k]]++;
        }


      }
    }


  }

  printf("c3 =");
	for(int i = 0; i<N; i++) {
		printf(" %d", c3[i]);
	}
	printf("\n");


}
