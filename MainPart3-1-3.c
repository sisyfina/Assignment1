#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "mmio.h"
#include <inttypes.h>

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

int mergeArrays(int *arr1, int *arr2, int n1, int n2)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int clash = 0;

  int *arr3 = (int *)malloc((n1+n2) * sizeof(int *));

  while (i<n1 && j<n2)
  {
    if(arr1[i] < arr2[j])
      {
        arr3[k] = arr1[i];
        i++;
      }
    else
      {
        arr3[k] = arr2[j];
        j++;
      }
      k++;
  }

  while(i<n1)
  {
    arr3[k] = arr1[i];
    k++;
    i++;
  }

  while(j<n2)
  {
    arr3[k] = arr2[j];
    k++;
    j++;
  }

  for(int fl=1; fl<n1+n2; fl++)
  {
    if(arr3[fl-1] == arr3[fl])
      clash++;
  }

  return clash;

};

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

	I = (int *) malloc(2*nz * sizeof(int));
	J = (int *) malloc(2*nz * sizeof(int));
	val = (int *) malloc(2*nz * sizeof(int));

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
      I[flag+1] = J[flag];
			J[flag+1] = I[flag];
			val[flag+1]=1;
			flag = flag+2;
			NNZ = NNZ+2;
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

	// V4 PART

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

  // C matrix
  //int n = N;
  //int nnz = NNZ;
  int *c = (int *)malloc(NNZ * sizeof(int));
	//int *row_index = (int *)malloc(nnz * sizeof(int));
	//int *col_index = (int *)malloc((n+1) * sizeof(int));
	//col_index[0] = 0;
	//col_index[n] = nnz;
/*
  for(int i=0; i<nnz; nnz++)
  {
    row_index[i] = ROW_INDEX[i];
    col_index[i] = COL_INDEX[i];
  }
*/
  int elem = 0;
  int c3idx = 0;

  // masking
  for(int i=0; i<N; i++)
  {
    int inumElem = COL_INDEX[i+1] - COL_INDEX[i];
    printf("i%dnumElem = %d\n", i, inumElem);
    if(inumElem > 0)
    {
      // find i's neighbours
      int *cineigh = (int *)malloc(inumElem * sizeof(int *));
      printf("cineigh =");
      for(int flag = 0; flag<inumElem; flag++)
      {
        int index = COL_INDEX[i] + flag;
        cineigh[flag] = ROW_INDEX[index];
        printf(" %d", cineigh[flag]);
      }
      printf("\n");

      // for column i, find j's of  i's elements
      for(int fl=0; fl<inumElem; fl++)
      {
        int index = COL_INDEX[i] + fl;
        int j = ROW_INDEX[index];

        // find j's neighbours
        int jnumElem = COL_INDEX[j+1] - COL_INDEX[j];
        if(jnumElem > 0)
        {
          printf("j%dnumElem = %d\n", j, jnumElem);
          int *cjneigh = (int *)malloc(jnumElem * sizeof(int *));
          printf("cjneigh =");
          for(int flag = 0; flag<jnumElem; flag++)
          {
            int index2 = COL_INDEX[j] + flag;
            cjneigh[flag] = ROW_INDEX[index2];
            printf(" %d", cjneigh[flag]);
          }
          printf("\n");

          c[elem] = mergeArrays(cineigh, cjneigh, inumElem, jnumElem);
          printf("c[%d] = %d\n", elem ,c[elem]);
          c3[i] = c3[i] + c[elem];
          elem++;
        }
      }
    }
    c3[i] = c3[i]/2;
  }

  printf("v =");
	for(int i = 0; i<NNZ; i++) {
		printf(" %d", c[i]);
	}
	printf("\n");

  printf("c3 =");
	for(int i = 0; i<N; i++) {
		printf(" %d", c3[i]);
	}
	printf("\n");



}
