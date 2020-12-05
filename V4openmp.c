#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "mmio.h"
#include <inttypes.h>
#include <omp.h>
#include <ctype.h>


//#define OMP_NESTED TRUE

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

int mergeArrays(uint32_t *arr1, uint32_t *arr2, uint32_t n1, uint32_t n2)
{
  uint32_t i = 0;
  uint32_t j = 0;
  uint32_t k = 0;
  uint32_t clash = 0;
  uint32_t fl;

  uint32_t *arr3 = (uint32_t *)malloc((n1+n2) * sizeof(uint32_t *));

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

  for(fl=1; fl<n1+n2; fl++)
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
	uint32_t M, N, nz, NNZ;
	uint32_t *I, *J;
	uint32_t *val;
  uint32_t triangles = 0;

/*
*	2 arguments for now
*	2nd argument: is the MatrixMarket FILE
*
*/

	if( argc == 3 ) {
		if ((f = fopen(argv[1], "r")) == NULL)
			exit(1);
  }
  else if( argc > 3 )
	{
  		printf("\nToo many arguments supplied.\n");
			exit(1);
  }
  else
	{
  		printf("\nArgument supply error.\n");
			exit(1);
  }

  int nthreads = (2<<atoi(argv[2]));

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

	I = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
	J = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
	val = (int *) malloc(2*nz * sizeof(int));

	/* Replace missing or double val column with 1s and change the fscanf to match pattern matrices*/

	uint32_t flag = 0;
	for (uint32_t i=0; i<nz; i++)
	{
		fscanf(f, "%u %u\n", &I[flag], &J[flag]);
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

	uint32_t *V = (int *)malloc(NNZ * sizeof(int));
	uint32_t *ROW_INDEX = (uint32_t *)malloc(NNZ * sizeof(uint32_t));
	uint32_t *COL_INDEX = (uint32_t *)malloc((N+1) * sizeof(uint32_t));
	COL_INDEX[0] = 0;
	COL_INDEX[N] = NNZ;

	coo2csc(ROW_INDEX, COL_INDEX, I, J, NNZ, N, 0);
  /*
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
  */

  // omp variables
  uint32_t i, j;
  uint32_t inumElem, jnumElem;
  uint32_t *cineigh, *cjneigh;
  uint32_t index, index2;
  uint32_t fl;

	// vector with number of triangles that corresponds to a node
	uint32_t *c3;
  c3 = (uint32_t *)malloc(N * sizeof(uint32_t *));

	// initialize c3
  #pragma omp parallel for shared(N,c3) private(i) num_threads(nthreads)
	for(i=0; i<N; i++)
		c3[i] = 0;

  // C matrix
  //uint32_t *c = (uint32_t *)malloc(NNZ * sizeof(uint32_t *));

  uint32_t elem = 0;
  uint32_t cint = 0;


  // allow nested parallelism
  //omp_set_nested(1);
  uint32_t chunk;
  chunk = N/nthreads;


  // masking
  #pragma omp parallel for \
  shared(N, COL_INDEX, ROW_INDEX, elem, c3) \
  private(i, inumElem, cineigh, flag, index, j, fl, index2, cint) \
  schedule(static, chunk) ordered
    for(i=0; i<N; i++)
    {
      uint32_t c3int = 0;
      uint32_t a, b, c, d, e, f, g;
      #pragma omp atomic read
      a = COL_INDEX[i+1];
      #pragma omp atomic read
      b = COL_INDEX[i];

      inumElem = a - b;
      printf("i%dnumElem = %d\n", i, inumElem);
      if(inumElem > 0)
      {
        // find i's neighbours
        cineigh = (uint32_t *)malloc(inumElem * sizeof(uint32_t *));

        //printf("ci%dneigh =", i);
        for(flag = 0; flag<inumElem; flag++)
        {
          index = b + flag;
          #pragma omp atomic read
          c = ROW_INDEX[index];
          #pragma omp atomic write
          cineigh[flag] = c;
          //printf(" %d", c);
        }
        //#pragma omp atomic read
        //&(*d) = &(*cineigh);
        //printf("\n");

        // for column i, find j's of  i's elements
    //    #pragma omp parallel for \
    //    shared(inumElem, COL_INDEX, ROW_INDEX) \
    //    private(fl, index, j, jnumElem, cjneigh, flag, index2)
        for(fl=0; fl<inumElem; fl++)
        {
          index = b + fl;
          #pragma omp atomic read
          j = ROW_INDEX[index];

          // find j's neighbours
          #pragma omp atomic read
          d = COL_INDEX[j+1];
          #pragma omp atomic read
          e = COL_INDEX[j];

          jnumElem = d - e;
          if(jnumElem > 0)
          {
            //printf("j%dnumElem = %d\n", j, jnumElem);
            cjneigh = (uint32_t *)malloc(jnumElem * sizeof(uint32_t *));

            printf("cj%dneigh =", j);

            for(flag = 0; flag<jnumElem; flag++)
            {
              int index2 = e + flag;
              #pragma omp atomic read
              f = ROW_INDEX[index2];
              #pragma omp atomic write
              cjneigh[flag] = f;
              printf(" %d", f);
            }
            printf("\n");

            #pragma omp ordered
            {
              #pragma omp critical
              {
                cint = mergeArrays(cineigh, cjneigh, inumElem, jnumElem);
                //c[elem] = cint;
                //printf("c[%u] = %u\n", elem ,c[elem]);
                #pragma omp atomic read
                g = c3int;
                g = g + cint;
                #pragma omp atomic write
                c3int = g;
                //printf("c3[%u] + c[%u] = %u\n", i, elem, c3int);
                //elem++;
                //printf("elem++ = %u\n", elem);
              }

            }

          }
        }
      }
      #pragma omp atomic write
        c3[i] = c3int/2;

      printf("Hi from thread %d! c3[%d] = %d\n", i, i, c3[i]);


    }


  // prohibid nested parallelism
  //omp_set_nested(0);

  /*printf("v =");
	for(int i = 0; i<NNZ; i++) {
		printf(" %d", c[i]);
	}
	printf("\n");
  */

  #pragma omp parallel for shared(N,c3,triangles) private(i) schedule(static, chunk) ordered
	for(i=0; i<N; i++)
  {
    triangles = triangles + c3[i];
    printf("Hi from thread %d! c3[%d] = %d\n", i, i, c3[i]);
  }


  triangles = triangles/3;

  printf("%u\n", triangles);

  printf("c3 =");
	for(i = 0; i<N; i++) {
		printf(" %u", c3[i]);
	}
	printf("\n");



}
