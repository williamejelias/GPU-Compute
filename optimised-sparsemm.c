#include "utils.h"

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
                        const COO, const COO, const COO,
                        COO *);

struct _dataEntry {
  int i, j;
  double data;
};


int cmpfunc(const void *p, const void *q)
{
    int l = ((struct _dataEntry *)p)->j;
    int r = ((struct _dataEntry *)q)->j;
    if ((l - r) == 0) {
      int li = ((struct _dataEntry *)p)->i;
      int ri = ((struct _dataEntry *)q)->i;
      return(li - ri);
    }
    return (l - r);
}


void checkA(const COO A, const int* Ams) {
  int outI;
  int dimX = A->m;

  for (outI = 0; outI < dimX; outI ++) {
    int Astart = Ams[outI];
    int Aend = Ams[outI + 1]-1;
    int rowIndex = Astart;
    // printf("Row %i has values: ", outI);
    while (rowIndex <= Aend) {
      int aJ = A -> coords[rowIndex].j;
      int aI = A-> coords[rowIndex].i;
      double aData = A->data[rowIndex];
      printf("%f ", aData);
      rowIndex ++;
    }
    // printf("\n");
  }
}

void checkB(const COO B, const int* Bns) {
  int outJ;
  int dimY = B->n;


  for (outJ = 0; outJ < dimY; outJ ++) {
    int Bstart = Bns[outJ];
    int Bend = Bns[outJ + 1]-1;
    int colIndex = Bstart;
    // printf("col %i has values: ", outJ);
    while (colIndex <= Bend) {
      int bJ = B -> coords[colIndex].j;
      int bI = B -> coords[colIndex].i;
      double bData = B->data[colIndex];
      printf("%f ", bData);
      colIndex ++;
    }
    // printf("\n");
  }
}



/* Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
  //printf("basic\n");
  // return basic_sparsemm(A, B, C);

  int Bm = B->m;
  int Bn = B->n;
  int BNZ = B->NZ;

  // printf("A NNZ: %i, m:%i, n:%i\n", A->NZ, A->m, A->n);
  // printf("B NNZ: %i, m:%i, n:%i\n", BNZ, Bm, Bn);

  // create array of structs for B
  struct _dataEntry *B_entries = NULL;
  B_entries = malloc( BNZ * sizeof(struct _dataEntry));
  if (B_entries == NULL) {
    fprintf(stderr, "Out of memory");
    exit(1);
  }

  // printf("Allocated B...\n");

  int b_entry_index;
  #pragma vector aligned
  #pragma acc parallel loop
  for (b_entry_index = 0; b_entry_index < BNZ; b_entry_index++) {
    B_entries[b_entry_index].i = B->coords[b_entry_index].i;
    B_entries[b_entry_index].j = B->coords[b_entry_index].j;
    B_entries[b_entry_index].data = B->data[b_entry_index];
  }
  // printf("Sorting B...\n");
  qsort(B_entries, BNZ, sizeof(struct _dataEntry), cmpfunc);

  int outputNNZ = 0;

  //--------------------------------------------------------------------------
  // Convert A CSR
  //--------------------------------------------------------------------------
  // printf("A to CSR format...\n");

  int *Ams = NULL;
  Ams = calloc((A->m + 1), sizeof(Ams));
  Ams[0] = 0;

  int *Bns = NULL;
  Bns = calloc((B->n + 1), sizeof(Bns));
  Bns[0] = 0;

  if (Ams == NULL || Bns == NULL) {
    fprintf(stderr, "Out of memory");
    exit(1);
  }

  int ANZ = A -> NZ;
  int AMIndex = 0;
  int rowsize = 0;

  // record current row within A
  int currentRow = A->coords[0].i;
  if (currentRow > 0) {
    int pad;
    // #pragma acc parallel loop
    #pragma vector aligned
    for (pad = 0; pad < currentRow; pad ++) {
      Ams[AMIndex + 1] = 0;
      AMIndex ++;
    }
  }

  // index in A of start of row
  int currentRowIndexA = 0;
  // index in A of end of row
  int endOfRowIndexA = 0;

  // iterate through all rows within A
  int Aindex;
  // #pragma acc parallel loop
  for (Aindex = 0; Aindex < ANZ; Aindex ++) {
    if (A -> coords[Aindex].i != currentRow) {
      Ams[AMIndex + 1] = Ams[AMIndex] + rowsize;
      AMIndex ++;

      // pad zeroes for completely empty rows
      int diff = A -> coords[Aindex].i - currentRow;
      if (diff > 1) {
        int pad;
	//#pragma acc parallel loop
        #pragma vector aligned
        for (pad = 0; pad < diff - 1; pad ++) {
          Ams[AMIndex + 1] = Ams[AMIndex];
          AMIndex ++;
        }
      }

      // move to next row
      currentRow = A->coords[Aindex].i;
      currentRowIndexA = Aindex;
      endOfRowIndexA = currentRowIndexA;

      rowsize = 1;

    } else if (A -> coords[Aindex].i == currentRow) {
      // increment index towards end of row
      endOfRowIndexA = Aindex;
      rowsize ++;
    }
    if (Aindex == ANZ - 1) {
      endOfRowIndexA = Aindex;
      Ams[AMIndex + 1] = ANZ;
      AMIndex ++;
    }
  }

  // pad end
  int pad;
  #pragma vector aligned
  // #pragma acc parallel loop
  for (pad = 0; pad < A->m-currentRow-1; pad ++) {
    Ams[AMIndex + 1] = Ams[AMIndex];
    AMIndex ++;
  }

  // printf("\n");
  // checkA(A, Ams);


  //--------------------------------------------------------------------------
  // Convert B CSC
  //--------------------------------------------------------------------------
  // printf("B to CSC format...\n");
  int BNIndex = 0;
  int colsize = 0;

  // record current column
  int currentColumn = B_entries[0].j;
  if (currentColumn > 0) {
    int pad;
    // #pragma acc parallel loop
    #pragma vector aligned
    for (pad = 0; pad < currentColumn; pad ++) {
      Bns[BNIndex + 1] = 0;
      BNIndex ++;
    }
  }
  // index in B of start of column
  int currentColumnIndexB = 0;
  // index in B of end of column
  int endOfColumnIndexB = 0;

  int Bindex;
  // #pragma acc parallel loop
  for (Bindex = 0; Bindex < BNZ; Bindex ++) {
    if (B_entries[Bindex].j != currentColumn) {
      Bns[BNIndex + 1] = Bns[BNIndex] + colsize;
      BNIndex ++;

      // pad zeroes for completely empty rows
      int diff = B_entries[Bindex].j - currentColumn;
      if (diff > 1) {
        int pad;
	//#pragma acc parallel loop
        #pragma vector aligned
        for (pad = 0; pad < diff - 1; pad ++) {
          Bns[BNIndex + 1] = Bns[BNIndex];
          BNIndex ++;
        }
      }

      // move to next column
      currentColumn = B_entries[Bindex].j;
      currentColumnIndexB = Bindex;
      endOfColumnIndexB = currentColumnIndexB;

      colsize = 1;

    } else if (B_entries[Bindex].j == currentColumn) {
      // increment index towards end of column
      endOfColumnIndexB = Bindex;
      colsize ++;
    }
    // found end of column
    if (Bindex == BNZ - 1) {
      endOfColumnIndexB = Bindex;
      Bns[BNIndex + 1] = Bns[BNIndex] + colsize;
      BNIndex ++;
    }
  }

  int pad2;
  // #pragma acc parallel loop
  #pragma vector aligned
  for (pad2 = 0; pad2 < Bn-currentColumn; pad2 ++) {
    Bns[BNIndex + 1] = Bns[BNIndex];
    BNIndex ++;
  }

  // printf("\n");
  // checkB(B, Bns);

  // printf("Allocating sparse...\n");
  // int outputNonZeroes = getNonZeroes(A, B, *C, Ams, Bns);
  int outputNonZeroes = getNonZeroes(A, Bn, B_entries, Ams, Bns);
  alloc_sparse(A -> m, Bn, outputNonZeroes, C);

  //--------------------------------------------------------------------------
  // MULTIPLY CSR CSC v1
  // --------------------------------------------------------------------------
  // printf("multiplying...\n");
  int outputI;
  int outputJ;
  int dimX = A->m;
  int dimY = Bn;

  // row size is the number of items in each row of A
  // this is also equal to the number of columns in B
  int rowSize = A->n;

  #pragma acc parallel loop
  for (outputI = 0; outputI < dimX; outputI ++) {
    #pragma acc loop
    for (outputJ = 0; outputJ < dimY; outputJ ++) {

      int Astart = Ams[outputI];
      int Aend = Ams[outputI + 1]-1;

      int Bstart = Bns[outputJ];
      int Bend = Bns[outputJ + 1]-1;

      int aJ;
      int bI;

      double result = 0.0;
      int rowIndex = Astart;
      int colIndex = Bstart;
      while (rowIndex <= Aend && colIndex <= Bend) {
        aJ = A -> coords[rowIndex].j;
        bI = B_entries[colIndex].i;
        if (aJ == bI) {
          result += A -> data[rowIndex] * B_entries[colIndex].data;
          rowIndex ++;
          colIndex ++;
        } else if (aJ < bI) {
          rowIndex ++;
        } else if (aJ > bI) {
          colIndex ++;
        }
      }

      if (result != 0.0) {
        // printf("NNZ %i/%i value of %f at coordinates (%d, %d)\n", outputNNZ, outputNonZeroes, result, outputI, outputJ);
        (*C) -> coords[outputNNZ].i = outputI;
        (*C) -> coords[outputNNZ].j = outputJ;
        (*C) -> data[outputNNZ] = result;
        outputNNZ ++;
      }
      // printf("output cell %i %i has value %f\n", outputI, outputJ, result);
    }
  }

  //--------------------------------------------------------------------------
  // May be better for PARALLEL MULTIPLY
  //--------------------------------------------------------------------------
  // printf("Multiplying...\n");
  // int outputI;
  // int outputJ;
  // int dimX = A->m;
  // int dimY = B->n;
  //
  // // row size is the number of items in each row of A
  // // this is also equal to the number of columns in B
  // int rowSize = A->n;
  //
  // //#pragma acc kernels
  // for (outputI = 0; outputI < dimX; outputI ++) {
  //   double *rowVals = NULL;
  //   rowVals = calloc(rowSize, sizeof(*rowVals));
  //   if (rowVals == NULL) {
  //     fprintf(stderr, "Out of memory");
  //     exit(1);
  //   }
  //
  //   int Astart = Ams[outputI];
  //   int Aend = Ams[outputI + 1]-1;
  //   // build dense row
  //   int rowIndex;
  //
  //   //#pragma acc parallel loop
  //   for (rowIndex = Astart; rowIndex <= Aend; rowIndex++ ) {
  //       rowVals[A -> coords[rowIndex].j] = A -> data[rowIndex];
  //   }
  //
  //   //#pragma acc parallel loop
  //   for (outputJ = 0; outputJ < dimY; outputJ ++) {
  //     double *colVals = NULL;
  //     colVals = calloc(rowSize, sizeof(*colVals));
  //     if (colVals == NULL) {
  //       fprintf(stderr, "Out of memory");
  //       exit(1);
  //     }
  //
  //     int Bstart = Bns[outputJ];
  //     int Bend = Bns[outputJ + 1]-1;
  //     // build dense col
  //     // All in each col
  //     int colIndex;
  //
  //     //#pragma acc parallel loop
  //     for (colIndex = Bstart; colIndex <= Bend; colIndex++ ) {
  //       colVals[B_entries[colIndex].i] = B_entries[colIndex].data;
  //     }
  //
  //     int rowcol;
  //     double result = 0.0;
  //
  //     #pragma acc parallel loop reduction(+:result)
  //     for (rowcol = 0; rowcol < rowSize; rowcol ++) {
  //       result += rowVals[rowcol] * colVals[rowcol];
  //     }
  //
  //     if (result != 0.0) {
  //       // printf("NNZ %i/%i value of %f at coordinates (%d, %d)\n", outputNNZ, outputNonZeroes, result, outputI, outputJ);
  //       (*C) -> coords[outputNNZ].i = outputI;
  //       (*C) -> coords[outputNNZ].j = outputJ;
  //       (*C) -> data[outputNNZ] = result;
  //       outputNNZ ++;
  //     }
  //     free(colVals);
  //   }
  //   free(rowVals);
  // }
  //


  free(Ams);
  free(Bns);
  free(B_entries);
  // printf("Multiplied.\n");
  // printf("output matrix has dimensions %i x %i\n", (*C)->m, (*C)->n);
  // printf("Output Matrix has dimensions %d x %d, and has %d Non-Zero Entries.\n\n", A -> m, B -> n, outputNNZ);
}



int getNonZeroes(const COO A, int Bn, struct _dataEntry* arr, const int* Ams, const int* Bns) {
  // printf("Calculating Non Zeroes...\n");
  int outNZ = 0;
  int outputI;
  int outputJ;
  int dimX = A->m;
  int dimY = Bn;

  // row size is the number of items in each row of A
  // this is also equal to the number of columns in B
  int rowSize = A->n;

  #pragma acc parallel loop
  for (outputI = 0; outputI < dimX; outputI ++) {
    #pragma acc loop
    for (outputJ = 0; outputJ < dimY; outputJ ++) {
      int Astart = Ams[outputI];
      int Aend = Ams[outputI + 1]-1;

      int Bstart = Bns[outputJ];
      int Bend = Bns[outputJ + 1]-1;

      int aJ;
      int bI;

      int addedVal = 0;
      int rowIndex = Astart;
      int colIndex = Bstart;
      while (rowIndex <= Aend && colIndex <= Bend) {
        aJ = A -> coords[rowIndex].j;
        bI = arr[colIndex].i;
        if (aJ == bI) {
          rowIndex ++;
          colIndex ++;
          addedVal = 1;
        } else if (aJ < bI) {
          rowIndex ++;
        } else if (aJ > bI) {
          colIndex ++;
        }
      }

      if (addedVal == 1) {
        outNZ ++;
      }
    }
  }

  return outNZ;
}


// Slower - replaced by addThreeMatrices2
void addThreeMatrices(COO m1, COO m2, COO m3, COO out)
{
  int m1index = 0;
  int m2index = 0;
  int m3index = 0;
  int outNZ = 0;

  // i and j markers for each matrix
  int m1I, m1J, m2I, m2J, m3I, m3J;

  // min i and j markers for entry in output
  int minI, minJ;

  // addition
  while (m1index < m1 -> NZ || m2index < m2 -> NZ || m3index < m3 -> NZ) {
    m1I = m1 -> coords[m1index].i;
    m1J = m1 -> coords[m1index].j;

    m2I = m2 -> coords[m2index].i;
    m2J = m2 -> coords[m2index].j;

    m3I = m3 -> coords[m3index].i;
    m3J = m3 -> coords[m3index].j;

    minI = m1 -> m;
    minJ = m1 -> n;

    float result = 0.0;

    if (m1index != m1 -> NZ) if (m1I < minI) minI = m1I;
    if (m2index != m2 -> NZ) if (m2I < minI) minI = m2I;
    if (m3index != m3 -> NZ) if (m3I < minI) minI = m3I;

    if ((m1I == minI) && (m1J < minJ)) minJ = m1J;
    if ((m2I == minI) && (m2J < minJ)) minJ = m2J;
    if ((m3I == minI) && (m3J < minJ)) minJ = m3J;

    out -> coords[outNZ].i = minI;
    out -> coords[outNZ].j = minJ;

    if (m1I == minI && m1J == minJ) {
      result += m1->data[m1index];
      m1index ++;
    }
    if (m2I == minI && m2J == minJ) {
      result += m2->data[m2index];
      m2index ++;
    }
    if (m3I == minI && m3J == minJ) {
      result += m3->data[m3index];
      m3index ++;
    }

    out -> data[outNZ] = result;
    outNZ ++;
  }
}


// sum all elements in 3 array
//#pragma acc routine worker
double sum_array(double* arr) {
  double result = 0.0;
  result = arr[0] + arr[1] + arr[2];
  return result;
}


// #pragma acc routine worker
void addThreeMatrices2(COO m1, COO m2, COO m3, COO* out) {
  int m1index = 0;
  int m2index = 0;
  int m3index = 0;
  int m1NZ = m1 -> NZ;
  int m2NZ = m2 -> NZ;
  int m3NZ = m3 -> NZ;
  int outNZ = 0;

  // i and j markers for each matrix
  int m1I, m1J, m2I, m2J, m3I, m3J;

  // min i and j markers for entry in output
  int minI, minJ;

  // size of each row and column
  int rowSize = m1 -> m;
  int colSize = m1 -> n;

  // Passover to find indices of rows and columns to add
  // object that has 3 arrays are of size outNZ
  // array of output i's
  // array of output j's
  // array vals which contains 3array of vals
  // 3 array of vals contains the m1, m2, m3 nz or zeroes
  // sum of this 3array becomes out -> data[nz]
  // initially give these arrays size of biggest NZ of input matrix
  int maxNZ = 0;
  maxNZ += m1NZ;
  maxNZ += m2NZ;
  maxNZ += m3NZ;

  int *outIs = NULL;
  outIs = malloc(maxNZ*sizeof(*outIs));

  int *outJs = NULL;
  outJs = malloc(maxNZ*sizeof(*outJs));

  double** additions = malloc(maxNZ * sizeof(double*));
  if(outIs == NULL || outJs == NULL || additions == NULL){
    fprintf(stderr, "Out of memory");
    exit(1);
  }
  int i;
  for(i = 0; i < maxNZ; i++){
    double* p = calloc(3, sizeof(double));
    if(p == NULL){
      fprintf(stderr, "Out of memory");
      exit(1);
    }
    additions[i] = p;
  }

  // addition
  // printf("m1 NZ: %d, m2 NZ: %i, m3 NZ: %i\n", m1 -> NZ, m2 -> NZ, m3 -> NZ);
  while (m1index < m1NZ || m2index < m2NZ || m3index < m3NZ) {
    minI = rowSize;
    minJ = colSize;
    m1I = minI;
    m1J = minJ;
    m2I = minI;
    m2J = minJ;
    m3I = minI;
    m3J = minJ;

    if (m1index < m1NZ){
      m1I = m1 -> coords[m1index].i;
      m1J = m1 -> coords[m1index].j;
      if (m1I <= minI) minI = m1I;
    }

    if (m2index < m2NZ) {
      m2I = m2 -> coords[m2index].i;
      m2J = m2 -> coords[m2index].j;
      if (m2I <= minI) minI = m2I;
    }

    if (m3index < m3NZ) {
      m3I = m3 -> coords[m3index].i;
      m3J = m3 -> coords[m3index].j;
      if (m3I <= minI) minI = m3I;
    }

    if ((m1I == minI) && (m1J <= minJ)) minJ = m1J;
    if ((m2I == minI) && (m2J <= minJ)) minJ = m2J;
    if ((m3I == minI) && (m3J <= minJ)) minJ = m3J;

    outIs[outNZ] = minI;
    outJs[outNZ] = minJ;

    if (m1index < m1NZ) {
      if (m1I == minI && m1J == minJ) {
        additions[outNZ][0] = m1->data[m1index];
        m1index ++;
      }
    }
    if (m2index < m2NZ) {
      if (m2I == minI && m2J == minJ) {
        additions[outNZ][1] = m2->data[m2index];
        m2index ++;
      }
    }
    if (m3index < m3NZ) {
      if (m3I == minI && m3J == minJ) {
        additions[outNZ][2] = m3->data[m3index];
        m3index ++;
      }
    }
    outNZ ++;
  }

  alloc_sparse(m1 -> m, m1 -> n, outNZ, out);

  int index;
  // CHECK THIS PRAGMA FOR EFFICIENCY
  // #pragma vector always
  #pragma acc parallel loop
  #pragma vector aligned
  for (index = 0; index < outNZ; index ++) {
    (*out) -> coords[index].i = outIs[index];
    (*out) -> coords[index].j = outJs[index];
    (*out) -> data[index] = sum_array(additions[index]);
  }

  int addition_index;
  #pragma acc parallel loop
  for(addition_index = 0; addition_index < maxNZ; addition_index++) {
    free(additions[addition_index]);
  }
  free(additions);
  free(outIs);
  free(outJs);
}


// ./sparsemm result.matrix ../small-matrices/DG2-ip-laplace-2D.matrix ../small-matrices/DG2-ip-laplace-2D.matrix  ../small-matrices/DG2-ip-laplace-2D.matrix  ../small-matrices/DG2-mass-2D.matrix  ../small-matrices/DG2-mass-2D.matrix  ../small-matrices/DG2-mass-2D.matrix


/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
                            const COO D, const COO E, const COO F,
                            COO *O) {
  // return basic_sparsemm_sum(A, B, C, D, E, F, O);

  // generate skeleton sparse matrices for G&H which will hold the sum of (A+B+C) and (D+E+F) respectively
  COO G;
  COO H;

  // printf("m1 NZ: %d, m2 NZ: %i, m3 NZ: %i\n", A -> NZ, B -> NZ, C -> NZ);

  // addThreeMatrices(A, B, C, G);
  // addThreeMatrices(D, E, F, H);
  addThreeMatrices2(A, B, C, &G);
  addThreeMatrices2(D, E, F, &H);

  // call optimised_sparsemm on G & H to output O
  optimised_sparsemm(G, H, O);
  free_sparse(&G);
  free_sparse(&H);
}
