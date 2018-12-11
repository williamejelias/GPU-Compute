#include "utils.h"

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
                        const COO, const COO, const COO,
                        COO *);


void sortColumnOrder(COO matrix)
{
    // organise B into column major order
    // bubble sort is O(N^2) so shouldnt affect runtime
    int i, j;
    for (i = 0; i < matrix -> NZ - 1; i ++) {
        for (j = 0; j < matrix -> NZ - i - 1; j ++) {
            if (matrix -> coords[j].j > matrix -> coords[j + 1].j) {
                // temp
                int itemp = matrix -> coords[j].i;
                int jtemp = matrix -> coords[j].j;
                double dtemp = matrix -> data[j];

                // swap
                matrix -> coords[j].i = matrix -> coords[j + 1].i;
                matrix -> coords[j].j = matrix -> coords[j + 1].j;
                matrix -> data[j] = matrix -> data[j + 1];

                // replace val
                matrix -> coords[j + 1].i = itemp;
                matrix -> coords[j + 1].j = jtemp;
                matrix -> data[j + 1] = dtemp;
            }
        }
    }
}

/* Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
    // printf("multiplying\n");

    // return basic_sparsemm(A, B, C);
    sortColumnOrder(B);

    int outputNNZ = 0;
    int outputNonZeroes = getNonZeroes(A, B, *C);
    alloc_sparse(A -> m, B -> n, outputNonZeroes, C);

    // printf("A NNZ: %i, m:%i, n:%i\n", A->NZ, A->m, A->n);
    // printf("B NNZ: %i, m:%i, n:%i\n", B->NZ, B->m, B->n);
    // printf("C NNZ: %i, m:%i, n:%i\n", outputNonZeroes, A->m, B->n);

    //--------------------------------------------------------------------------
    // Convert A CSR
    //--------------------------------------------------------------------------
    int *Ams = NULL;
    Ams = malloc((A->m + 1)*sizeof(*Ams));
    Ams[0] = 0;

    int *Bns = NULL;
    Bns = malloc((B->n + 1)*sizeof(*Bns));
    Bns[0] = 0;

    if (Ams == NULL || Bns == NULL) {
      fprintf(stderr, "Out of memory");
      exit(1);
    }


    int AMIndex = 0;
    int rowsize = 0;

    // record current row within A
    int currentRow = A->coords[0].i;
    if (currentRow > 0) {
      int pad;
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
    for (Aindex = 0; Aindex < A->NZ; Aindex ++) {
      if ((A -> coords[Aindex].i != currentRow) || (Aindex == A -> NZ - 1)) {
        // pad zeroes for completely empty rows
        int diff = A -> coords[Aindex].i - currentRow;
        if (diff > 1) {
          int pad;
          for (pad = 0; pad < diff - 1; pad ++) {
            Ams[AMIndex + 1] = 0;
            AMIndex ++;
          }
        }
        // found end of row
        if (Aindex == A -> NZ - 1) {
          endOfRowIndexA = Aindex;
          rowsize ++;
        }

        // move to next row
        currentRow = A->coords[Aindex].i;

        Ams[AMIndex + 1] = Ams[AMIndex] + rowsize;
        AMIndex ++;

        currentRowIndexA = Aindex;
        endOfRowIndexA = currentRowIndexA;

        rowsize = 1;

      } else if (A -> coords[Aindex].i == currentRow) {
        // increment index towards end of row
        endOfRowIndexA = Aindex;
        rowsize ++;
      }
    }


    //--------------------------------------------------------------------------
    // Convert B CSC
    //--------------------------------------------------------------------------
    int BNIndex = 0;
    int colsize = 0;

    // record current column
    int currentColumn = B -> coords[0].j;
    if (currentColumn > 0) {
      int pad;
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
    for (Bindex = 0; Bindex < B -> NZ; Bindex ++) {
      if (B -> coords[Bindex].j != currentColumn || Bindex == B -> NZ - 1) {
        // pad zeroes for completely empty rows
        int diff = B -> coords[Bindex].j - currentColumn;
        if (diff > 1) {
          int pad;
          for (pad = 0; pad < diff - 1; pad ++) {
            Bns[BNIndex + 1] = 0;
            BNIndex ++;
          }
        }
        // found end of column
        if (Bindex == B -> NZ - 1) {
          endOfColumnIndexB = Bindex;
          colsize ++;
        }

        // move to next column
        currentColumn = B -> coords[Bindex].j;

        Bns[BNIndex + 1] = Bns[BNIndex] + colsize;
        BNIndex ++;

        currentColumnIndexB = Bindex;
        endOfColumnIndexB = currentColumnIndexB;

        colsize = 1;
      } else if (B -> coords[Bindex].j == currentColumn) {
        // increment index towards end of column
        endOfColumnIndexB = Bindex;
        colsize ++;
      }
    }
    // printf("total B NNZ:%i out of %i", sum2, B->NZ);




    //--------------------------------------------------------------------------
    // MULTIPLY CSR CSC v1
    //--------------------------------------------------------------------------
    // int AM_counter = 0;
    // // iterate over rows;
    // int Arow;
    // int Acol;
    // double Aval;
    //
    // int BN_counter = 0;
    // // iterate over rows;
    // int Brow;
    // int Bcol;
    // double Bval;
    //
    // // All Columns
    // for (BN_counter = 0; BN_counter < B ->n; BN_counter ++) {
    //     Bcol = BN_counter;
    //     int Bstart = Bns[BN_counter];
    //     int Bend = Bns[BN_counter + 1]-1;
    //
    //     // All in each col
    //     int colIndex;
    //     for (colIndex = Bstart; colIndex <= Bend; colIndex++ ) {
    //         Brow = B -> coords[colIndex].i;
    //         Bval = B -> data[colIndex];
    //
    //         // All Rows
    //         for (AM_counter = 0; AM_counter < A ->m; AM_counter ++) {
    //             Arow = AM_counter;
    //             int Astart = Ams[AM_counter];
    //             int Aend = Ams[AM_counter + 1]-1;
    //
    //             // All in each row
    //             int rowIndex;
    //             double result = 0.0;
    //             for (rowIndex = Astart; rowIndex <= Aend; rowIndex++ ) {
    //                 Acol = A -> coords[rowIndex].j;
    //                 Aval = A -> data[rowIndex];
    //
    //                 if (Acol == Brow) {
    //                     result += Aval * Bval;
    //                 }
    //             }
    //
    //             if (result != 0.0) {
    //                 (*C) -> coords[outputNNZ].i = Arow;
    //                 (*C) -> coords[outputNNZ].j = Bcol;
    //                 (*C) -> data[outputNNZ] = result;
    //                 outputNNZ ++;
    //                 printf("NNZ %i/%ivalue of %f at coordinates (%d, %d)\n", outputNNZ, outputNonZeroes, result, Arow, Bcol);
    //             }
    //         }
    //     }
    // }





    //--------------------------------------------------------------------------
    // May be better for PARALLEL MULTIPLY
    //--------------------------------------------------------------------------
    int outputI;
    int outputJ;
    int dimX = A->m;
    int dimY = B->n;

    // row size is the number of items in each row of A
    // this is also equal to the number of columns in B
    int rowSize = A->n;

    //#pragma acc kernels
    for (outputI = 0; outputI < dimX; outputI ++) {
      double *rowVals = NULL;
      rowVals = calloc(rowSize, sizeof(*rowVals));
      if (rowVals == NULL) {
        fprintf(stderr, "Out of memory");
        exit(1);
      }

      int Astart = Ams[outputI];
      int Aend = Ams[outputI + 1]-1;
      // build dense row
      int rowIndex;

      //#pragma acc parallel loop
      for (rowIndex = Astart; rowIndex <= Aend; rowIndex++ ) {
          rowVals[A -> coords[rowIndex].j] = A -> data[rowIndex];
      }

      //#pragma acc parallel loop
      for (outputJ = 0; outputJ < dimY; outputJ ++) {
        double *colVals = NULL;
        colVals = calloc(rowSize, sizeof(*colVals));
        if (colVals == NULL) {
          fprintf(stderr, "Out of memory");
          exit(1);
        }

        int Bstart = Bns[outputJ];
        int Bend = Bns[outputJ + 1]-1;
        // build dense col
        // All in each col
        int colIndex;

        //#pragma acc parallel loop
        for (colIndex = Bstart; colIndex <= Bend; colIndex++ ) {
          colVals[B -> coords[colIndex].i] = B -> data[colIndex];
        }

        int rowcol;
        double result = 0.0;

        #pragma acc parallel loop reduction(+:result)
        for (rowcol = 0; rowcol < rowsize; rowcol ++) {
          result += rowVals[rowcol] * colVals[rowcol];
        }

        if (result != 0.0) {
          // printf("NNZ %i/%i value of %f at coordinates (%d, %d)\n", outputNNZ, outputNonZeroes, result, outputI, outputJ);
          (*C) -> coords[outputNNZ].i = outputI;
          (*C) -> coords[outputNNZ].j = outputJ;
          (*C) -> data[outputNNZ] = result;
          outputNNZ ++;
        }
        free(colVals);
      }
      free(rowVals);
    }

    free(Ams);
    free(Bns);












    //--------------------------------------------------------------------------
    // Original MULTIPLY
    //--------------------------------------------------------------------------
    // int outputNNZ = 0;
    // int outputNonZeroes = getNonZeroes(A, B, *C);
    // alloc_sparse(A -> m, B -> n, outputNonZeroes, C);
    //
    // // printf("%i, %i", A->m, B->n);
    // printf("A NNZ: %i, m:%i, n:%i\n", A->NZ, A->m, A->n);
    //
    // // record current row within A
    // int currentRow = A->coords[0].i;
    //
    // // index in A of start of row
    // int currentRowIndexA = 0;
    //
    // // index in A of end of row
    // int endOfRowIndexA = 0;
    //
    // // iterate through all values within A
    // int Aindex;
    // for (Aindex = 0; Aindex < A->NZ; Aindex ++) {
    //     // printf("%i\n",A -> coords[Aindex].i );
    //     if ((A -> coords[Aindex].i != currentRow) || (Aindex == A -> NZ - 1)) {
    //         // found end of row
    //         if (Aindex == A -> NZ - 1) {
    //             endOfRowIndexA = Aindex;
    //         }
    //
    //         // iterate over all columns of matrix B
    //         // record current column
    //         int currentColumn = B -> coords[0].j;
    //
    //         // index in B of start of column
    //         int currentColumnIndexB = 0;
    //
    //         // index in B of end of column
    //         int endOfColumnIndexB = 0;
    //
    //         int Bindex;
    //         for (Bindex = 0; Bindex < B -> NZ; Bindex ++) {
    //             if (B -> coords[Bindex].j != currentColumn || Bindex == B -> NZ - 1) {
    //                 // found end of column
    //                 if (Bindex == B -> NZ - 1) {
    //                     endOfColumnIndexB = Bindex;
    //                 }
    //
    //                 double result = 0.0;
    //                 // iterate between start of row and end of row for column in B
    //                 int inRowIndexA;
    //                 for (inRowIndexA = currentRowIndexA; inRowIndexA <= endOfRowIndexA; inRowIndexA++) {
    //                     // for each non zero value within the column, multiply with respective row elements and set respective row element within C
    //
    //                     double aData = A -> data[inRowIndexA];
    //                     int aIValue = A->coords[inRowIndexA].i;
    //                     int aJValue = A->coords[inRowIndexA].j;
    //
    //                     double bData;
    //                     int bIValue;
    //                     int bJValue;
    //
    //                     int inColumnIndexB;
    //                     for (inColumnIndexB = currentColumnIndexB; inColumnIndexB <= endOfColumnIndexB; inColumnIndexB++) {
    //                         bData = B -> data[inColumnIndexB];
    //                         bIValue = B->coords[inColumnIndexB].i;
    //                         bJValue = B->coords[inColumnIndexB].j;
    //
    //                         if (aJValue == bIValue) {
    //                             result += aData * bData;
    //                         }
    //                     }
    //                 }
    //
    //                 if (result != 0.0) {
    //                     (*C) -> coords[outputNNZ].i = currentRow;
    //                     (*C) -> coords[outputNNZ].j = currentColumn;
    //                     (*C) -> data[outputNNZ] = result;
    //                     outputNNZ ++;
    //                     // printf("New output NNZ value of %f at coordinates (%d, %d)\n", result, currentRow, currentColumn);
    //                 }
    //
    //                 currentColumn = B -> coords[Bindex].j;
    //                 currentColumnIndexB = Bindex;
    //                 endOfColumnIndexB = currentColumnIndexB;
    //             } else if (B -> coords[Bindex].j == currentColumn) {
    //                 // increment index towards end of column
    //                 endOfColumnIndexB = Bindex;
    //             }
    //         }
    //         // move to next row
    //         currentRow = A->coords[Aindex].i;
    //         currentRowIndexA = Aindex;
    //         endOfRowIndexA = currentRowIndexA;
    //
    //     } else if (A -> coords[Aindex].i == currentRow) {
    //         // increment index towards end of row
    //         endOfRowIndexA = Aindex;
    //     }
    // }







    //--------------------------------------------------------------------------
    // Iterate CSR/CSC
    //--------------------------------------------------------------------------
    // int AM_counter;
    // // iterate over rows;
    // int row;
    // int col;
    // double val;
    // for (AM_counter = 0; AM_counter < A ->m; AM_counter ++) {
    //     row = AM_counter;
    //     int start = Ams[AM_counter];
    //     int end = Ams[AM_counter + 1]-1;
    //     // printf("row %i has NNZ index %i to %i\n", AM_counter, start, end);
    //
    //     // All Rows
    //     int rowIndex;
    //     for (rowIndex = start; rowIndex <= end; rowIndex++ ) {
    //         col = A -> coords[rowIndex].j;
    //         val = A -> data[rowIndex];
    //         // printf("(%i, %i) has value %f (A->NZ index = %i)\n", row, col, val, rowIndex);
    //     }
    //
    // }
    // free(Ams);



    // printf("\nAllocated: %i, Actual: %i\n", outputNonZeroes, outputNNZ);
    // printf("multiplied\n");
    // printf("output matrix has dimensions %i x %i\n", (*C)->m, (*C)->n);
    // printf("Output Matrix has dimensions %d x %d, and has %d Non-Zero Entries.\n", A -> m, B -> n, outputNNZ);
}



int getNonZeroes(const COO A, const COO B, COO *C) {
    int outputNNZ = 0;
    // record current row within A
    int currentRow = A->coords[0].i;

    // index in A of start of row
    int currentRowIndexA = 0;

    // index in A of end of row
    int endOfRowIndexA = 0;

    // iterate through all values within A
    int Aindex;
    for (Aindex = 0; Aindex < A->NZ; Aindex ++) {

        if ((A -> coords[Aindex].i != currentRow) || (Aindex == A -> NZ - 1)) {
            // found end of row
            if (Aindex == A -> NZ - 1) {
                endOfRowIndexA = Aindex;
            }

            // iterate over all columns of matrix B
            // record current column
            int currentColumn = B -> coords[0].j;

            // index in B of start of column
            int currentColumnIndexB = 0;

            // index in B of end of column
            int endOfColumnIndexB = 0;

            int Bindex;
            for (Bindex = 0; Bindex < B -> NZ; Bindex ++) {
                if (B -> coords[Bindex].j != currentColumn || Bindex == B -> NZ - 1) {
                    // found end of column
                    if (Bindex == B -> NZ - 1) {
                        endOfColumnIndexB = Bindex;
                    }

                    int addedVal = 0;
                    // iterate between start of row and end of row for column in B
                    int inRowIndexA;
                    for (inRowIndexA = currentRowIndexA; inRowIndexA <= endOfRowIndexA; inRowIndexA++) {
                        // for each non zero value within the column, multiply with respective row elements and set respective row element within C

                        int aJValue = A->coords[inRowIndexA].j;
                        int bIValue;

                        int inColumnIndexB;
                        for (inColumnIndexB = currentColumnIndexB; inColumnIndexB <= endOfColumnIndexB; inColumnIndexB++) {
                            bIValue = B->coords[inColumnIndexB].i;

                            if (aJValue == bIValue) {
                                addedVal = 1;
                            }
                        }
                    }

                    if (addedVal == 1) {
                        outputNNZ ++;
                    }

                    currentColumn = B -> coords[Bindex].j;
                    currentColumnIndexB = Bindex;
                    endOfColumnIndexB = currentColumnIndexB;
                } else if (B -> coords[Bindex].j == currentColumn) {
                    // increment index towards end of column
                    endOfColumnIndexB = Bindex;
                }
            }

            // move to next row
            currentRow = A->coords[Aindex].i;
            currentRowIndexA = Aindex;
            endOfRowIndexA = currentRowIndexA;
        } else if (A -> coords[Aindex].i == currentRow) {
            // increment index towards end of row
            endOfRowIndexA = Aindex;
        }
    }

    return outputNNZ;
}





void addThreeMatrices(COO m1, COO m2, COO m3, COO out)
{
    int m1index = 0;
    int m2index = 0;
    int m3index = 0;
    int outNZ = 0;

    int m1I;
    int m1J;
    int m2I;
    int m2J;
    int m3I;
    int m3J;

    int minI;
    int minJ;

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
double sum_array(double* arr) {
    double result = 0.0;
    result = arr[0] + arr[1] + arr[2];
    return result;
}


void add_three_matrices2(COO m1, COO m2, COO m3, COO* out) {
    int m1index = 0;
    int m2index = 0;
    int m3index = 0;
    int outNZ = 0;

    int m1I;
    int m1J;
    int m2I;
    int m2J;
    int m3I;
    int m3J;

    int minI;
    int minJ;

    // Passover to find indices of rows and columns to add
    // object that has 3 arrays are of size outNZ
    // array of output i's
    // array of output j's
    // array vals which contains 3array of vals
    // 3 array of vals contains the m1, m2, m3 nz or zeroes
    // sum of this 3array becomes out -> data[nz]
    // initially give these arrays size of biggest NZ of input matrix
    int maxNZ = 0;
    maxNZ += m1 -> NZ;
    maxNZ += m2 -> NZ;
    maxNZ += m3 -> NZ;

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
    while (m1index < m1 -> NZ || m2index < m2 -> NZ || m3index < m3 -> NZ) {
        minI = m1 -> m;
        minJ = m1 -> n;
        m1I = minI;
        m1J = minJ;
        m2I = minI;
        m2J = minJ;
        m3I = minI;
        m3J = minJ;

        if (m1index < m1 -> NZ){
          m1I = m1 -> coords[m1index].i;
          m1J = m1 -> coords[m1index].j;
          if (m1I <= minI) {
            minI = m1I;
          }
        }

        if (m2index < m2 -> NZ) {
          m2I = m2 -> coords[m2index].i;
          m2J = m2 -> coords[m2index].j;
          if (m2I <= minI) {
            minI = m2I;
          }
        }

        if (m3index < m3 -> NZ) {
          m3I = m3 -> coords[m3index].i;
          m3J = m3 -> coords[m3index].j;
          if (m3I <= minI) {
            minI = m3I;
          }
        }

        if ((m1I == minI) && (m1J <= minJ)) minJ = m1J;
        if ((m2I == minI) && (m2J <= minJ)) minJ = m2J;
        if ((m3I == minI) && (m3J <= minJ)) minJ = m3J;

        outIs[outNZ] = minI;
        outJs[outNZ] = minJ;

        if (m1index < m1 -> NZ) {
          if (m1I == minI && m1J == minJ) {
            additions[outNZ][0] = m1->data[m1index];
            m1index ++;
          }
        }

        if (m2index < m2 -> NZ) {
          if (m2I == minI && m2J == minJ) {
            additions[outNZ][1] = m2->data[m2index];
            m2index ++;
          }
        }

        if (m3index < m3 -> NZ) {
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
    for (index = 0; index < outNZ; index ++) {
        (*out) -> coords[index].i = outIs[index];
        (*out) -> coords[index].j = outJs[index];
        (*out) -> data[index] = sum_array(additions[index]);
    }

    int addition_index;
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
                            COO *O)
{

    // return basic_sparsemm_sum(A, B, C, D, E, F, O);

    // generate skeleton sparse matrices for G&H which will hold the sum of (A+B+C) and (D+E+F) respectively
    COO G;
    COO H;

    // printf("m1 NZ: %d, m2 NZ: %i, m3 NZ: %i\n", A -> NZ, B -> NZ, C -> NZ);

    // addThreeMatrices(A, B, C, G);
    // addThreeMatrices(D, E, F, H);
    add_three_matrices2(A, B, C, &G);
    add_three_matrices2(D, E, F, &H);

    // call optimised_sparsemm on G & H to output O
    optimised_sparsemm(G, H, O);
    free_sparse(&G);
    free_sparse(&H);
}
