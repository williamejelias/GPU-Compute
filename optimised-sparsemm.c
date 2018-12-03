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
    printf("multiplying\n");

    // return basic_sparsemm(A, B, C);

    sortColumnOrder(B);

    int outputNNZ = 0;
    alloc_sparse(A -> m, B -> n, A -> m * B -> n, C);



    // convert A to CSR form
    // NZs
    // double *ANZ = NULL;
    // ANZ = malloc(A -> NZ *sizeof(*ANZ));
    //
    // // AIs
    // double *AIs = NULL;
    // AIs = malloc((A -> m + 1)*sizeof(*AIs));
    //
    // // AJs
    // double *AJs = NULL;
    // AJs = malloc(A -> NZ *sizeof(*AJs));
    //
    // if(ANZ == NULL || AIs == NULL || AJs == NULL){
    //     fprintf(stderr, "Out of memory...\n");
    //     exit(1);
    // }


    // printf("%i, %i", A->m, B->n);

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

                    double result = 0.0;
                    // iterate between start of row and end of row for column in B
                    int inRowIndexA;
                    for (inRowIndexA = currentRowIndexA; inRowIndexA <= endOfRowIndexA; inRowIndexA++) {
                        // for each non zero value within the column, multiply with respective row elements and set respective row element within C

                        double aData = A -> data[inRowIndexA];
                        int aIValue = A->coords[inRowIndexA].i;
                        int aJValue = A->coords[inRowIndexA].j;

                        double bData;
                        int bIValue;
                        int bJValue;

                        int inColumnIndexB;
                        for (inColumnIndexB = currentColumnIndexB; inColumnIndexB <= endOfColumnIndexB; inColumnIndexB++) {
                            bData = B -> data[inColumnIndexB];
                            bIValue = B->coords[inColumnIndexB].i;
                            bJValue = B->coords[inColumnIndexB].j;

                            if (aJValue == bIValue) {
                                result += aData * bData;
                            }
                        }
                    }

                    if (result != 0.0) {
                        (*C) -> coords[outputNNZ].i = currentRow;
                        (*C) -> coords[outputNNZ].j = currentColumn;
                        (*C) -> data[outputNNZ] = result;
                        outputNNZ ++;
                        // printf("New output NNZ value of %f at coordinates (%d, %d)\n", result, currentRow, currentColumn);
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

    // allocate correct NNZ value
    // printf("\n%i\n", outputNNZ);
    (*C)-> NZ = outputNNZ;
    printf("multiplied\n");
    // printf("Output Matrix has dimensions %d x %d, and has %d Non-Zero Entries.\n", A -> m, B -> n, outputNNZ);


    // free(ANZ);
    // free(AIs);
    // free(AJs);
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
    result = result + arr[0] + arr[1] + arr[2];
    return result;
}


void add_three_matrices2(COO m1, COO m2, COO m3, COO* out) {
    printf("Adding\n");
    int m1index = 0;
    int m2index = 0;
    int m3index = 0;
    int outNZ = 0;

    int m1I = m1 -> m;
    int m1J = m1 -> n;
    int m2I = m2 -> m;
    int m2J = m2 -> n;
    int m3I = m3 -> m;
    int m3J = m3 -> n;

    int minI;
    int minJ;

    // Passover to find indices of rows and columns to add

    // object that has 3 arrays
    // all arrays are of size outNZ
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

    double *outIs = NULL;
    outIs = malloc(maxNZ*sizeof(*outIs));
    if(outIs == NULL){
        fprintf(stderr, "Out of memory");
        exit(1);
    }

    double *outJs = NULL;
    outJs = malloc(maxNZ*sizeof(*outJs));
    if(outJs == NULL){
        fprintf(stderr, "Out of memory");
        exit(1);
    }

    double** additions = malloc(maxNZ * sizeof(double*));
    if(outIs == NULL || outJs == NULL || additions == NULL){
        fprintf(stderr, "Out of memory");
        exit(1);
    }
    int i;
    for(i = 0; i < maxNZ; i++){
        double* p = malloc(3 * sizeof(double));
        if(p == NULL){
            fprintf(stderr, "Out of memory");
            exit(1);
        }
        additions[i] = p;
    }

    // addition
    while (m1index < m1 -> NZ || m2index < m2 -> NZ || m3index < m3 -> NZ) {
        // if (m1index < m1 -> NZ) {
            m1I = m1 -> coords[m1index].i;
            m1J = m1 -> coords[m1index].j;
        // } else {
        //   printf("m1 index at max");
        // }

        // if (m2index < m2 -> NZ) {
            m2I = m2 -> coords[m2index].i;
            m2J = m2 -> coords[m2index].j;
        // } else {
        //   printf("m2 index at max");
        // }

        // if (m3index < m3 -> NZ) {
            m3I = m3 -> coords[m3index].i;
            m3J = m3 -> coords[m3index].j;
        // } else {
        //   printf("m3 index at max");
        // }

        minI = m1 -> m;
        minJ = m1 -> n;

        printf("\nM1 %d, %d\n", m1index, m1 -> NZ);
        printf("M2 %d, %d\n", m2index, m2 -> NZ);
        printf("M3 %d, %d\n", m3index, m3 -> NZ);

        if (m1index < m1 -> NZ) if (m1I < minI) minI = m1 -> coords[m1index].i;
        if (m2index < m2 -> NZ) if (m2I < minI) minI = m2 -> coords[m2index].i;
        if (m3index < m3 -> NZ) if (m3I < minI) minI = m3 -> coords[m3index].i;

        if ((m1I == minI) && (m1J < minJ)) minJ = m1 -> coords[m1index].j;
        if ((m2I == minI) && (m2J < minJ)) minJ = m2 -> coords[m2index].j;
        if ((m3I == minI) && (m3J < minJ)) minJ = m3 -> coords[m3index].j;

        outIs[outNZ] = minI;
        outJs[outNZ] = minJ;

        if (m1I == minI && m1J == minJ) {
            additions[outNZ][0] = m1->data[m1index];
            m1index ++;
        }

        if (m2I == minI && m2J == minJ) {
            additions[outNZ][1] = m2->data[m2index];
            m2index ++;
        }

        if (m3I == minI && m3J == minJ) {
            additions[outNZ][2] = m3->data[m3index];
            m3index ++;
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

    for(i = 0; i < maxNZ; i++){
        free(additions[i]);
    }
    free(additions);
    free(outIs);
    free(outJs);
    printf("Added\n");
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

    // addThreeMatrices(A, B, C, G);
    // addThreeMatrices(D, E, F, H);
    add_three_matrices2(A, B, C, &G);
    add_three_matrices2(D, E, F, &H);

    // call optimised_sparsemm on G & H to output O
    optimised_sparsemm(G, H, O);
    free_sparse(&G);
    free_sparse(&H);
}
