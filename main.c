//
//  main.c
//  matrix1
//
//  Created by Jonathan Hattab on 17/02/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// GET MATRIX & BANDS
double * getMatrix(int dim){
    double *matrix = NULL;
    matrix = malloc(dim * dim * sizeof(double));
    if (matrix == NULL) exit(0);
    
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
//            double coeff = rand() % 11;
            double coeff = i*dim + j;
            matrix[i*dim + j] = coeff;
        }
    }
    return matrix;
}
double * getVector(int dim){
    double *vector = NULL;
    vector = malloc(dim * sizeof(double));
    if (vector == NULL) exit(0);
    
    for(int i = 0; i < dim; i++){
        double coeff = rand() % 11;
        // double coeff = i;
        vector[i] = coeff;
    }
    return vector;
}
double * getBand(double * matrix, int nbRows, int nbColumns, int startRow, int rowsInBand){
    double * band = NULL;
    band = malloc(rowsInBand * nbColumns * sizeof(double));
    if (band == NULL) exit(0);
    
    for(int i = 0; i < rowsInBand; i++){
        for(int j = 0; j < nbColumns; j++){
            band[i * nbColumns + j] = matrix[(i + startRow) * nbColumns + j];
        }
    }
    return band;
}
double * rotateMatrix(double * matrix, int dim){
    double * rotatedMatrix = NULL;
    rotatedMatrix = malloc(dim * dim * sizeof(double));
    if (rotatedMatrix == NULL) exit(0);
    
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            rotatedMatrix[i*dim + j] = matrix[j*dim+i];
        }
    }
    return rotatedMatrix;
}

// PRINT MATRIX
void printMatrix(double * matrix, int nbRows, int nbColumns){
    for(int i=0; i<nbRows; i++){
        for(int j=0; j<nbColumns; j++){
            printf("%f\t", matrix[i*nbColumns + j]);
        }
        printf("\n");
    }
}

// COMPUTE ROW * VECTOR PRODUCT
double rowVectorProduct(double * row, double * vector, int nbColumns){
    double result = 0.0;
    for(int i = 0; i < nbColumns; i++) result += row[i]*vector[i];
    return result;
}
double * matrixVectorProduct(double * matrix, double * vector, int nbRows, int nbColumns){
    double *result = NULL;
    result = malloc(nbRows * sizeof(double));
    if (result == NULL) exit(0);
    
    for(int i = 0; i < nbRows; i++){
        result[i] = 0;
        for(int j = 0; j < nbColumns; j++){
            result[i] += matrix[i*nbColumns + j]*vector[j];
        }
    }
    
    return result;
}
double * submatrixSubvectorProduct(double * submatrix, double * subVector, int nbRows, int nbColumns){
    double * result = NULL;
    result = malloc(nbRows * sizeof(double));
    if (result == NULL) exit(0);
    
    for(int j = 0; j < nbRows; j++){
        result[j] = 0;
    }
    for(int i = 0; i < nbColumns; i++){
        for(int j = 0; j < nbRows; j++){
            result[i] += submatrix[i*nbRows + j] * subVector[i];
        }
    }
    return result;
}
double * submatrixSubmatrixProduct(double * submatrix1, double * submatrix2, int dim, int bandWidth, int myrank){
    double * result = NULL;
    result = malloc(dim * dim * sizeof(double));
    if (result == NULL) exit(0);
    
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            result[i*dim + j] = 0;
        }
    }
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            for(int k = 0; k < bandWidth; k++){
                result[i*dim + j] += submatrix1[j + k * dim] * submatrix2[i + k * dim];
                if(myrank == 0){
                    //printf("\n(%d, %d, (%d)) : %f * %f = %f", i, j, k, submatrix1[j*dim+k], submatrix2[i*dim + k], result[i*dim + j]);
                }
            }
        }
    }
    return result;
}

// INT POWER
int int_pow(int base, int exp){
    int result = 1;
    for(int i = 0; i < exp; i++) result *= base;
    return result;
}

// MERGING
void mergeMatrixes(double * matrix1, double * matrix2, int dim){
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            matrix1[i*dim + j] += matrix2[i*dim + j];
        }
    }
}

void propagateResult(int step, int myrank, int size, double * local_result, int dim, MPI_Status * status, int startStep){
    if(myrank % int_pow(2, step) == 0){
        if(myrank + int_pow(2, step-1) < size){
            double * received_result = NULL;
            received_result = malloc(dim * dim * sizeof(double));
            if(received_result == NULL) exit(0);
            
            MPI_Recv(received_result, dim * dim, MPI_DOUBLE, myrank + int_pow(2, step-1), step +  startStep, MPI_COMM_WORLD, status);
            
            mergeMatrixes(local_result, received_result, dim);
            
            free(received_result);
        }
        
        if(int_pow(2, step) < size){
            return propagateResult(step + 1, myrank, size, local_result, dim, status, startStep);
        }
        return;
    }
    else{
        MPI_Send(local_result, dim * dim, MPI_DOUBLE, myrank - int_pow(2, step-1), step + startStep, MPI_COMM_WORLD);
        free(local_result);
        return;
    }
}

int main(int argc, char * argv[]) {
    // The Dimension of the original square matrices
    int dim = 10;
    
    int myrank;
    int size;
    
    MPI_Status status;
    MPI_Init (&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
    
    // We can't use more processors than dim
    if(size > dim) size = dim;
    if(myrank >= size){
        MPI_Finalize();
        return 0;
    }
    
    double * local_column_band  = NULL;
    double * local_row_band     = NULL;
    int local_rowsInBand = 0;
    
    /* ========================= PREPARATION DE BASE DE L'EXERCICE =========================== //
        Le Processeur n°0 :
            Crée les deux matrices de travail.
            Découpe les matrice en bandes ± équitables mais les bandes correspondantes de chaque matrice sont de même taille.
            Envoie la hauteur des bandes et les deux bandes aux autres processeurs.
        Les autres processeurs :
            Reçoivent la hauteur de leur bandes et leurs deux bandes.
    */
    if (myrank == 0){
        srand((unsigned)time(NULL));
        
        double * primaryMatrix1 = getMatrix(dim);
        double * matrix2        = getMatrix(dim);
        
        // Uncomment to view the global matrix and the global vector
        /*
        printf("Matrix 1 = \n");
        printMatrix(primaryMatrix1, dim, dim);
        printf("\nMatrix 2 = \n");
        printMatrix(matrix2, dim, dim);
        */
        
        // In order to split the first matrix in column bands using the same function as for row bands, we rotate it.
        double * matrix1 = rotateMatrix(primaryMatrix1, dim);
        free(primaryMatrix1);
         
        int startRow = 0;
        
        for(int bandNumber = 0; bandNumber < size; bandNumber++){
            int rowsLeft = dim - startRow;

            int rowsInBand = ceil((double)rowsLeft / (double)(size - bandNumber));
            
            double * columnBand = getBand(matrix1, dim, dim, startRow, rowsInBand);
            double * rowBand    = getBand(matrix2, dim, dim, startRow, rowsInBand);
            
            // Uncomment to view the coresponding band of the matrix and the associated sub-vector
            /*
            printf("\nSub-Division n°%d\n", bandNumber);
            printf("\tColumn Band : \n");
            printMatrix(columnBand, rowsInBand, dim);
            printf("\tRow Band : \n");
            printMatrix(rowBand, rowsInBand, dim);
            */

            if(bandNumber == 0){
                local_rowsInBand    = rowsInBand;
                local_column_band   = columnBand;
                local_row_band      = rowBand;
            }
            else{
                MPI_Send(&rowsInBand,   1,              MPI_INT,    bandNumber, 1, MPI_COMM_WORLD);
                MPI_Send(columnBand,  rowsInBand * dim, MPI_DOUBLE, bandNumber, 2, MPI_COMM_WORLD);
                MPI_Send(rowBand,     rowsInBand * dim, MPI_DOUBLE, bandNumber, 3, MPI_COMM_WORLD);
                
                free(columnBand);
                free(rowBand);
            }
            
            startRow += rowsInBand;
        }
        
        free(matrix1);
        free(matrix2);
    }
    else{
        MPI_Recv(&local_rowsInBand, 1,                        MPI_INT,    0, 1, MPI_COMM_WORLD, &status);
        
        local_column_band   = malloc(dim * local_rowsInBand * sizeof(double));
        local_row_band      = malloc(dim * local_rowsInBand * sizeof(double));
        if(local_column_band == NULL || local_row_band == NULL) exit(0);
        MPI_Recv(local_column_band, local_rowsInBand * dim,   MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(local_row_band,    local_rowsInBand * dim,   MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
    }
    
    // ======================================= EXERCICE ====================================== //
    // Uncomments to view the local band, the local vector and the full vector stored in the processor 1
    /*
    if(myrank == 1){
         printf("\nLocal Band : \n");
         printMatrix(local_band, local_rowsInBand, dim);
         printf("\nLocal Vector : \n");
         printMatrix(local_vector, local_rowsInBand, 1);
         printf("\nFull Vector : \n");
         printMatrix(fullVector, dim, 1);
    }
     */
    
    // PHASE DE CALCUL
    double * result = submatrixSubmatrixProduct(local_column_band, local_row_band, dim, local_rowsInBand, myrank);

    /*
    if(myrank == 0){
        printf("\nLocal rows in band : %d", local_rowsInBand);
        printf("\nLocal Result : \n");
        printMatrix(result, dim, dim);
    }
    */
    
    
    free(local_column_band);
    free(local_row_band);
    
    // RENVOI DU RESULTAT AU PROCESSEUR 0
    propagateResult(1, myrank, size, result, dim, &status, 3);
    
    if(myrank == 0){
        // The result is transposed so we have to rotate the matrix back.
        double * computedMatrix = NULL;
        computedMatrix = malloc(dim * dim * sizeof(double));
        if(computedMatrix == NULL) exit(0);
        
        computedMatrix = rotateMatrix(result, dim);
        free(result);
        
        // Printing the result
        printf("\nFinished Computing Result : \n");
        printMatrix(computedMatrix, dim, dim);
        free(computedMatrix);
    }
    
    MPI_Finalize();
    return 0;
}
