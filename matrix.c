//
//  matrix.c
//  matrix1
//
//  Created by Jonathan Hattab on 26/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#include "matrix.h"

// CUT AND ROTATE MATRIX
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
double * submatrixSubvectorProduct(double * submatrix, double * subVector, int dim, int bandWidth){
    double * result = NULL;
    result = malloc(dim * sizeof(double));
    if (result == NULL) exit(0);
    
    for(int i = 0; i < dim; i++){
        result[i] = 0;
    }
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < bandWidth; j++){
            result[i] += submatrix[i + j * dim] * subVector[j];
        }
    }
    return result;
}
double * submatrixSubmatrixProduct(double * submatrix1, double * submatrix2, int dim, int bandWidth){
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
            }
        }
    }
    return result;
}

// MERGING
void mergeMatrixes(double * matrix1, double * matrix2, int nbRows, int nbColumns){
    for(int i = 0; i < nbColumns; i++){
        for(int j = 0; j < nbRows; j++){
            matrix1[i*nbRows + j] += matrix2[i*nbRows + j];
        }
    }
}

// NORME
double computeError(double * vector1, double * vector2, int dim){
    double error = 0;
    for(int i = 0; i < dim; i++){
        error += (vector1[i] - vector2[i]) * (vector1[i] - vector2[i]);
    }
    
    return error;
}

// PROPAGATE RESULT
void propagateResult(int step, int myrank, int size, double * local_result, int nbRows, int nbColumns, MPI_Status * status, int startStep){
    if(myrank % int_pow(2, step) == 0){
        if(myrank + int_pow(2, step-1) < size){
            double * received_result = NULL;
            received_result = malloc(nbRows * nbColumns * sizeof(double));
            if(received_result == NULL) exit(0);
            
            MPI_Recv(received_result, nbRows * nbColumns, MPI_DOUBLE, myrank + int_pow(2, step-1), step +  startStep, MPI_COMM_WORLD, status);
            
            mergeMatrixes(local_result, received_result, nbRows, nbColumns);
            
            free(received_result);
        }
        
        if(int_pow(2, step) < size){
            return propagateResult(step + 1, myrank, size, local_result, nbRows, nbColumns, status, startStep);
        }
        return;
    }
    else{
        MPI_Send(local_result, nbRows * nbColumns, MPI_DOUBLE, myrank - int_pow(2, step-1), step + startStep, MPI_COMM_WORLD);
        free(local_result);
        return;
    }
}