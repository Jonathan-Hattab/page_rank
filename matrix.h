//
//  matrix.h
//  matrix1
//
//  Created by Jonathan Hattab on 26/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#ifndef matrix_h
#define matrix_h

#include <stdio.h>
#include <stdlib.h>

double * getBand(double * matrix, int nbRows, int nbColumns, int startRow, int rowsInBand);
double * rotateMatrix(double * matrix, int dim);

void printMatrix(double * matrix, int nbRows, int nbColumns);

double rowVectorProduct(double * row, double * vector, int nbColumns);
double * matrixVectorProduct(double * matrix, double * vector, int nbRows, int nbColumns);
double * submatrixSubvectorProduct(double * submatrix, double * subVector, int nbRows, int nbColumns);
double * submatrixSubmatrixProduct(double * submatrix1, double * submatrix2, int dim, int bandWidth, int myrank);

void mergeMatrixes(double * matrix1, double * matrix2, int dim);

#endif /* matrix_h */
