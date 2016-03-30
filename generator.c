//
//  generator.c
//  matrix1
//
//  Created by Jonathan Hattab on 26/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#include "generator.h"

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
double * getPageRankMatrix(int dim){
    double *matrix = NULL;
    matrix = malloc(dim * dim * sizeof(double));
    if (matrix == NULL) exit(0);
    
    for(int i = 0; i < dim; i++){
        int total = 0;
        for(int j = 0; j < dim; j++){
            double coeff = rand() % 11;
            if(coeff < 7) matrix[i*dim + j] = 0;
            else{
                matrix[i*dim + j] = 1;
                total++;
            }
        }
        if(total == 0) continue;
        for(int j = 0; j < dim; j++) matrix[i*dim + j] /= total;
        
    }
    return matrix;
}