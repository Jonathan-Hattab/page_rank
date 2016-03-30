//
//  main.c
//  matrix1
//
//  Created by Jonathan Hattab on 17/02/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "generator.h"
#include "matrix.h"
#include "page.h"

int main(int argc, char * argv[]) {
    // The Dimension of the original square matrices
    int dim = 30;
    
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
    int local_bandWidth = 0;
    
    /* ========================= PREPARATION DE BASE DE L'EXERCICE =========================== //
        Le Processeur n°0 :
            Crée la matrice de page rank.
            Découpe la matrice en bandes verticales ± équitables.
            Envoie la largeurs des bandes et les bandes aux autres processeurs.
        Les autres processeurs :
            Reçoivent la largeur de leur bande et leur bande.
    */
    if (myrank == 0){
        srand((unsigned)time(NULL)); // Initialization for the random function
        
        double * matrix = getPageRankMatrix(dim);
        
        // Uncomment to view the page rank matrix
        
        printf("Page Rank Matrix (see it transposed) = \n");
        printMatrix(matrix, dim, dim);
        
        
        int startRow = 0;
        
        for(int bandNumber = 0; bandNumber < size; bandNumber++){
            int rowsLeft = dim - startRow;

            int bandWidth = ceil((double)rowsLeft / (double)(size - bandNumber));
            
            double * columnBand = getBand(matrix, dim, dim, startRow, bandWidth);
            
            // Uncomment to view the coresponding band of the matrix and the associated sub-vector
            /*
            printf("\nSub-Division n°%d\n", bandNumber);
            printf("\tColumn Band : \n");
            printMatrix(columnBand, bandWidth, dim);
            */

            if(bandNumber == 0){
                local_bandWidth     = bandWidth;
                local_column_band   = columnBand;
            }
            else{
                MPI_Send(&bandWidth,   1,              MPI_INT,    bandNumber, 1, MPI_COMM_WORLD);
                MPI_Send(columnBand,  bandWidth * dim, MPI_DOUBLE, bandNumber, 2, MPI_COMM_WORLD);
                
                free(columnBand);
            }
            
            startRow += bandWidth;
        }
        
        free(matrix);
    }
    else{
        MPI_Recv(&local_bandWidth, 1,                        MPI_INT,    0, 1, MPI_COMM_WORLD, &status);
        
        local_column_band   = malloc(dim * local_bandWidth * sizeof(double));
        if(local_column_band == NULL) exit(0);
        MPI_Recv(local_column_band, local_bandWidth * dim,   MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
    }
    
    // ======================================= EXERCICE ====================================== //
    // Uncomments to view the local band, the local vector and the full vector stored in the processor 0
    /*
    printf("\n[%d] Local Column Band : \n", myrank);
    printMatrix(local_column_band, local_bandWidth, dim);
    */
    
    double * result = NULL;

    if(myrank == 0){
        result = malloc(dim * sizeof(double));
        if(result == NULL) exit(0);
        for(int i = 0; i < dim; i++) result[i] = 1;
    }
    
    double * local_subVector = NULL;
    
    double precision = 0.001;
    double error = 0;
    int iteration = 0;
    int maxIterations = 100;
    
    // Approach solution with iterative method
    do{
        iteration++;
        
        if(myrank == 0){
            int external_bandWidth = 0;
            int startRow = 0;
            
            for(int i = 0; i < size; i++){
                if(i == myrank) external_bandWidth = local_bandWidth;
                else{
                    MPI_Recv(&external_bandWidth, 1, MPI_INT, i, 2 + i, MPI_COMM_WORLD, &status);
                }
                
                double * subVector = getBand(result, dim, 1, startRow, external_bandWidth);
                
                if(i == myrank){
                    local_subVector = subVector;
                }
                else{
                    MPI_Send(subVector, external_bandWidth, MPI_DOUBLE, i, 2 + i + size, MPI_COMM_WORLD);
                    free(subVector);
                }
                
                startRow += external_bandWidth;
            }
        }else{
            MPI_Send(&local_bandWidth, 1,                  MPI_INT,    0, 2 + myrank,          MPI_COMM_WORLD);
            
            local_subVector = malloc(local_bandWidth * sizeof(double));
            if(local_subVector == NULL) exit(0);
            MPI_Recv(local_subVector,   local_bandWidth,   MPI_DOUBLE, 0, 2 + myrank + size,   MPI_COMM_WORLD, &status);
        }
        
        // PHASE DE CALCUL
        double * previousResult = result;
        
        result = submatrixSubvectorProduct(local_column_band, local_subVector, dim, local_bandWidth);
        /*
        printf("\nSub Result n°%d at iteration n°%d :\n", myrank, iteration);
        printMatrix(result, dim, 1);
        fflush(stdout);
        */
        propagateResult(1, myrank, size, result, dim, 1, &status, 2 + 2 * size);
        
        if(myrank == 0){
            error = computeError(result, previousResult, dim);
            free(previousResult);
        }
        
        MPI_Bcast(&error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if(myrank == 0){
            /*
            printf("\nResult at iteration n°%d :\n", iteration);
            printMatrix(result, dim, 1);
            printf("\nError at iteration n°%d : %f\n", iteration, error);
            fflush(stdout);
            */
        }
    }while(error > precision  && iteration <= maxIterations);
    
    free(local_column_band);
    
    // Sort ranking in parallel
    MPI_Datatype page_type;
    getPageDataType(&page_type);
    struct Page * local_ranking = NULL;
    
    if(myrank == 0){
        // Transform score list to Page list
        struct Page * ranking = NULL;
        ranking = malloc(dim * (sizeof(struct Page)));
        if(ranking == NULL) exit(0);
        
        for(int i = 0; i < dim; i++){
            ranking[i].score = result[i];
            ranking[i].page_index = i;
        }
        
        free(result);
        
        int external_bandWidth = 0;
        int startRow = 0;
        
        for(int i = 0; i < size; i++){
            if(i == myrank) external_bandWidth = local_bandWidth;
            else{
                MPI_Recv(&external_bandWidth, 1, MPI_INT, i, 2 + i, MPI_COMM_WORLD, &status);
            }
            
            struct Page * subRanking = &ranking[startRow];
            
            if(i == myrank){
                local_ranking = subRanking;
            }
            else{
                MPI_Send(subRanking, external_bandWidth, page_type, i, 2 + i + size, MPI_COMM_WORLD);
            }
            
            startRow += external_bandWidth;
        }
    }else{
        MPI_Send(&local_bandWidth, 1, MPI_INT, 0, 2 + myrank, MPI_COMM_WORLD);
        
        local_ranking = malloc(local_bandWidth * sizeof(struct Page));
        if(local_ranking == NULL) exit(0);
        MPI_Recv(local_ranking, local_bandWidth, page_type, 0, 2 + myrank + size, MPI_COMM_WORLD, &status);
    }
    
    MPI_Type_free(&page_type);
    
    struct Page * sortedRanking = sortRanking(local_ranking, local_bandWidth);
    
    /*
    printf("\n[%d] Finished Computing Local Result : \n", myrank);
    printRanking(sortedRanking, local_bandWidth);
    */
     
    sortedRanking = propagateRanking(1, myrank, size, sortedRanking, local_bandWidth, &status, 2 + 2 * size);
    
    if(myrank == 0){
        printf("\nFinished Computing Result : \n");
        printRanking(sortedRanking, dim);
        
        free(sortedRanking);
    }
    
    MPI_Finalize();
    return 0;
}
