//
//  page.c
//  matrix1
//
//  Created by Jonathan Hattab on 27/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#include "page.h"

// PRINT RANKING
void printRanking(struct Page * ranking, int dim){
    for(int i = 0; i < dim; i++){
        printf("(%d)\t", (i + 1));
        printf("%f\t", ranking[i].score);
        printf("Page n°%d", ranking[i].page_index);
        printf("\n");
    }
}


// MPI DATA TYPE
void getPageDataType(MPI_Datatype * page_type){
    struct Page page;
    page.page_index = 1;
    page.score = 1.0;
    
    int blocks[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
    
    MPI_Aint baseAddr, pageIndexAddr, scoreAddr;
    MPI_Get_address(&page,            &baseAddr);
    MPI_Get_address(&page.page_index, &pageIndexAddr);
    MPI_Get_address(&page.score,      &scoreAddr);
    
    MPI_Aint displacements[2];
    displacements[0] = pageIndexAddr - baseAddr;
    displacements[1] = scoreAddr - baseAddr;
    MPI_Type_create_struct(2, blocks, displacements, types, page_type);
    MPI_Type_commit(page_type);
}

//SORT
struct Page * sortRanking(struct Page * ranking, int dim){
    if(dim == 1) return ranking;
    else{
        int half = (int) dim / 2.0;
        //printf("\nPart 1 : \n");
        //printRanking(ranking, half);
        //printf("\nPart 2 : \n");
        //printRanking(&(ranking[half]), dim - half);
        
        struct Page * ranking1 = sortRanking(ranking, half);
        struct Page * ranking2 = sortRanking(&(ranking[half]), dim - half);
        struct Page * result = mergeRanking(ranking1, ranking2, half, dim - half);
        return result;
    }
}
struct Page * mergeRanking(struct Page * ranking1, struct Page * ranking2, int size1, int size2){
    struct Page * result = NULL;
    result = malloc((size1 + size2) * sizeof(struct Page));
    if(result == NULL) exit(0);
    
    int j = 0;
    int index = 0;
    for(int i = 0; i < size1; i++){
        while(j < size2 && ranking2[j].score > ranking1[i].score){
            result[index] = ranking2[j];
            index++;
            j++;
        }
        result[index] = ranking1[i];
        index++;
    }
    while(j < size2){
        result[index] = ranking2[j];
        index++;
        j++;
    }
    
    /*
     printf("\nPart 1\n");
     printRanking(ranking1, size1);
     
     printf("\nPart 2\n");
     printRanking(ranking2, size2);
     
     printf("\n Result\n");
     printRanking(result, size1 + size2);
     */
    
    if(size1 > 1) free(ranking1);
    if(size2 > 1) free(ranking2);
    
    return result;
}

struct Page * propagateRanking(int step, int myrank, int size, struct Page * local_ranking, int size1, MPI_Status * status, int startStep){
    MPI_Datatype page_type;
    getPageDataType(&page_type);
    
    if(myrank % int_pow(2, step) == 0){
        int size2 = 0;
        if(myrank + int_pow(2, step-1) < size){
            MPI_Recv(&size2, 1, MPI_INT, myrank + int_pow(2, step-1), step + startStep, MPI_COMM_WORLD, status);
            
            struct Page * external_ranking = NULL;
            external_ranking = malloc(size2 * sizeof(struct Page));
            if(external_ranking == NULL) exit(0);
            
            MPI_Recv(external_ranking, size2, page_type, myrank + int_pow(2, step-1), step +  startStep + size, MPI_COMM_WORLD, status);
            MPI_Type_free(&page_type);
            
            local_ranking = mergeRanking(local_ranking, external_ranking, size1, size2);
        }
        
        if(int_pow(2, step) < size){
            return propagateRanking(step + 1, myrank, size, local_ranking, size1 + size2, status, startStep);
        }
        return local_ranking;
    }
    else{
        MPI_Send(&size1,        1,      MPI_INT,    myrank - int_pow(2, step-1), step + startStep,          MPI_COMM_WORLD);
        MPI_Send(local_ranking, size1,  page_type,  myrank - int_pow(2, step-1), step + startStep + size,   MPI_COMM_WORLD);
        free(local_ranking);
        MPI_Type_free(&page_type);
        return NULL;
    }
}