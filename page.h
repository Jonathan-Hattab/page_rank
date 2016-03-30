//
//  page.h
//  matrix1
//
//  Created by Jonathan Hattab on 27/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#ifndef page_h
#define page_h

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "tools.h"

struct Page{
    int page_index;
    double score;
};

void printRanking(struct Page * ranking, int dim);
void getPageDataType(MPI_Datatype * page_type);
struct Page * sortRanking(struct Page * ranking, int dim);
struct Page * mergeRanking(struct Page * ranking1, struct Page * ranking2, int size1, int size2);

struct Page * propagateRanking(int step, int myrank, int size, struct Page * local_ranking, int size1, MPI_Status * status, int startStep);

#endif /* page_h */
