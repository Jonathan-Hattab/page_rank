//
//  generator.h
//  matrix1
//
//  Created by Jonathan Hattab on 26/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#ifndef generator_h
#define generator_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double * getMatrix(int dim);
double * getVector(int dim);
double * getPageRankMatrix(int dim);

#endif // generator_h
