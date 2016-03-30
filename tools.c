//
//  tools.c
//  matrix1
//
//  Created by Jonathan Hattab on 28/03/2016.
//  Copyright © 2016 Jonathan Hattab. All rights reserved.
//

#include "tools.h"

// INT POWER
int int_pow(int base, int exp){
    int result = 1;
    for(int i = 0; i < exp; i++) result *= base;
    return result;
}