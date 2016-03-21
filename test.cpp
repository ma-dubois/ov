#include <stdio.h>

extern "C"{
#include "ov.h"
}

#define AUI 1.889644746

int main()
{
    real r[49];
    int i,j;


    calc_R_overlap(r, 0,0,1,
    2, 1, 1., 0, 1, 0,
    3, 2, 1., 0, 1, 0 );

    for(i=0;i<3;i++) {
        for(j=0;j<5;j++)
        printf("%.4f ", r[i*5+j]);
        printf("\n");
    }
    
    return 0;
}
