//
//  simulation_Xr.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <string>

#include "GTAcpp_selfish.h"

using namespace std;




int weightedSampling(long *rnum1, int *reactionArray){ // this reaction arry should include all the reactions
    int x=0;
    int j=0;
    long temp = 0;
    for (j=0; j < 4; j++){
        if (*rnum1 > temp && *rnum1 <= temp + reactionArray [j]) x =j;
        temp = temp + reactionArray[j];
    }
    return x;
}


double calcTau(int *reactionArray, double *rnum2){
    double tau = 0; // for update the time;
    double sum =  0;
    int i=0;
    for (i =0; i < 4; i++){
        sum = sum + (double)reactionArray[i];
    }
    //printf("sum of array is %lf\n", sum);
    //printf("log sampel is %lf\n", log(sample));
    if(*rnum2 == 0){
        tau = (-log(*rnum2+0.0001))/(sum);
    }
    else{
        tau = (-log(*rnum2))/(sum);
    }
    return tau;
}


// determine logic condition for while loop
// no popluation size is negative
int isNegativePopSize(long *popSizeArray){ // this pop size array should include all the populations
    int a=0;
    //int x=0;
    int positivePopSize= 0;
    if (popSizeArray[0] < 0 || popSizeArray[1] < 0){
    		return (0);
    }
    else {
    		return (1);
    }
}

// no reaction is negative
int isNegativeReaction(int *reactionArray){
    int c=0;
    //int y;
    int positiveReaction = 0;
    for (c=0; c < 4; c++){
        if (reactionArray[c] >=0) positiveReaction = positiveReaction +1;
    }
    if (positiveReaction == 4) {
        //y =1;
        //printf("%d\n", positiveReaction);
        //printf("Y = %d\n", y);
        return (1);
    }
    else {
        return (0);
    }
}
