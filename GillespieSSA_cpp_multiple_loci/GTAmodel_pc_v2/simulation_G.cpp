//
//  simulation_Xr.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <string>
#include "GTAcpp_G.h"

using namespace std;

int reactionSystem(string s){
    if (s.compare("Xa") ==0) {return 1;}
    else if (s.compare("Xn") ==0) {return 2;}
    else if (s.compare("Xr") ==0) {return 3;}
    else if (s.compare("XaXn") ==0) {return 4;}
    else if (s.compare("XaXr") ==0) {return 5;}
    else if (s.compare("XnXr") ==0) {return 6;}
    else if (s.compare("XaXnXr") ==0) {return 7;}
    else return 0;
}

int matrixIndex(int rowNum, int colNum){ // reaction status of lysis influence Xn and Xr pop size
    int i =0;
    int j =0;
    int k =0;
    int n =0;
    int x =0;
    int y =0;
    int l1 =0;
    int l2 =0;
    int index =0;
    int times =0;
    struct matrix {
        int rowNum[MAX_genotypeNum];
    };
    struct matrix geneMatrix[NUM_LOCI];
    for (i =0; i <NUM_LOCI; i++){
        l1 = (int)pow(2,i);
        l2 = (int)pow(2,(i+1));
        int seq1[l1];
        int seq2[l1];
        int seq3[l2];
        j = 2*i+1;
        k = 2*i+2;
        for(n=0; n < l1; n++){
            seq1[n] = j;
            seq2[n] = k;
            seq3[n] = seq1[n];
            seq3[(l1+n)] = seq2[n];
        }
        times = (int)(pow(2, NUM_LOCI)/pow(2,(i+1)));
        for (x=0; x < times;x++){
            for (y=0; y < l2; y++){
                geneMatrix[i].rowNum[(x*l2+y)] = seq3[y];
            }
        }
        
    }
    /*int a,b;
    for (a =0; a < NUM_LOCI; a++){

    		for (b = 0; b < (times*l2+l2); b++){
    			printf("%d\t",geneMatrix[a].rowNum[b]);
    		}
    		printf("\n");
    }
    */
    index = geneMatrix[rowNum].rowNum[colNum];
    return index;
}


int weightedSampling(long *rnum1, int *reactionArray){ // this reaction arry should include all the reactions
    int x=0;
    int j=0;
    long temp = 0;
    for (j=0; j < (numOfgeonReact+numofgeneReact); j++){
        if (*rnum1 > temp && *rnum1 <= temp + reactionArray [j]) x =j;
        temp = temp + reactionArray[j];
    }
    return x;
}

double calcTau(int *reactionArray, double *rnum2){
    double tau = 0; // for update the time;
    double sum =  0;
    int i=0;
    for (i =0; i < (numOfgeonReact+numofgeneReact); i++){
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
    for (a=0; a < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); a++){
        if (popSizeArray[a] >= 0) {
            positivePopSize = positivePopSize +1;
        }
    }
    
    if (positivePopSize == (genotypeNum*NUM_OF_POPTYPE+genetypeNum)){
        //x =1;
        //printf("%d\n", positivePopSize);
        //printf("X = %d\n", x);
        return (1);
    }
    else {
        return (0);
    }
}

// no reaction is negative
int isNegativeReaction(int *reactionArray){
    int c=0;
    //int y;
    int positiveReaction = 0;
    for (c=0; c < (numOfgeonReact+numofgeneReact); c++){
        if (reactionArray[c] >=0) positiveReaction = positiveReaction +1;
    }
    if (positiveReaction == (numOfgeonReact+numofgeneReact)) {
        //y =1;
        //printf("%d\n", positiveReaction);
        //printf("Y = %d\n", y);
        return (1);
    }
    else {
        return (0);
    }
}
