//
//  simulation_Xr.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <string>

#include "GTAcpp_recomb.h"

using namespace std;

int reactionSystem(string s){
    if (s.compare("Xa") ==0) {return 1;}
    else if (s.compare("Xn") ==0) {return 2;}
    else if (s.compare("Xr") ==0) {return 3;}
    else if (s.compare("Xnr") ==0) {return 4;}
    else if (s.compare("XaXn") ==0) {return 5;}
    else if (s.compare("XaXr") ==0) {return 6;}
    else if (s.compare("XaXnr") ==0) {return 7;}
    else if (s.compare("XnXr") ==0) {return 8;}
    else if (s.compare("XnXnr") ==0) {return 9;}
    else if (s.compare("XrXnr") ==0) {return 10;}
    else if (s.compare("XaXnXr") ==0) {return 11;}
    else if (s.compare("XaXnXnr") ==0) {return 12;}
    else if (s.compare("XaXrXnr") ==0) {return 13;}
    else if (s.compare("XnXrXnr") ==0) {return 14;}
    else if (s.compare("XaXnXrXnr") ==0) {return 15;}
    else return 0;
}


int *initPopSize(int totalPopSize){

	static int popSizeArray[8];

	int a;
	for(a=0; a<8; a++){
		popSizeArray[a] = 0;
	}
	if (NUM_STRAIN == 1){

	    popSizeArray[0] = totalPopSize*PROP_STRAIN1*(1 - PROP_WILDTYPE);
	    popSizeArray[1] = totalPopSize*PROP_STRAIN1*PROP_WILDTYPE;
	}

	if (NUM_STRAIN == 2){

	    if (PROP_STRAIN1==0){
	    	PROP_STRAIN1 = 1/((double)NUM_STRAIN);
	    }
	    popSizeArray[0] = totalPopSize*PROP_STRAIN1*(1 - PROP_WILDTYPE);
	    popSizeArray[1] = totalPopSize*PROP_STRAIN1*PROP_WILDTYPE;
	    popSizeArray[2] = totalPopSize*(1-PROP_STRAIN1)*(1 - PROP_WILDTYPE);
	    popSizeArray[3] = totalPopSize*(1-PROP_STRAIN1)*PROP_WILDTYPE;

	}

	if (NUM_STRAIN == 3){

	    if (PROP_STRAIN1==0 && PROP_STRAIN2==0){
	    		PROP_STRAIN1 = 1/((double)NUM_STRAIN);
	    		PROP_STRAIN2 = 1/((double)NUM_STRAIN);
	    }
	    popSizeArray[0] = totalPopSize*PROP_STRAIN1*(1 - PROP_WILDTYPE);
	    popSizeArray[1] = totalPopSize*PROP_STRAIN1*PROP_WILDTYPE;
	    popSizeArray[2] = totalPopSize*PROP_STRAIN2*(1 - PROP_WILDTYPE);
	    popSizeArray[3] = totalPopSize*PROP_STRAIN2*PROP_WILDTYPE;
	    popSizeArray[4] = totalPopSize*(1-PROP_STRAIN1-PROP_STRAIN2)*(1 - PROP_WILDTYPE);
	    popSizeArray[5] = totalPopSize*(1-PROP_STRAIN1-PROP_STRAIN2)*PROP_WILDTYPE;
	}

	if (NUM_STRAIN == 4){

	    if (PROP_STRAIN1==0 && PROP_STRAIN2==0 && PROP_STRAIN3){
	    		PROP_STRAIN1 = 1/((double)NUM_STRAIN);
	    		PROP_STRAIN2 = 1/((double)NUM_STRAIN);
	    		PROP_STRAIN3 = 1/((double)NUM_STRAIN);
	    }

	    popSizeArray[0] = totalPopSize*PROP_STRAIN1*(1 - PROP_WILDTYPE);
	    popSizeArray[1] = totalPopSize*PROP_STRAIN1*PROP_WILDTYPE;
	    popSizeArray[2] = totalPopSize*PROP_STRAIN2*(1 - PROP_WILDTYPE);
	    popSizeArray[3] = totalPopSize*PROP_STRAIN2*PROP_WILDTYPE;
	    popSizeArray[4] = totalPopSize*PROP_STRAIN3*(1 - PROP_WILDTYPE);
	    popSizeArray[5] = totalPopSize*PROP_STRAIN3*PROP_WILDTYPE;
	    popSizeArray[6] = totalPopSize*(1-PROP_STRAIN1-PROP_STRAIN2-PROP_STRAIN3)*(1 - PROP_WILDTYPE);
	    popSizeArray[7] = totalPopSize*(1-PROP_STRAIN1-PROP_STRAIN2-PROP_STRAIN3)*PROP_WILDTYPE;
	}
	return popSizeArray;
}

int *lineNumOfLysisDeath(string s){

	static int rowNumArray[6];
	int a;

	for(a = 0; a<6; a++){
		rowNumArray[a] = 0;
	}

    if (s.compare("Xn") ==0) {
    		rowNumArray[0] = 4;
    		rowNumArray[1] = 9;
    }
    else if (s.compare("Xr") ==0 || s.compare("Xnr") ==0) {
		rowNumArray[0] = 6;
		rowNumArray[1] = 13;
    }
    else if (s.compare("XaXn") ==0) {
		rowNumArray[0] = 12;
		rowNumArray[1] = 17;
    }
    else if (s.compare("XaXr") ==0 || s.compare("XaXnr") ==0) {
		rowNumArray[0] = 14;
		rowNumArray[1] = 21;
    }
    else if (s.compare("XnXr") ==0 || s.compare("XnXnr") ==0) {
		rowNumArray[0] = 4;
		rowNumArray[1] = 9;
		rowNumArray[2] = 16;
		rowNumArray[3] = 23;
    }

    else if (s.compare("XrXnr") ==0) {
		rowNumArray[0] = 6;
		rowNumArray[1] = 13;
		rowNumArray[2] = 20;
		rowNumArray[3] = 27;
    }
    else if (s.compare("XaXnXr") ==0 || s.compare("XaXnXnr") ==0) {
		rowNumArray[0] = 12;
		rowNumArray[1] = 17;
		rowNumArray[2] = 24;
		rowNumArray[3] = 31;
    }
    else if (s.compare("XaXrXnr") ==0) {
		rowNumArray[0] = 14;
		rowNumArray[1] = 21;
		rowNumArray[2] = 28;
		rowNumArray[3] = 35;
    }
    else if (s.compare("XnXrXnr") ==0) {
		rowNumArray[0] = 4;
		rowNumArray[1] = 9;
		rowNumArray[2] = 16;
		rowNumArray[3] = 23;
		rowNumArray[4] = 30;
		rowNumArray[5] = 37;
    }
    else if (s.compare("XaXnXrXnr") ==0) {
		rowNumArray[0] = 12;
		rowNumArray[1] = 17;
		rowNumArray[2] = 24;
		rowNumArray[3] = 31;
		rowNumArray[4] = 38;
		rowNumArray[5] = 45;
    }

	return rowNumArray;
}

int *lineNumOfGrowth2(string s) {
	static int rowNumArray[4];
		int a;

		for(a = 0; a<4; a++){
			rowNumArray[a] = 0;
		}

		if (s.compare("Xn") ==0) {
	    		rowNumArray[0] = 0; //growth2
	    		rowNumArray[1] = 5; //growth2
	    }
	    else if (s.compare("Xnr") ==0) {
	    		rowNumArray[0] = 0; //growth2
	    		rowNumArray[1] = 7; //growth2
	    }
	    else if (s.compare("XaXn") ==0) {
			rowNumArray[0] = 8; //growth2
			rowNumArray[1] = 13; //growth2
	    }
	    else if (s.compare("XaXnr") ==0) {
			rowNumArray[0] = 8; // growth2
			rowNumArray[1] = 15; // grwoth2
	    }
	    else if (s.compare("XnXr") ==0) {
			rowNumArray[0] = 0;
			rowNumArray[1] = 5;
	    }
	    else if (s.compare("XnXnr") ==0) {
			rowNumArray[0] = 0;
			rowNumArray[1] = 5;
			rowNumArray[2] = 10;
			rowNumArray[3] = 17;
	    }
	    else if (s.compare("XrXnr") ==0) {
			rowNumArray[0] = 14;
			rowNumArray[1] = 21;
	    }
	    else if (s.compare("XaXnXr") ==0) {
			rowNumArray[0] = 8;
			rowNumArray[1] = 13;
	    }
	    else if (s.compare("XaXnXnr") ==0) {
			rowNumArray[0] = 8;
			rowNumArray[1] = 13;
			rowNumArray[2] = 18;
			rowNumArray[3] = 25;
	    }
	    else if (s.compare("XaXrXnr") ==0) {
			rowNumArray[0] = 21;
			rowNumArray[1] = 28;
	    }
	    else if (s.compare("XnXrXnr") ==0) {
			rowNumArray[0] = 0;
			rowNumArray[1] = 5;
			rowNumArray[2] = 24;
			rowNumArray[3] = 31;
	    }
	    else if (s.compare("XaXnXrXnr") ==0) {
			rowNumArray[0] = 8;
			rowNumArray[1] = 13;
			rowNumArray[2] = 32;
			rowNumArray[3] = 39;
	    }

		return rowNumArray;

}


int *lineNumOfRecomb(string s){

	static int rowNumArray[8];
	int a;

	for(a = 0; a<8; a++){
		rowNumArray[a] = 0;
	}

    if (s.compare("Xr") ==0) {
		rowNumArray[0] = 4; //recombp
		rowNumArray[1] = 5; //recombn
		rowNumArray[2] = 11; //recombp
		rowNumArray[3] = 12; //recombn
    }
    else if (s.compare("Xnr") ==0) {
    		rowNumArray[0] = 4; //recombp
    		rowNumArray[1] = 5; //recombn
		rowNumArray[2] = 11; //recombp
		rowNumArray[3] = 12; //recombn
    }
    else if (s.compare("XaXr") ==0) {
		rowNumArray[0] = 12;
		rowNumArray[1] = 13;
		rowNumArray[2] = 19;
		rowNumArray[3] = 20;
    }
    else if (s.compare("XaXnr") ==0) {
		rowNumArray[0] = 12;
		rowNumArray[1] = 13;
		rowNumArray[2] = 19;
		rowNumArray[3] = 20;
    }
    else if (s.compare("XnXr") ==0) {
		rowNumArray[0] = 14;
		rowNumArray[1] = 15;
		rowNumArray[2] = 21;
		rowNumArray[3] = 22;
    }
    else if (s.compare("XnXnr") ==0) {
		rowNumArray[0] = 14;
		rowNumArray[1] = 15;
		rowNumArray[2] = 21;
		rowNumArray[3] = 22;
    }
    else if (s.compare("XrXnr") ==0) {
		rowNumArray[0] = 4;
		rowNumArray[1] = 5;
		rowNumArray[2] = 11;
		rowNumArray[3] = 12;
		rowNumArray[4] = 18;
		rowNumArray[5] = 19;
		rowNumArray[6] = 25;
		rowNumArray[7] = 26;
    }
    else if (s.compare("XaXnXr") ==0) {
		rowNumArray[0] = 22;
		rowNumArray[1] = 23;
		rowNumArray[2] = 29;
		rowNumArray[3] = 30;
    }
    else if (s.compare("XaXnXnr") ==0) {
		rowNumArray[0] = 22;
		rowNumArray[1] = 23;
		rowNumArray[2] = 29;
		rowNumArray[3] = 30;
    }
    else if (s.compare("XaXrXnr") ==0) {
		rowNumArray[0] = 12;
		rowNumArray[1] = 13;
		rowNumArray[2] = 19;
		rowNumArray[3] = 20;
		rowNumArray[4] = 25;
		rowNumArray[5] = 26;
		rowNumArray[6] = 32;
		rowNumArray[7] = 33;

    }
    else if (s.compare("XnXrXnr") ==0) {
		rowNumArray[0] = 14;
		rowNumArray[1] = 15;
		rowNumArray[2] = 21;
		rowNumArray[3] = 22;
		rowNumArray[4] = 28;
		rowNumArray[5] = 29;
		rowNumArray[6] = 35;
		rowNumArray[7] = 36;
    }
    else if (s.compare("XaXnXrXnr") ==0) {
		rowNumArray[0] = 22;
		rowNumArray[1] = 23;
		rowNumArray[2] = 29;
		rowNumArray[3] = 30;
		rowNumArray[4] = 36;
		rowNumArray[5] = 37;
		rowNumArray[6] = 43;
		rowNumArray[7] = 44;
    }
	return rowNumArray;
}

int weightedSampling(long *rnum1, int *reactionArray){ // this reaction arry should include all the reactions
    int x=0;
    int j=0;
    int temp = 0;
    for (j=0; j < TOTAL_NUM_REACTIONS; j++){
        if (*rnum1 > temp && *rnum1 <= temp + reactionArray [j]){
        	x =j;
        }
        temp = temp + reactionArray[j];
    }
    return x;
}


double calcTau(int *reactionArray, double *rnum2){
    double tau = 0; // for update the time;
    double sum =  0;
    int i=0;
    for (i =0; i < TOTAL_NUM_REACTIONS; i++){
        sum = sum + (double)reactionArray[i];
        //printf("reactions[%d] \t %d \n", i, reactionArray[i]);
    }
    //printf("sum of array is %lf\n", sum);
    //printf("log sample is %lf\n", log(sample));
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
int isNegativePopSize(int *popSizeArray){ // this pop size array should include all the populations
    int a;
    //int y=0;
    int positivePopSize = 0;
    for (a=0; a< TOTAL_NUM_POPULATION; a++){
    		if(popSizeArray[a] >= 0){
    			positivePopSize = positivePopSize +1;
    			//printf("popSize[%d]\t%d\n", a, popSizeArray[a]);
    		}
    }
    //printf("positivePopSize is %d\n", positivePopSize);
    if (positivePopSize==0){
        //y =0;
        //printf("positive pop is %d\n", positivePopSize);
        //printf("Y = %d\n", y);
    		return(0);
    }
    else {
        //y =1;
        //printf("positive pop is %d\n", positivePopSize);
        //printf("Y = %d\n", y);
    		return(1);
    }
}

// no reaction is negative
int isNegativeReaction(int *reactionArray){
    int c=0;
    //int y;
    int positiveReaction = 0;
    for (c=0; c < TOTAL_NUM_REACTIONS; c++){
        if (reactionArray[c]>=0) {
        		positiveReaction = positiveReaction +1;
        		//printf("reaction[%d]\t%d\n", c, reactionArray[c]);
        }
    }
    //printf("positive reaction is %d\n", positiveReaction);
    if (positiveReaction == TOTAL_NUM_REACTIONS) {
        //y =1;
        //printf("positive reaction is %d\n", positiveReaction);
        //printf("Y = %d\n", y);
        return (1);
    }
    else {
    	    //y =0;
        //printf("positive reaction is %d\n", positiveReaction);
        //printf("Y = %d\n", y);
        return (0);
    }
}


