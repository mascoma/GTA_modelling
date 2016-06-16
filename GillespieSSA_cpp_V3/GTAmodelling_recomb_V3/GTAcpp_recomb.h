//
//  GTAcppmodelling_selfish
//
//  Created by Xin Chen on 5/31/16.
//  Copyright (c) 2016 Xin Chen. All rights reserved.
//

#ifndef __GTAcppmodelling__GTAcpp_selfish__
#define __GTAcppmodelling__GTAcpp_selfish__

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;


#define TIME_INIT 0
#define MAX_NUM_ARRAY 50

extern int CARRYING_CAPACITY;
extern int INIT_TOTALPOPSIZE;
extern double PROP_WILDTYPE; // genotype 0 is wildtype
extern double GROWTH_RATE;
extern double MUTATION;
extern double NUTRITIONAL_VALUE;
extern double UPTAKE_RATE;
extern double INCORPORATION_RATE;
extern double LYSIS;
extern double SELECTION;
extern double FITNESS;
extern int BURST_SIZE;
extern double NUM_LOCI;
extern bool CONSOLE;
extern int TIME_MAX;
extern string OUTPUT;
// parameters built the system
extern int NUM_STRAIN;
extern double PROP_STRAIN1;
extern double PROP_STRAIN2;
extern double PROP_STRAIN3;
extern string STRAIN_NAME;
extern int TOTAL_NUM_REACTIONS;
extern int TOTAL_NUM_POPULATION;
extern int NUM_LYSIS;
extern int NUM_G2;
extern int NUM_RECOMB;
extern int SAMPLING_FREQ;  // frequency of storing results in the output file


int growth(int, int *totalPopSize, double);
int growth2(int, int *totalPopSize, long *numOfGTA, double);
int lysisDeath(int);
int death(double, double, int);
int mut(double, int);
int recomb(long *numOfGTA, int);
long GTAIncrease(int);
long GTADecrease(long *numOfGTA);

int weightedSampling(long *rnum1, int *array);
double calcTau(int *reactionArray, double *rnum2);
int isNegativePopSize(int *popSizeArray);
int isNegativeReaction(int *reactionArray);
int *reactArray(int *popSize, int *totalPopSize, double, double, long *numOfGTA1, long *numOfGTA0, long *totalGTA, string s);
int reactionSystem(string s);
int *lineNumOfLysisDeath(string s);
int *lineNumOfGrowth2(string s);
int *lineNumOfRecomb(string s);
int *initPopSize(int);
long  *GTAdynamic(int *lineNum, int *popSizeArray, long *GTA1, long *GTA0);


struct statusMatrix {
    int reactStat[MAX_NUM_ARRAY][MAX_NUM_ARRAY];
};

statusMatrix getArray(string index);




#endif /* defined(__GTAcppmodelling__GTAcpp_Xr__) */
