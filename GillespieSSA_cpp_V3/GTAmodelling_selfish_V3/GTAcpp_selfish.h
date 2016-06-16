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





extern long CARRYING_CAPACITY;
extern long INIT_TOTALPOPSIZE;
//extern int NUM_LOCI;
extern double PROP_WILDTYPE; // non GTA group
extern double GROWTH_RATE;
extern double UPTAKE_RATE;
extern double INCORPORATION_RATE;
extern double LYSIS;
extern int BURST_SIZE;
extern double LARGE_HEAD_FREQ;
extern bool CONSOLE;
extern int TIME_MAX;
extern string OUTPUT;


int growth_popSize1(long, long, long *totalPopSize, double, double);
int growth_popSize2(long, long, long *totalPopSize, double);
int lysis(double, double);
int transfer(double, int, double, double, double, int, int);

int weightedSampling(long *rnum1, int *array);
double calcTau(int *reactionArray, double *rnum2);
int isNegativePopSize(long *popSizeArray);
int isNegativeReaction(int *reactionArray);








#endif /* defined(__GTAcppmodelling__GTAcpp_Xr__) */
