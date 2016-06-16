//
//  reactions_selfish.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <cstdio>
#include <cmath>
#include <algorithm>

#include "GTAcpp_selfish.h"

int growth_popSize1(long popSize1, long carryCapacity, long *totalPopSize, double growth_rate , double lysis){ //logistic growth of GTA+ pop
    int delta_growth = 0;
    double carryCap = (double)carryCapacity; // has to be float type to get divide result not equal to 0
    double popRatio = 0;
    popRatio = ((double)*totalPopSize)/carryCap;
    delta_growth = (1-popRatio)*growth_rate*(1-lysis)*popSize1;
    //printf("delta growth is %li\n", delta_growth);
    return (delta_growth);
}

int growth_popSize2(long popSize2, long carryCapacity, long *totalPopSize, double growth_rate){ //logistic growth of GTA- pop
    int delta_growth = 0;
    double carryCap = (double)carryCapacity; // has to be float type to get divide result not equal to 0
    double popRatio = 0;
    popRatio = ((double)*totalPopSize)/carryCap;
    delta_growth = (1-popRatio)*growth_rate*popSize2;
    //printf("delta growth is %li\n", delta_growth);
    return (delta_growth);
}


int lysis(double popSize1, double lysis){ // lysis decrease the number of bacteria cells
 
    double lysisNum = 0;
    int lyNum = 0;
    lysisNum = lysis*popSize1;
    lyNum = (int)lysisNum;
    return(lyNum);
}


int transfer(double lysis, int burst, double infect, double recomb, double large_head_freq, int popSize1, int popSize2){
	int transferNum = 0;
	transferNum = lysis*burst*infect*recomb*large_head_freq*popSize1*popSize2;
	return(transferNum);
}




