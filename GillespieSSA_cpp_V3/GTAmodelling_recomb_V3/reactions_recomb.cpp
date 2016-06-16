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

#include "GTAcpp_recomb.h"

int growth(int popSize, int *totalPopSize, double lysis){ //logistic growth
    int delta_growth = 0;
    double carryCap = (double)CARRYING_CAPACITY; // has to be float type to get divide result not equal to 0
    double popRatio = 0;
    popRatio = ((double)*totalPopSize)/carryCap;
    delta_growth = popSize*(1-popRatio)*GROWTH_RATE*(1-lysis);
    //printf("delta growth is %li\n", delta_growth);
    return (delta_growth);
}

int growth2(int popSize, int *totalPopSize, long *numOfGTA, double lysis){ // growth because of nutro supply
    int delta_growth = 0;
    double carryCap = (double)CARRYING_CAPACITY; // has to be float type to get divide result not equal to 0
    double popRatio = 0;
    popRatio = ((double)*totalPopSize)/carryCap;
    delta_growth = popSize*(1-popRatio)*(GROWTH_RATE + NUTRITIONAL_VALUE*UPTAKE_RATE**numOfGTA)*(1-lysis);
    //printf("delta growth is %li\n", delta_growth);
    return (delta_growth);
}

int lysisDeath(int popSize){ // lysis decrease the number of bacteria cells
 
    double lysisNum = 0;
    int lyNum = 0;
    lysisNum = LYSIS*popSize;
    lyNum = (int)lysisNum;
    return(lyNum);
}

int death(double fitness, double lysis, int popSize){
	// fitness related death
	int delta_death = 0;
	delta_death = popSize*SELECTION*(1-fitness)*(1-lysis);
	return(delta_death);
}

int mut(double lysis, int popSize){ // mutation transfer genoytype from 0 to 1
	int delta_mutp = 0;
	delta_mutp = popSize*MUTATION*(1/NUM_LOCI)*(1-lysis);
	return(delta_mutp);
}


int recomb(long *numOfGTA, int popSize){ // recomb from 0 to 1
	int delta_recomb = 0;
	delta_recomb = popSize**numOfGTA*UPTAKE_RATE*INCORPORATION_RATE*(1-LYSIS);
	if (delta_recomb>=popSize || delta_recomb>=*numOfGTA){
		return(delta_recomb);
	}
	else if (popSize>*numOfGTA){
		delta_recomb = *numOfGTA;
		return(delta_recomb);
	}
	else if (*numOfGTA>=popSize){
		delta_recomb = popSize;
		return(delta_recomb);
	}
}


long GTAIncrease(int popSize){ // increasing of certain type of GTA (contain fragment 0 or 1)
	int delta_GTA = 0;
	delta_GTA = popSize*LYSIS*BURST_SIZE*(1/NUM_LOCI);
	return(delta_GTA);
}


long GTADecrease(long *numOfGTA){ //
	long delta_GTA = 0;
	delta_GTA = *numOfGTA*UPTAKE_RATE;
	//printf("number of GTA decreasing is %li \n", delta_GTA);
	return(delta_GTA);
}

