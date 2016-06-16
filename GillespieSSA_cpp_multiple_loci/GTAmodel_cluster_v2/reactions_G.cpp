//
//  reactions_Xr.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "GTAcpp_G.h"

int growth(long popSize, long *totalPopSize, double lysis = LYSIS, int totalGeneSize=0){ //if Xa, lysis =0
    int delta_growth = 0;
    double carryCapacity = (double)CARRYING_CAPACITY; // has to be float type to get divide result not equal to 0
    double popRatio = 0;
    popRatio = ((double)*totalPopSize)/carryCapacity;
    delta_growth = (1-popRatio)*(GROWTH_RATE+NUTRITIONAL_VALUE*UPTAKE_RATE*totalGeneSize)*(1-lysis)*popSize;
    //printf("delta growth is %li\n", delta_growth);
    return (delta_growth);
}



int * genotypeArrayCreator(int genotypeIndex, int **genotypePool){
 
    //static int genotypeArray[NUM_LOCI];
    static int genotypeArray[MAX_NUM_LOCI];
    int a=0;
    for (a=0; a < NUM_LOCI; a++) {
        genotypeArray[a] = genotypePool[genotypeIndex][a];
    }
    return genotypeArray;
}


int genotypeFitness(int genotype[NUM_LOCI]){
    int i;
    int sum = 0;
    for (i=0; i < NUM_LOCI; i++) {
        sum = sum + genotype[i];
    }
    return (sum);
}

int death(int genotypeFitness, long popSize, double lysis =LYSIS){ // if Xa lysis =0
    int d_death = 0;
    d_death=(INTRI_DEATH_RATE + SELECTION*(NUM_LOCI-genotypeFitness))*(1-LYSIS)*popSize;
    return (d_death);
}


int mutationNegative(double popSize, double lysis = LYSIS){ // if Xa lysis =0
    double mutation_negative = 0;
    int mut_return = 0;
    mutation_negative = NUM_LOCI*MUTATION*(1-LYSIS)*popSize;
    mut_return = (int)mutation_negative;
    return(mut_return);
}


int genotypeDistance(int targetgenotype[NUM_LOCI], int genotype2[NUM_LOCI]){
    int genotypeDist = 0;
    int temp =0;
    int i;
    for (i=0; i < NUM_LOCI; ++i){
        temp = targetgenotype[i] - genotype2[i];
        genotypeDist = genotypeDist + abs(temp);
    }
    //printf("genotypeDistance is %d\n", genotypeDist);
    if (genotypeDist == 1){
        return (1);
    }
    else {
        return (0);
    }
}

// this function is to first determine if the distance between genotype2 and our target genotype equal to 1
// if so, then genotype2 could potentially change to our target genotype through a single mutation


double mutationPositive(int targetgenotype[NUM_LOCI], int comparedgenotype[NUM_LOCI],double popSize, double lysis =LYSIS){ // if Xa lysis =0
    int genotype2=0;
    double delta_mutatioin_genotyp2 = 0;
    genotype2 = genotypeDistance(targetgenotype, comparedgenotype);
    if (genotype2){
        delta_mutatioin_genotyp2 = popSize * MUTATION*(1-lysis);
        //printf("delta_mutatioin: %lf\n", delta_mutatioin_genotyp2);
    }
    else {
        delta_mutatioin_genotyp2 = 0;
        //printf("delta_mutatioin: %lf\n", delta_mutatioin_genotyp2);
    }
    return(delta_mutatioin_genotyp2);
}


int * genotypeMutationPositive(long *popSizeArray, int **genotypePool, double lysis =LYSIS){ // if Xa lysis =0
    //static int targetGenotype[NUM_LOCI];
    //static int comparedGenotype[NUM_LOCI];
    static int targetGenotype[MAX_NUM_LOCI];
    static int comparedGenotype[MAX_NUM_LOCI];
    int a=0;
    int b=0;
    int c=0;
    int d=0;
    double genotypePopSize= 0;
    //static long genotypeMutationPositiveNumber[genotypeNum];
    static int genotypeMutationPositiveNumber[MAX_genotypeNum];
    double mutationIncrease = 0;
    for (a=0; a < genotypeNum; a++){
        int *gtArray1; // a pointer for targetgenotype array
        double mutationPositiveNumber = 0;
        gtArray1 = genotypeArrayCreator(a, genotypePool);
        // need to use the every generation of popSize array
        for (b=0; b < NUM_LOCI; b++){
            targetGenotype[b] = *(gtArray1+b);
            //printf("tragetGT[%d]: %d\n", b, targetGenotype[b]);
        }
        for (c=0; c < genotypeNum; c++){
            int *gtArray2; // a pointer for comparedgenotype array
            gtArray2 = genotypeArrayCreator(c, genotypePool);
            for (d =0; d < NUM_LOCI; d++){
                comparedGenotype[d] = *(gtArray2+d);
                //printf("comparedGT[%d]: %d\n", d, comparedGenotype[d]);
            }
            //targetGenotype comparing comparedGenotype
            //genotypeDistance(targetGenotype[NUM_LOCI], comparedGenotype[NUM_LOCI]);
            genotypePopSize = (double)popSizeArray[c];
            //printf("compared genotype [%d] popsize is %lf\n",c, (double)popSizeArray[c]);
            mutationIncrease = mutationPositive(targetGenotype, comparedGenotype, genotypePopSize, lysis);
            //printf("%lf\n", *mutationIncrease);
            mutationPositiveNumber = mutationPositiveNumber + mutationIncrease;
        }
        //printf("%mutation of genotype[%d] %lf\n",a, *mutationPositiveNumber);
        //printf("%mutation of genotype[%d] %li\n",a, (long)*mutationPositiveNumber);
        genotypeMutationPositiveNumber[a] = (int)mutationPositiveNumber;
        //printf("genotypeMutationPositiveNumber[%d]: %li\n",a, genotypeMutationPositiveNumber[a]);
    }
    return genotypeMutationPositiveNumber;
}

int geneTakenUp (double numOfGene, long *totalPopSize, long popSizeXa) { // influence the number of different type of GTA particles
    double takenUpNum = 0;
    int takeupNum = 0;
    double temp;
    temp = (double)*totalPopSize - (double)popSizeXa;
    takenUpNum = numOfGene*temp*UPTAKE_RATE;
    if (takenUpNum <= temp && takenUpNum <= numOfGene){
    		takeupNum = (int)takenUpNum;
    }
    else if (takenUpNum > temp || takenUpNum > numOfGene){
        takeupNum = std::min((int)numOfGene, (int)temp);
    }
    	return(takeupNum);
}

int lysis(double popSize){ // lysis influence number of bacteria cells
    //as well as number of GTA particles containing different gene loci in environment
    double lysisNum = 0;
    int lyNum = 0;
    lysisNum = LYSIS*popSize;
    lyNum = (int)lysisNum;
    return(lyNum);
}

double genetypeSize1(int targetgenotype[NUM_LOCI], int comparedgenotype[NUM_LOCI], long *popSizeArray){
    // to determine if genotype distance is 1, which locus is different
    // then, based on that to selection corresponding gene fragment to calculate recombination
    static int genetypeArray[MAX_genotypeNum];
    int j=0;
    double genetypeSize = 0;
    for (j=genotypeNum; j < (genotypeNum+genetypeNum); j++){
        genetypeArray[(j-genotypeNum)] = popSizeArray[j];
        //printf("genetypesize [%d], %li\n",(j-genotypeNum), genetypeArray[(j-genotypeNum)]);
    }
    int i=0;
    for (i = 0; i < NUM_LOCI; i++){
        if (targetgenotype[i] != comparedgenotype[i]){
            if (comparedgenotype[i] == 0){
                genetypeSize = (double)genetypeArray[(2*i)];
            }
            else if (comparedgenotype[i] ==1){
                genetypeSize = (double)genetypeArray[(2*i+1)];
            }
            else {printf("something wrong\n");}
        }
    }
    //printf("genetype size of comparedgenotype is %d\n", genetypeSize);
    return (genetypeSize);
}

double genotypeSize(int comparedgenotype[NUM_LOCI], long *popSizeArray){
    static long genotypeArray[MAX_genotypeNum];
    int k=0;
    double genotypeSize = 0 ;
    for (k=0; k < genotypeNum; k++){
        genotypeArray[k] = popSizeArray[k];
    }
    int i=0;
    int temp=0;
    int genotypeIndex = 0;
    for(i=0; i < NUM_LOCI; ++i){
        temp = comparedgenotype[i]*pow(2, i);
        genotypeIndex = genotypeIndex + temp;
    }
    genotypeSize = (double)genotypeArray[genotypeIndex];
    //printf("genotype size of comparedgenotype is %lf\n", *genotypeSize);
    return (genotypeSize);
}

double recombPositive(int targetgenotype[NUM_LOCI], int comparedgenotype[NUM_LOCI],long *popSizeArray, double lysis=LYSIS){
    int genotype2 =0;
    double temp;
    double delta_recomb_genotyp2 = 0;
    genotype2 = genotypeDistance(targetgenotype, comparedgenotype);
    //printf("dist is %d\n", genotype2);
    if (genotype2){
        double genotypePopSize=0;
        double genetypePopSize=0;
        genotypePopSize = genotypeSize(comparedgenotype, popSizeArray);
        //printf("genotypeSize is %lf\n", genotypePopSize);
        genetypePopSize = genetypeSize1(targetgenotype, comparedgenotype, popSizeArray);
        //printf("genetypeSize is %lf\n", genetypePopSize);
         temp = genotypePopSize * genetypePopSize * UPTAKE_RATE * INCORPORATION_RATE * (1-lysis);
         if (temp <= genotypePopSize && temp <= genetypePopSize){
        	 	 delta_recomb_genotyp2 = temp;
         }
         else if (temp > genetypePopSize || temp > genotypePopSize){
        	 	 delta_recomb_genotyp2 = std::min(genetypePopSize, genotypePopSize);
         }

        //printf("delta_recomb: %lf\n", delta_recomb_genotyp2);
    }
    else{
        delta_recomb_genotyp2 = 0;
        //printf("delta_recomb: %lf\n", *delta_recomb_genotyp2);
    }
    return(delta_recomb_genotyp2);
}

int * genotypeRecombPositive (long *popSizeArray, int **genotypePool, double lysis =LYSIS){
    int *gtArray1; // a pointer for targetgenotype array
    int *gtArray2; // a pointer for comparedgenotype array
    //static int targetGenotype[NUM_LOCI];
    //static int comparedGenotype[NUM_LOCI];
    static int targetGenotype[MAX_NUM_LOCI];
    static int comparedGenotype[MAX_NUM_LOCI];
    int a=0;
    int b=0;
    int c=0;
    int d=0;
    double recombIncrease = 0;
    //int gt = numberOfGenotype(NUM_LOCI); // number of genotype
    //static long genotypeRecombPositiveNumber[genotypeNum];
    static int genotypeRecombPositiveNumber[MAX_genotypeNum];
    for (a=0; a < genotypeNum; a++){
        double recombPositiveNumber = 0;
        gtArray1 = genotypeArrayCreator(a, genotypePool);
        // need to use the every generation of popSize array
        for (b=0; b < NUM_LOCI; b++){
            targetGenotype[b] = *(gtArray1+b);
            //printf("tragetGT[%d]: %d\n", b, targetGenotype[b]);
        }
        for (c=0; c < genotypeNum; c++){
            gtArray2 = genotypeArrayCreator(c, genotypePool);
            for (d =0; d < NUM_LOCI; d++){
                comparedGenotype[d] = *(gtArray2+d);
                //printf("comparedGT[%d]: %d\n", d, comparedGenotype[d]);
            }
            //targetGenotype comparing comparedGenotype
            //genotypeDistance(targetGenotype[NUM_LOCI], comparedGenotype[NUM_LOCI]);
            //printf("compared genotype [%d] popsize is %lf\n",c, (double)genotypePopSizeArray[c]);
            recombIncrease = recombPositive(targetGenotype, comparedGenotype, popSizeArray, lysis);
            //printf("%lf\n", *recombIncrease);
            recombPositiveNumber = recombPositiveNumber + recombIncrease;
        }
        //printf("recomb of genotype[%d] %lf\n",a, recombPositiveNumber);
        //printf("recomb of genotype[%d] %d\n",a, (int)recombPositiveNumber);
        genotypeRecombPositiveNumber[a] = (int)recombPositiveNumber;
        //printf("genotypeRecombPositiveNumber[%d]: %li\n",a, genotypeRecombPositiveNumber[a]);
    }
    return genotypeRecombPositiveNumber;
}

double genetypeSize2(int genotypeIndex, long *popSizeArray, int **genotypePool){
    // to calculate the genetype size given a genotype,
    //for example, genotype 11 will have genotypeSize = genetype _1 + genetype 1_
    int * gtArray;// a pointer for a genotype array
    static int genotypeArray[MAX_NUM_LOCI];
    int i=0;
    int j=0;
    int k=0;
    gtArray = genotypeArrayCreator(genotypeIndex, genotypePool);
    for (i=0; i < NUM_LOCI; ++i){
        genotypeArray[i] = *(gtArray+i);
    }
    static int genetypeArray[MAX_genotypeNum];
    int temp = 0;
    double genetypeSize = 0;
    for (k=genotypeNum; k < (genotypeNum+genetypeNum); k++){
        genetypeArray[(k-genotypeNum)] = popSizeArray[k];
    }
    for (j = 0; j < NUM_LOCI; j++){
        if (genotypeArray[j] == 0){
            temp = genetypeArray[(2*j)];
        }
        else if (genotypeArray[j] ==1){
            temp = genetypeArray[(2*j+1)];
        }
        else {printf("something wrong\n");}
        genetypeSize = genetypeSize + temp;
    }
    return(genetypeSize);
}

int recombNegative(int genotypeIndex, long *popSizeArray, int **genotypePool, double lysis =LYSIS){
    int j=0;
    double temp = 0;
    double totalGenetypeSize = 0;
    double genetypeSize_gt; // genetype size for one specific genotype
    double negativeRecombNum = 0;
    int negRecombNum = 0;
    double temp1 =0;

    for (j=genotypeNum; j <(genotypeNum+genetypeNum); j++){
        temp = (double)popSizeArray[j];
        totalGenetypeSize = temp + totalGenetypeSize;
    }
    //printf("total genetype size %lf\n", *totalGenetypeSize);
    genetypeSize_gt = genetypeSize2(genotypeIndex, popSizeArray, genotypePool);
    //printf("genetypeSize %lf\n", *genetypeSize_gt);
    temp1 = totalGenetypeSize - genetypeSize_gt;
    if (negativeRecombNum <= popSizeArray[genotypeIndex] && negativeRecombNum <= temp1) {
    		negativeRecombNum = (double)popSizeArray[genotypeIndex]*temp1* UPTAKE_RATE * INCORPORATION_RATE *(1-lysis);
    }
    else if (negativeRecombNum > popSizeArray[genotypeIndex] || negativeRecombNum > temp1  ){
    		temp = popSizeArray[genotypeIndex];
    		negativeRecombNum = std::min(temp1, temp);
    }
    //printf("popSizeArray[%d]: %li\n", genotypeIndex, popSizeArray[genotypeIndex]);
    //printf("negativeRecombNum %lf\n", *negativeRecombNum);
    negRecombNum = (int)negativeRecombNum;
    return (negRecombNum);
}
