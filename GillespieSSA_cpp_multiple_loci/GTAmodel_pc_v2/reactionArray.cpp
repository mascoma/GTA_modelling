//
//  reactionArray.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/23/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <iostream>
#include "GTAcpp_G.h"


int * calcReaction(long *geno_gene, long *totalPopSize, int **genotypePool, string popType, long popSizeXa){
    //static long reactions[(NUM_OF_PROCESSES*genotypeNum+NUM_OF_PROCESSES_f*numofgeneReact)];
	static int reactions[MAX_NUM_ARRAY];
	int a=0;
    int totalGeneSize =0;
    int temp = 0;
    int i=0;
    int delta_growth;
    int *locusIndex;
    //static int gtArray[NUM_LOCI];
    static int gtArray[MAX_NUM_LOCI];
    int j=0;
    int fitness=0;
    int delta_death;
    int negative_mutation;
    int *positiveMutIndex;
    int y=0;
    int takenUp;
    int r=0;
    int lysisNum;
    int negative_recomb;
    int * positiveRecIndex;
    int popIndex =0;
    int k=0;
    //static long genotypeSize[genotypeNum];
    static long genotypeSize[MAX_genotypeNum];
    if (popType.compare("Xa") == 0) { popIndex = 1;}
    else if (popType.compare("Xn")== 0) {popIndex = 2;
    } else if (popType.compare("Xr") == 0) {popIndex = 3;}
    
    switch(popIndex) {
        case 1:
           
            for(a=genotypeNum; a <(genotypeNum+genetypeNum); a++){
                temp = temp + geno_gene[a];
            }
            totalGeneSize = temp;
            for (i=0; i < genotypeNum; i++){
                delta_growth = growth(geno_gene[i],totalPopSize, 0, totalGeneSize);
                locusIndex = genotypeArrayCreator(i, genotypePool);
                for (j=0; j<NUM_LOCI; j++){
                    gtArray[j]=*(locusIndex+j);
                }
                fitness = genotypeFitness(gtArray);
                delta_death = death(fitness, geno_gene[i], 0);
                negative_mutation = mutationNegative((double)geno_gene[i], 0);
                //lysisNum = lysis((double)geno_gene[i]);
                //negative_recomb = recombNegative(i, geno_gene, genotypePool, LYSIS);
                reactions[(NUM_OF_PROCESSES*i)] = delta_growth;
                reactions[(NUM_OF_PROCESSES*i+1)] = delta_death;
                reactions[(NUM_OF_PROCESSES*i+3)] = negative_mutation;
                reactions[(NUM_OF_PROCESSES*i+4)] = 0;
                reactions[(NUM_OF_PROCESSES*i+6)] = 0;
            }
            

            for (k=0; k < (genotypeNum); k++){
                genotypeSize[k] = geno_gene[k];
            }
            positiveMutIndex = genotypeMutationPositive(genotypeSize, genotypePool, 0);
     
            for (y=0; y < genotypeNum; y++) {
                reactions[(NUM_OF_PROCESSES*y+2)] = *(positiveMutIndex+y);
                reactions[(NUM_OF_PROCESSES*y+5)] = 0;
            }
            
            for (r=genotypeNum; r < (genotypeNum+genetypeNum); r++){
                reactions[((r-genotypeNum)+genotypeNum*7)] = 0;
            }
            
            break;
 
        case 2:
            for(a=genotypeNum; a <(genotypeNum+genetypeNum); a++){
                temp = temp + geno_gene[a];
            }
            totalGeneSize = temp;
            for (i=0; i < genotypeNum; i++){
                delta_growth = growth(geno_gene[i],totalPopSize, LYSIS, totalGeneSize);
                locusIndex = genotypeArrayCreator(i, genotypePool);
                for (j=0; j<NUM_LOCI; j++){
                    gtArray[j]=*(locusIndex+j);
                }
                fitness = genotypeFitness(gtArray);
                delta_death = death(fitness, geno_gene[i], LYSIS);
                negative_mutation = mutationNegative((double)geno_gene[i], LYSIS);
                lysisNum = lysis((double)geno_gene[i]);
                //negative_recomb = recombNegative(i, geno_gene, genotypePool, LYSIS);
                reactions[(NUM_OF_PROCESSES*i)] = delta_growth;
                reactions[(NUM_OF_PROCESSES*i+1)] = delta_death;
                reactions[(NUM_OF_PROCESSES*i+3)] = negative_mutation;
                reactions[(NUM_OF_PROCESSES*i+4)] = lysisNum;
                reactions[(NUM_OF_PROCESSES*i+6)] = 0;
            }
            
            for (k=0; k < (genotypeNum); k++){
                genotypeSize[k] = geno_gene[k];
            }
            positiveMutIndex = genotypeMutationPositive(genotypeSize, genotypePool, LYSIS);
            //positiveRecIndex = genotypeRecombPositive(geno_gene, genotypePool, LYSIS);
            for (y=0; y < genotypeNum; y++) {
                reactions[(NUM_OF_PROCESSES*y+2)] = *(positiveMutIndex+y);
                reactions[(NUM_OF_PROCESSES*y+5)] = 0;
            }
            for (r=genotypeNum; r < (genotypeNum+genetypeNum); r++){
                takenUp = geneTakenUp((double)geno_gene[r], totalPopSize, popSizeXa);
                reactions[((r-genotypeNum)+genotypeNum*NUM_OF_PROCESSES)] = takenUp;
            }

            break;
            
        case 3:
            for(a=genotypeNum; a <(genotypeNum+genetypeNum); a++){
                temp = temp + geno_gene[a];
            }
            totalGeneSize = temp;
            for (i=0; i < genotypeNum; i++){
                delta_growth = growth(geno_gene[i],totalPopSize, LYSIS, totalGeneSize);
                locusIndex = genotypeArrayCreator(i, genotypePool);
                for (j=0; j<NUM_LOCI; j++){
                    gtArray[j]=*(locusIndex+j);
                }
                fitness = genotypeFitness(gtArray);
                delta_death = death(fitness, geno_gene[i], LYSIS);
                negative_mutation = mutationNegative((double)geno_gene[i], LYSIS);
                lysisNum = lysis((double)geno_gene[i]);
                negative_recomb = recombNegative(i, geno_gene, genotypePool, LYSIS);
                reactions[(NUM_OF_PROCESSES*i)] = delta_growth;
                reactions[(NUM_OF_PROCESSES*i+1)] = delta_death;
                reactions[(NUM_OF_PROCESSES*i+3)] = negative_mutation;
                reactions[(NUM_OF_PROCESSES*i+4)] = lysisNum;
                reactions[(NUM_OF_PROCESSES*i+6)] = negative_recomb;
            }
 
            for (k=0; k < (genotypeNum); k++){
                genotypeSize[k] = geno_gene[k];
            }
            positiveMutIndex = genotypeMutationPositive(genotypeSize, genotypePool, LYSIS);
            positiveRecIndex = genotypeRecombPositive(geno_gene, genotypePool, LYSIS);
            for (y=0; y < genotypeNum; y++) {
                reactions[(NUM_OF_PROCESSES*y+2)] = *(positiveMutIndex+y);
                reactions[(NUM_OF_PROCESSES*y+5)] = *(positiveRecIndex+y);
            }
            for (r=genotypeNum; r < (genotypeNum+genetypeNum); r++){
                takenUp = geneTakenUp((double)geno_gene[r], totalPopSize, popSizeXa);
                reactions[((r-genotypeNum)+genotypeNum*NUM_OF_PROCESSES)] = takenUp;
            }
            
            break;
            
        default:
            printf("argument of population type is missing!\n");
            break;
    }
 
    //int x=0;
    //for (x=0; x < numOfReactions; ++x){
    //printf("reaction[%d]: %li\n", x, reactions[x]);
    //}
    return reactions;
}


int *reactArray(long *geno_gene, long *totalPopSize, int **genotypePool, string s){
    int reactSysIndex = 0;
    int i = 0;
    int a = 0;
    int b = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    long popSizeXa = 0;
    reactSysIndex = reactionSystem(s);
    //static long reactArray[(numOfgeonReact+numofgeneReact)];
    static int reactArray[MAX_NUM_ARRAY];
    int *reactIndex1;
    int *reactIndex2;
    int *reactIndex3;
    long geno_gene_Xa[(genotypeNum+genetypeNum)];
    long geno_gene_Xn[(genotypeNum+genetypeNum)];
    long geno_gene_Xr[(genotypeNum+genetypeNum)];

    switch (reactSysIndex) {
        case 1:
            reactIndex1 = calcReaction(geno_gene, totalPopSize, genotypePool, "Xa", 0);
            for (i =0; i< (numOfgeonReact+numofgeneReact); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            break;
            
        case 2:
            reactIndex1 = calcReaction(geno_gene, totalPopSize, genotypePool, "Xn", 0);
            for (i =0; i< (numOfgeonReact+numofgeneReact); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            break;
            
        case 3:
            reactIndex1 = calcReaction(geno_gene, totalPopSize, genotypePool, "Xr", 0);
            for (i =0; i< (numOfgeonReact+numofgeneReact); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            break;
            
        case 4:
            for (a =0; a<genotypeNum; a++){
                geno_gene_Xa[a]= *(geno_gene+a);
                geno_gene_Xn[a] = *(geno_gene+genotypeNum+a);
                popSizeXa = popSizeXa + geno_gene[a];
            }
            
            for (b = 0; b <genetypeNum; b++){
                geno_gene_Xa[b+genotypeNum] = 0;
                geno_gene_Xn[b+genotypeNum] = *(geno_gene+genotypeNum*2+b);
            }
            
            reactIndex1 = calcReaction(geno_gene_Xa, totalPopSize, genotypePool, "Xa", 0);
            for (i =0; i< (NUM_OF_PROCESSES*genotypeNum); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            
            reactIndex2 = calcReaction(geno_gene_Xn, totalPopSize, genotypePool, "Xn", popSizeXa);
            for (j = 0; j<(NUM_OF_PROCESSES*genotypeNum+genetypeNum*NUM_OF_PROCESSES_f); j++){
                reactArray[(NUM_OF_PROCESSES*genotypeNum+j)] = *(reactIndex2 +j);
            }
            
            break;
        
        case 5:
            for (a =0; a<genotypeNum; a++){
                geno_gene_Xa[a]= *(geno_gene+a);
                geno_gene_Xr[a] = *(geno_gene+genotypeNum+a);
                popSizeXa = popSizeXa + geno_gene[a];
            }
            
            for (b = 0; b <genetypeNum; b++){
                geno_gene_Xa[b+genotypeNum] = 0;
                geno_gene_Xr[b+genotypeNum] = *(geno_gene+genotypeNum*2+b);
            }
            
            reactIndex1 = calcReaction(geno_gene_Xa, totalPopSize, genotypePool, "Xa", 0);
            for (i =0; i< (NUM_OF_PROCESSES*genotypeNum); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            
            reactIndex2 = calcReaction(geno_gene_Xr, totalPopSize, genotypePool, "Xr", popSizeXa);
            for (j = 0; j<(NUM_OF_PROCESSES*genotypeNum+genetypeNum*NUM_OF_PROCESSES_f); j++){
                reactArray[(NUM_OF_PROCESSES*genotypeNum+j)] = *(reactIndex2 +j);
            }

            break;
            
        case 6:
            for (a =0; a<genotypeNum; a++){
                geno_gene_Xn[a]= *(geno_gene+a);
                geno_gene_Xr[a] = *(geno_gene+genotypeNum+a);
            }
            
            for (b = 0; b <genetypeNum; b++){
                geno_gene_Xn[b+genotypeNum] = *(geno_gene+genotypeNum*2+b);
                geno_gene_Xr[b+genotypeNum] = *(geno_gene+genotypeNum*2+b);
            }
            
            reactIndex1 = calcReaction(geno_gene_Xn, totalPopSize, genotypePool, "Xn", 0);
            for (i =0; i< (NUM_OF_PROCESSES*genotypeNum); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            
            reactIndex2 = calcReaction(geno_gene_Xr, totalPopSize, genotypePool, "Xr", 0);
            for (j =0; j<(NUM_OF_PROCESSES*genotypeNum); j++){
                reactArray[(NUM_OF_PROCESSES*genotypeNum+j)] = *(reactIndex2 +j);
            }
            for (k = 0; k< (NUM_OF_PROCESSES_f*genetypeNum); k++){
            
                reactArray[(NUM_OF_PROCESSES*genotypeNum*2+k)] = *(reactIndex1+NUM_OF_PROCESSES*genotypeNum+k) + *(reactIndex2+NUM_OF_PROCESSES*genotypeNum+k);
            }
            break;
        
        case 7:
            for (a =0; a<genotypeNum; a++){
                geno_gene_Xa[a]= *(geno_gene+a);
                geno_gene_Xn[a] = *(geno_gene+genotypeNum+a);
                geno_gene_Xr[a] = *(geno_gene+genotypeNum*2+a);
                popSizeXa = popSizeXa + geno_gene[a];
            }
            for (b = 0; b <genetypeNum; b++){
                geno_gene_Xa[b+genotypeNum] = 0;
                geno_gene_Xn[b+genotypeNum] = *(geno_gene+genotypeNum*3+b);
                geno_gene_Xr[b+genotypeNum] = *(geno_gene+genotypeNum*3+b);
            }
            reactIndex1 = calcReaction(geno_gene_Xa, totalPopSize, genotypePool, "Xa", 0);
            for (i =0; i< (NUM_OF_PROCESSES*genotypeNum); i++){
                reactArray[i] = *(reactIndex1 +i);
            }
            
            reactIndex2 = calcReaction(geno_gene_Xn, totalPopSize, genotypePool, "Xn", popSizeXa);
            for (j = 0; j<(NUM_OF_PROCESSES*genotypeNum); j++){
                reactArray[(NUM_OF_PROCESSES*genotypeNum+j)] = *(reactIndex2 +j);
            }
            
            reactIndex3 = calcReaction(geno_gene_Xr, totalPopSize, genotypePool, "Xr", popSizeXa);
            for (k = 0; k<(NUM_OF_PROCESSES*genotypeNum); k++){
                reactArray[(NUM_OF_PROCESSES*genotypeNum*2+k)] = *(reactIndex3 +k);
            }
            for (l = 0; l<  (NUM_OF_PROCESSES_f*genetypeNum); l++){
                
                reactArray[(NUM_OF_PROCESSES*genotypeNum*3+l)] = *(reactIndex1+NUM_OF_PROCESSES*genotypeNum+l) + *(reactIndex2+NUM_OF_PROCESSES*genotypeNum+l) + *(reactIndex3+NUM_OF_PROCESSES*genotypeNum+l);
            }

            
            break;
            
        default:
            printf("argument of population type is missing!\n");
            break;
    }
    return reactArray;
}
