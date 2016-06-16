//
//  GTAcpp_Xr.h
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#ifndef __GTAcppmodelling__GTAcpp_Xr__
#define __GTAcppmodelling__GTAcpp_Xr__

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

#define MAX_NUM_ARRAY 500
#define MAX_NUM_LOCI 6
#define MAX_genotypeNum 100
#define TIME_INIT 0
#define GENERATIONFREQ 1000000 // frequency of storing results in the output file
#define NUM_OF_PROCESSES 7
//#define NUM_OF_PROCESSES_Xn 5
//#define NUM_OF_PROCESSES_Xr 7
#define NUM_OF_PROCESSES_f  1

/*#define CARRYING_CAPACITY 1000000000
#define INIT_TOTALPOPSIZE 900000000
#define NUM_LOCI 2
#define INTRI_DEATH_RATE 0.001
#define MUTATION 0.00000001
#define SELECTION 0.02
#define GROWTH_RATE 0.1
#define UPTAKE_RATE 0.0000000001
#define COST_COMPETENCE 0.001
#define NUTRITIONAL_VALUE 0.01
#define INCORPORATION_RATE 0.5
#define LYSIS 0.0005
//#define NUM_OF_STEPS 500
//#define CONSOLEINTERNAL 0
//#define CENSUSINTERNAL 0
//#define VERBOSE 1
//#define IGNORENEGATIVESTATE 1
#define TAU 0.00001  // STEP LENGTH
#define TIME_MAX 20000
//char method = 'D';
//int maxWallTime = 1.0/0.0; // positive infinity
#define genotypeNum 4 // number of genotypes = 2^number of loci
#define genetypeNum 4 // number of genetypes = 2*number of loci  genetype order in _0, 0_, _1, 1_
#define NUM_OF_PROCESSES 7
//#define NUM_OF_PROCESSES_Xn 5
//#define NUM_OF_PROCESSES_Xr 7
#define NUM_OF_PROCESSES_f  1

#define NUM_OF_POPTYPE 3
#define numOfgeonReact 84         // NUM_OF_PROCESSES * genotypeNum * NUM_OF_POPTYPE
#define numofgeneReact 4          // NUM_OF_PROCESSES_f * genetypeNUM


#define OUTPUT "simulation_XaXnXr.txt"
#define REACTIONSYSTEM "XaXnXr"
*/

extern long CARRYING_CAPACITY;
extern long INIT_TOTALPOPSIZE;
extern int NUM_LOCI;
extern double PROP_WILDTYPE;
extern double INTRI_DEATH_RATE;
extern double MUTATION;
extern double SELECTION;
extern double GROWTH_RATE;
extern double UPTAKE_RATE;
extern double NUTRITIONAL_VALUE;
extern double INCORPORATION_RATE;
extern double LYSIS;
extern bool CONSOLE;
extern int TIME_MAX;
extern int NUM_OF_POPTYPE;
extern string OUTPUT;
extern string REACTIONSYSTEM;

extern int genotypeNum;
extern int genetypeNum;
extern int numOfgeonReact;
extern int numofgeneReact;

int growth(long,long *totalPopSize, double, int);
int genotypeFitness(int *genotype);
int death(int, long, double);
int mutationNegative(double, double);
int genotypeDistance(int targetgenotype[NUM_LOCI], int genotype2[NUM_LOCI]);
double mutationPositive(int targetgenotype[NUM_LOCI], int comparedgenotype[NUM_LOCI],double popSize, double);
int * genotypeArrayCreator(int genotypeIndex, int **genotypePool);
int * genotypeMutationPositive (long *popSizeArray, int **genotypePool, double);
int reactionSystem(string s);
int matrixIndex(int rowNum, int colNum);
int * calcReaction(long *popSize, long *totalPopSize, int **genotypePool, string s);
int weightedSampling(long *rnum1,int *array);
double calcTau(int *reactionArray, double *rnum2);
int isNegativePopSize(long *popSizeArray);
int isNegativeReaction(int *reactionArray);
int geneTakenUp (double, long *totalPopSize, long popSizeXa);
int lysis(double);
double genetypeSize1(int targetgenotype[NUM_LOCI], int comparedgenotype[NUM_LOCI], long *popSizeArray);
double genotypeSize(int comparedgenotype[NUM_LOCI], long *popSizeArray);
double recombPositive(int targetgenotype[NUM_LOCI], int comparedgenotype[NUM_LOCI],long *popSizeArray, double);
int * genotypeRecombPositive (long *popSizeArray, int **genotypePool, double);
double genetypeSize2(int genotypeIndex, long *popSizeArray, int **genotypePool);
int recombNegative(int genotypeIndex, long *popSizeArray, int **genotypePool, double);
int *reactArray(long *geno_gene, long *totalPopSize, int **genetypePool, string s);
string headerName(string);
//genetype status during the lysis: so far just create it manually
// for Num of loci = 2,
//genotype 1 (00) increase genetype 1, 3 (_0, 0_)
//genotype 2 (01) increase genetype 2, 3 (_1, 0_)
//genotype 3 (10) increase genetype 1, 4 (_0, 1_)
//genotype 4 (11) increase genetype 2, 4 (_1, 1_)

struct statusMatrix {
    int reactStat[MAX_NUM_ARRAY][MAX_NUM_ARRAY];
};

statusMatrix getArray();




#endif /* defined(__GTAcppmodelling__GTAcpp_Xr__) */
