//
//  GTAcpp_Xr.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include "GTAcpp_G.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#define SIMPLE_SPRNG		/* simple interface                        */
#include "sprng_cpp.h"          /* SPRNG header file                       */

#define SEED 985456376
using namespace std;

long CARRYING_CAPACITY;
long INIT_TOTALPOPSIZE;
int NUM_LOCI;
double PROP_WILDTYPE;
double INTRI_DEATH_RATE;
double MUTATION;
double SELECTION;
double GROWTH_RATE;
double UPTAKE_RATE;
double NUTRITIONAL_VALUE;
double INCORPORATION_RATE;
double LYSIS;
bool CONSOLE;
int TIME_MAX;
int NUM_OF_POPTYPE;
string OUTPUT;
string REACTIONSYSTEM;
int genotypeNum;
int genetypeNum;
int numOfgeonReact;
int numofgeneReact;


int main(int argc, char *argv[]){
    //loop counter
    int a=0;
    int b=0;
    int e=0;
    int f=0;
    int h=0;
    int i=0;
    int j=0;
    int k=0;
    int l=0;
    int m=0;
    int o=0;
    int r=0;
    int s=0;
    int t=0;
    //int u=0;
    int v=0;
    int w=0;
    int x=0;
    int z=0;
    // read file to load the parameter values
    
    string INPUT;
    string varName;
    string line;
    INPUT = argv[1];
    ifstream input_file(INPUT.c_str());
    if(!input_file.is_open()){
        cout <<"error, no such file!"<<endl;
        return 1;
    }
    cout << "input parameters as follow:" <<endl;
    getline(input_file, line);
    istringstream aa(line);
    aa >> varName >> CARRYING_CAPACITY;
    cout << varName << "\t" << CARRYING_CAPACITY << endl;
    
    getline(input_file, line);
    istringstream bb(line);
    bb >> varName >> INIT_TOTALPOPSIZE;
    cout << varName << "\t" << INIT_TOTALPOPSIZE <<endl;
    
    getline(input_file, line);
    istringstream cc(line);
    cc >> varName >> NUM_LOCI;
    cout << varName << "\t" << NUM_LOCI <<endl;
    
    getline(input_file, line);
    istringstream rr(line);
    rr >> varName >> PROP_WILDTYPE;
    cout << varName << "\t" << PROP_WILDTYPE <<endl;

    getline(input_file, line);
    istringstream dd(line);
    dd >> varName >> INTRI_DEATH_RATE;
    cout << varName << "\t" << INTRI_DEATH_RATE <<endl;
    
    getline(input_file, line);
    istringstream ee(line);
    ee >> varName >> MUTATION;
    cout << varName << "\t" << MUTATION <<endl;
    
    getline(input_file, line);
    istringstream ff(line);
    ff >> varName >> SELECTION;
    cout << varName << "\t" << SELECTION <<endl;
    
    getline(input_file, line);
    istringstream gg(line);
    gg >> varName >> GROWTH_RATE;
    cout << varName << "\t" << GROWTH_RATE <<endl;
    
    getline(input_file, line);
    istringstream hh(line);
    hh >> varName >> UPTAKE_RATE;
    cout << varName << "\t" << UPTAKE_RATE <<endl;

    getline(input_file, line);
    istringstream jj(line);
    jj >> varName >> NUTRITIONAL_VALUE;
    cout << varName << "\t" << NUTRITIONAL_VALUE <<endl;
    
    getline(input_file, line);
    istringstream kk(line);
    kk >> varName >> INCORPORATION_RATE;
    cout << varName << "\t" << INCORPORATION_RATE <<endl;
    
    getline(input_file, line);
    istringstream ll(line);
    ll >> varName >> LYSIS;
    cout << varName << "\t" << LYSIS <<endl;
    
    getline(input_file, line);
    istringstream mm(line);
    mm >> varName >> CONSOLE;
    cout << varName << "\t" << CONSOLE <<endl;
    
    getline(input_file, line);
    istringstream nn(line);
    nn >> varName >> TIME_MAX;
    cout << varName << "\t" << TIME_MAX <<endl;
    
    getline(input_file, line);
    istringstream oo(line);
    oo >> varName >> NUM_OF_POPTYPE;
    cout << varName << "\t" << NUM_OF_POPTYPE <<endl;
    
    getline(input_file, line);
    istringstream pp(line);
    pp >> varName >> OUTPUT;
    cout << varName << "\t" << OUTPUT.c_str() <<endl;
    
    getline(input_file, line);
    istringstream qq(line);
    qq >> varName >> REACTIONSYSTEM;
    cout << varName << "\t" << REACTIONSYSTEM.c_str() <<endl;
    
    genotypeNum =(int)pow(2, NUM_LOCI);
    cout << "genotypeNum" << "\t" << genotypeNum <<endl;
    genetypeNum = 2*NUM_LOCI;
    cout << "genetypeNum" << "\t" << genetypeNum <<endl;
    numOfgeonReact = genotypeNum * NUM_OF_PROCESSES * NUM_OF_POPTYPE;
    cout << "numOfgeonReact" << "\t" << numOfgeonReact <<endl;
    numofgeneReact = genetypeNum * NUM_OF_PROCESSES_f;
    cout << "numofgeneReact" << "\t" << numofgeneReact <<endl;

    // initialize population size
    long *geno_gene_0 = (long*)malloc((genotypeNum*NUM_OF_POPTYPE+genetypeNum)*sizeof(long));
    double p = 0;
    long *popSizeArray = (long *)malloc((genotypeNum*NUM_OF_POPTYPE+genetypeNum)*sizeof(long));
    long *totalPopSize = (long *)malloc(sizeof(long));
    long *temp1 = (long *)malloc(sizeof(long));
    *temp1=0;
    for (z= 0; z< (genotypeNum*NUM_OF_POPTYPE+genetypeNum); z++){
        geno_gene_0[z] = 0;
    }
    
    for (a =0; a < NUM_OF_POPTYPE; a++) {
        geno_gene_0[genotypeNum*a] = (INIT_TOTALPOPSIZE*PROP_WILDTYPE)/NUM_OF_POPTYPE;
    
    }
    
    p = (double)(1-PROP_WILDTYPE)/(NUM_LOCI*NUM_OF_POPTYPE);
    for (b=0; b < NUM_OF_POPTYPE; b++){
        for (w=0; w< NUM_LOCI; w++){
            geno_gene_0[((int)pow(2,w)+b*genotypeNum)] = INIT_TOTALPOPSIZE*p;
        }
    }
    
  
    for (s=0; s < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); s++){
        popSizeArray[s] = geno_gene_0[s];
    }
    for (v=0; v < (genotypeNum*NUM_OF_POPTYPE); v++){
        *temp1 = *temp1 + popSizeArray[v];
    }
    *totalPopSize = *temp1;
    //printf("initial pop size is: %li\n", *totalPopSize);
    
    // initialize time
    double time_0 = TIME_INIT;
    double *gtime = (double *)malloc(sizeof(double));
    *gtime = time_0;
    
         //int aa, bb;
    //for(aa = 0; aa < numOfReactions; ++aa){
    //		for(bb = 0; bb < (genotypeNum+genetypeNum); ++bb){
    //			printf("reactionStatus is %d\n", reactionsStat[aa].reactionStatus[bb]);
    //	}
    //	}
    
    // generate the reaction status matrix
    statusMatrix reactionStatus = getArray();

    /*
    int c, d;
    for (c=0; c <(numOfgeonReact+numofgeneReact); c++){
    		   printf("[%d]\n", c);
           for (d=0; d < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); d++){
               printf("%d\t",reactionStatus.reactStat[c][d]);
               //if(reactionStatus.reactStat[c][55]== 1){printf("reactions[%d][%d] = 1\n", c, d);}

           }
           printf("\n");
       }
     */
    // generate a matrix to store all the genotypes (** genoetypePool)
    int **genotypePool = (int **)malloc(sizeof *genotypePool * genotypeNum);
    if (genotypePool){
        for (e = 0; e < genotypeNum; e++){
            genotypePool[e] = (int *)malloc(sizeof *genotypePool[e] * NUM_LOCI);
        }
    }
    
    for (j=0; j < genotypeNum; j++){
        f = j;
        for (h=0; h < NUM_LOCI; h++){
            genotypePool[j][h] = f%2;
            //printf("%d\n", f);
            f = f/2;
            //printf("genotype[%d][%d]: %d\t",j,h, genotypePool[j][h]);
        }
        //printf("\n");
    }
    //printf("genotype[%d][%d]: %d\n",0,0, genotypePool[0][0]);

   	int *reactionIndex;
    //int numOfReactions = NUM_OF_PROCESSES_Xr * genotypeNum + genetypeNum;
    //int *reactionArray_0 = (int *) malloc((numOfgeonReact+numofgeneReact)*sizeof(int));
    int *reactionArray =  (int *)malloc((numOfgeonReact+numofgeneReact)*sizeof(int));
    reactionIndex = reactArray(popSizeArray, totalPopSize, genotypePool, REACTIONSYSTEM);
    for (i=0; i < (numOfgeonReact+numofgeneReact); i++){
        //printf("There are %d of reactions!\n", numOfReactions);
        //printf("reaction_0[%d]: %li\n",i, *(reactionIndex+i));
    		reactionArray[i] = *(reactionIndex+i);
    }
    /*
    for(u=0; u < (numOfgeonReact+numofgeneReact); u++){
        reactionArray[u]=reactionArray_0[u];
        //printf("reactionArray[%d]: %d\n", u, reactionArray[u]);
    }
    */
    // show the initial status
    
    printf("Starting Stochastic simulation...\n");
    printf("The system includes %d population type(s) %s\n", NUM_OF_POPTYPE, REACTIONSYSTEM.c_str());
    printf("Generation\tTime\t\tPopulationsize of %d Genotype and %d Genetype \n", genotypeNum, genetypeNum);
    printf("[0]\ttime[0]= %e\t", time_0);
    for (r=0; r < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); r++){
        printf("\t%li\t", geno_gene_0[r]);
    }
    printf("\n");
    
    string header;
    header = headerName(REACTIONSYSTEM.c_str());
    //write the results in files
    FILE *fp;  // pointer to a file type
    fp = fopen(OUTPUT.c_str(), "w"); // Change to match your path
    fprintf(fp,"Generations\tTime\t");
    fprintf(fp, "%s\n",header.c_str());
    fprintf(fp,"0\t%e\t", time_0);
    for (r=0; r < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); r++){
        fprintf(fp,"%li\t", geno_gene_0[r]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    
    free(geno_gene_0);
    //free(reactionArray_0);
    
    // start stochastic loop
    
    long long *n = (long long*)malloc(sizeof(long long));
    *n=1;
    int *randomArray = (int *)malloc((genotypeNum*NUM_OF_POPTYPE+genetypeNum)*sizeof(int));
    long *g = (long *)malloc(sizeof(long));
    int *lineNum = (int *)malloc(sizeof(int));
    double *delta_tau = (double *)malloc(sizeof(double));
    long *rnum1 = (long *)malloc(sizeof(long));
    double *rnum2 = (double *)malloc(sizeof(double));
    long *reactionsum = (long *)malloc(sizeof(long));
    srand(time_t(NULL));
    //init_sprng(SEED,SPRNG_DEFAULT, DEFAULT_RNG_TYPE);
    //print_sprng();
    fp = fopen(OUTPUT.c_str(), "a"); // Change to match your pathq
    while((*gtime < TIME_MAX)&&(isNegativePopSize(popSizeArray))&&(isNegativeReaction(reactionArray))){

        // random number generator
        *reactionsum=0;
        for (l =0; l < (numOfgeonReact+numofgeneReact); l++){
            //printf("%li\n", array[i]);
            *reactionsum = *reactionsum + reactionArray[l];
        }
        //printf("reaction sum is %li\n", *reactionsum);
        //*rnum1 = rand();
        *rnum1 = isprng();
        *rnum1 = *rnum1 % *reactionsum +1;
        //*rnum2 = (double)rand()/(double)RAND_MAX;
        *rnum2 = sprng();
        //printf("%li\n", *rnum1);
        //printf("%lf\n", *rnum2);
        // ssa.d algorithm to return selected reaction status array
        *lineNum = weightedSampling(rnum1, reactionArray);
        //if (*lineNum == 172||*lineNum == 284)
        //{printf("selected lineNum is %d\n", *lineNum);}
        for (k=0; k < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); k++){
            randomArray[k]=reactionStatus.reactStat[*lineNum][k];
            //if(reactionStatus.reactStat[*lineNum][55]== 1){printf("reactions[55] = 1\n");}
            //if(randomArray[55]== 1){printf("reactions[55] = 1\n");}

        }
        // ssa.d algorithm to return tau
        *delta_tau = calcTau(reactionArray, rnum2);
        //printf("changes of tau is: %e\n", delta_tau);
        *gtime = *gtime + *delta_tau;
        
        // update the popsize array and time
        for (t=0; t < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); t++){
            popSizeArray[t] = popSizeArray[t] +randomArray[t];
        }
        *temp1=0;
        for (o=0; o < genotypeNum*NUM_OF_POPTYPE; o++){
            *temp1 = *temp1 + popSizeArray[o];
        }
        *totalPopSize = *temp1;
        
        if (CONSOLE){
        
            if (*n % GENERATIONFREQ == 0){
                *g = (long)(*n/GENERATIONFREQ);
                //fprintf(fp, "Generation\tTime\t\tPopulationsize of %d Genotype \n", genotypeNum);
                printf("[%li]\ttime[%li]= %e\t",*g, *g, *gtime);
                for (m=0; m < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); m++){
                    printf("\t%li\t", popSizeArray[m]);
                }
                printf("\n");
            }
        }
        
        //write the results in files
        if (*n % GENERATIONFREQ == 0){
    		    		*g = (long)(*n/GENERATIONFREQ);
            //fprintf(fp, "Generation\tTime\t\tPopulationsize of %d Genotype \n", genotypeNum);
            fprintf(fp,"%li\t%e\t",*g, *gtime);
            for (m=0; m < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); m++){
                fprintf(fp,"%li\t", popSizeArray[m]);
            }
            fprintf(fp, "\n");
        }
        
    	   // update reaction array
        int *reactionupdate;//pointers to index new reaction values
        reactionupdate = reactArray(popSizeArray, totalPopSize, genotypePool, REACTIONSYSTEM.c_str());
        for (x=0; x < (numOfgeonReact+numofgeneReact); x++){
            //printf("update reaction[%d]: %d\n", x, *(reactionupdate+x));
            reactionArray[x]=*(reactionupdate+x);
        }
        *n = *n+1;
        //printf("%lli\n", n);
    }
    
    //	int p, q;
    //p = isNegativePopSize(popSizeArray);
    //printf("p = %d\n", p);
    //q = isNegativeReaction(reactionArray);
    //printf("q = %d\n", q);
    fclose(fp);
    return(0);
}
