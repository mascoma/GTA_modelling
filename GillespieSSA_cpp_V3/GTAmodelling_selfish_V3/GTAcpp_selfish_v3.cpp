/*GTAcpp_selfish_v3.cpp
Jun 1, 2016
Xin
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "GTAcpp_selfish.h"

//#define SIMPLE_SPRNG		/* simple interface                        */
//#include "sprng_cpp.h"          /* SPRNG header file                       */

//#define SEED 985456376
using namespace std;

long CARRYING_CAPACITY;
long INIT_TOTALPOPSIZE;
double PROP_WILDTYPE;
double GROWTH_RATE;
double UPTAKE_RATE;
double INCORPORATION_RATE;
double LYSIS;
double LARGE_HEAD_FREQ;
int BURST_SIZE;
bool CONSOLE;
int TIME_MAX;
int GENERATIONFREQ; // frequency of storing results in the output file

string OUTPUT;




int main(int argc, char *argv[])
{
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
    cc >> varName >> PROP_WILDTYPE;
    cout << varName << "\t" << PROP_WILDTYPE <<endl;


    getline(input_file, line);
    istringstream dd(line);
    dd >> varName >> GROWTH_RATE;
    cout << varName << "\t" << GROWTH_RATE <<endl;

    getline(input_file, line);
    istringstream ee(line);
    ee >> varName >> UPTAKE_RATE;
    cout << varName << "\t" << UPTAKE_RATE <<endl;


    getline(input_file, line);
    istringstream ff(line);
    ff >> varName >> INCORPORATION_RATE;
    cout << varName << "\t" << INCORPORATION_RATE <<endl;

    getline(input_file, line);
    istringstream gg(line);
    gg >> varName >> LYSIS;
    cout << varName << "\t" << LYSIS <<endl;


    getline(input_file, line);
    istringstream hh(line);
    hh >> varName >> BURST_SIZE;
    cout << varName << "\t" << BURST_SIZE <<endl;

    getline(input_file, line);
    istringstream ii(line);
    ii >> varName >> LARGE_HEAD_FREQ;
    cout << varName << "\t" << LARGE_HEAD_FREQ <<endl;

    getline(input_file, line);
    istringstream jj(line);
    jj >> varName >> CONSOLE;
    cout << varName << "\t" << CONSOLE <<endl;

    getline(input_file, line);
    istringstream kk(line);
    kk >> varName >> TIME_MAX;
    cout << varName << "\t" << TIME_MAX <<endl;


    getline(input_file, line);
    istringstream ll(line);
    ll >> varName >> GENERATIONFREQ;
    cout << varName << "\t" << GENERATIONFREQ <<endl;


    getline(input_file, line);
    istringstream mm(line);
    mm >> varName >> OUTPUT;
    cout << varName << "\t" << OUTPUT.c_str() <<endl;

    // setup reaction status matrix

    int reactStat[4][2];
    reactStat[0][0] = 1;  // growth_pop1
    reactStat[1][0] = -1; // lysis
    reactStat[2][0] = 1; // transfer +
    reactStat[3][0] = 0; // growth_pop2
    reactStat[0][1] = 0; // growth_pop1
    reactStat[1][1] = 0; // lysis
    reactStat[2][1] = -1; // transfer -
    reactStat[3][1] = 1; // growth_pop2

    // setup the init pop sizen and time

    long popSize1;
    long popSize2;
    long *popSizeArray = (long *)malloc(2*sizeof(long));
    long *totalPopSize = (long *)malloc(sizeof(long));
    double time_0 = TIME_INIT;
    double *gtime = (double *)malloc(sizeof(double));
    *gtime = time_0;

    *totalPopSize = INIT_TOTALPOPSIZE;
    popSize1 = INIT_TOTALPOPSIZE*(1 - PROP_WILDTYPE);
    popSize2 = INIT_TOTALPOPSIZE*PROP_WILDTYPE;
    popSizeArray[0] = popSize1;
    popSizeArray[1] = popSize2;

    // setup reaction array, all the reaction functions are in the reaction file

    int *reactions =  (int *)malloc(4*sizeof(int));
    int delta_growth1;
    int lyNum;
    int transferNum;
    int delta_growth2;


    delta_growth1 = growth_popSize1(popSize1, CARRYING_CAPACITY,totalPopSize, GROWTH_RATE, LYSIS);
    delta_growth2 = growth_popSize2(popSize2, CARRYING_CAPACITY, totalPopSize, GROWTH_RATE);
    lyNum = lysis(popSize1, LYSIS);
    transferNum = transfer(LYSIS, BURST_SIZE, UPTAKE_RATE, INCORPORATION_RATE, LARGE_HEAD_FREQ, popSize1, popSize2);

    reactions[0] = delta_growth1;
    reactions[1] = lyNum;
    reactions[2] = transferNum;
    reactions[3] = delta_growth2;

    printf("reactions are %d, %d, %d, %d \n", reactions[0], reactions[1], reactions[2], reactions[3]);
    // show the initial status

    printf("Starting Stochastic simulation...\n");
    printf("The system includes 2 population type(s)\n");
    printf("Generation\tTime\tGTA+\tGTA-\n");
    printf("[0]\t%e\t%li\t%li\n", time_0, popSize1, popSize2);

    printf("\n");

    //write the results in files
    FILE *fp;  // pointer to a file type
    fp = fopen(OUTPUT.c_str(), "w"); // Change to match your path
    fprintf(fp,"Generations\tTime\tGTA+\tGTA-\n");

    fprintf(fp,"0\t%e\t%li\t%li\n", time_0, popSize1, popSize2);

    fprintf(fp, "\n");
    fclose(fp);


	int l;
	int k;
	int *lineNum = (int *)malloc(sizeof(int));
    long long *n = (long long*)malloc(sizeof(long long));
    *n=1;
    int *randomArray = (int *)malloc(2*sizeof(int));
    long *g = (long *)malloc(sizeof(long));
    double *delta_tau = (double *)malloc(sizeof(double));
    long *rnum1 = (long *)malloc(sizeof(long));
    double *rnum2 = (double *)malloc(sizeof(double));
    long *reactionsum = (long *)malloc(sizeof(long));
    srand(time_t(NULL));
    //init_sprng(SEED,SPRNG_DEFAULT, DEFAULT_RNG_TYPE);
    //print_sprng();
    fp = fopen(OUTPUT.c_str(), "a"); // Change to match your pathq

    while((*gtime < TIME_MAX)&&(isNegativePopSize(popSizeArray))&&(isNegativeReaction(reactions))){
    		// random sampling the reaction status matrix row, weighted by reaction array
    		*reactionsum=0;
    		for (l = 0; l < 4; l++){
    			//printf("%li\n", array[i]);
    			*reactionsum = *reactionsum + reactions[l];
    		}
    		//printf("sum of all the reactions is %li \n", *reactionsum);
    		*rnum1 = rand();
    		//*rnum1 = isprng();

    		*rnum1 = *rnum1 % *reactionsum +1;
    		//printf("random number 1 is %li \n", *rnum1);
    		*lineNum = weightedSampling(rnum1, reactions);
    		//printf("reactions are %d, %d, %d, %d \n", reactions[0], reactions[1], reactions[2], reactions[3]);
    		//printf("random line is %d \n", *lineNum);
    		for (k=0; k < 2; k++){
    			randomArray[k]=reactStat[*lineNum][k];
    		}

    		// random sampling for tau
    		*rnum2 = (double)rand()/(double)RAND_MAX;

    		*delta_tau = calcTau(reactions, rnum2);

    		// update population size and tau

    		int t;
    		for (t=0; t < 2; t++){
    			popSizeArray[t] = popSizeArray[t] +randomArray[t];
    		}



    		*totalPopSize = popSizeArray[0] + popSizeArray[1];


    		*delta_tau = calcTau(reactions, rnum2);
    		//printf("changes of tau is: %e\n", delta_tau);
    		*gtime = *gtime + *delta_tau;


    		// turn on/off CONSOLE
    		if (CONSOLE){

    			if (*n % GENERATIONFREQ == 0){
    				*g = (long)(*n/GENERATIONFREQ);
    				//fprintf(fp, "Generation\tTime\t\tPopulationsize of %d Genotype \n", genotypeNum);
    				printf("[%li]\ttime[%li]= %e\t%li\t%li\n",*g, *g, *gtime, popSizeArray[0], popSizeArray[1]);
    			}
    		}

    		//write the results in files
    		if (*n % GENERATIONFREQ == 0){
    			*g = (long)(*n/GENERATIONFREQ);
    			//fprintf(fp, "Generation\tTime\t\tPopulationsize of %d Genotype \n", genotypeNum);
    			fprintf(fp, "[%li]\ttime[%li]= %e\t%li\t%li\n",*g, *g, *gtime, popSizeArray[0], popSizeArray[1]);
    		}

    		// update reaction array
    		popSize1 = popSizeArray[0];
    		popSize2 = popSizeArray[1];
    	    delta_growth1 = growth_popSize1(popSize1, CARRYING_CAPACITY,totalPopSize, GROWTH_RATE, LYSIS);
    	    delta_growth2 = growth_popSize2(popSize2, CARRYING_CAPACITY, totalPopSize, GROWTH_RATE);
    	    lyNum = lysis(popSize1, LYSIS);
    	    transferNum = transfer(LYSIS, BURST_SIZE, UPTAKE_RATE, INCORPORATION_RATE, LARGE_HEAD_FREQ, popSize1, popSize2);

    	    reactions[0] = delta_growth1;
    	    reactions[1] = lyNum;
    	    reactions[2] = transferNum;
    	    reactions[3] = delta_growth2;
    		*n = *n+1;
    	}
    fclose(fp);
    return(0);


}
