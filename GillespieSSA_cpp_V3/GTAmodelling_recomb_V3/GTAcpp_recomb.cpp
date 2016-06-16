/*test.cpp
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
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end

#include "GTAcpp_recomb.h"

// read file to load the parameter values
int CARRYING_CAPACITY;
int INIT_TOTALPOPSIZE;
double PROP_WILDTYPE; // genotype 0 is wildtype
double GROWTH_RATE;
double MUTATION;
double NUTRITIONAL_VALUE;
double UPTAKE_RATE;
double INCORPORATION_RATE;
double LYSIS;
double SELECTION;
double FITNESS;
int BURST_SIZE;
double NUM_LOCI;
bool CONSOLE;
int TIME_MAX;
string OUTPUT;
// parameters built the system
int NUM_STRAIN;
double PROP_STRAIN1;
double PROP_STRAIN2;
double PROP_STRAIN3;
string STRAIN_NAME;
int TOTAL_NUM_REACTIONS;
int TOTAL_NUM_POPULATION;
int NUM_LYSIS;
int NUM_G2;
int NUM_RECOMB;
int SAMPLING_FREQ;  // frequency of storing results in the output file

//#define SIMPLE_SPRNG		/* simple interface                        */
//#include "sprng_cpp.h"          /* SPRNG header file                       */

//#define SEED 985456376
using namespace std;

int main(int argc, char *argv[]){
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
     ee >> varName >> MUTATION;
     cout << varName << "\t" << MUTATION <<endl;

     getline(input_file, line);
     istringstream ff(line);
     ff >> varName >> NUTRITIONAL_VALUE;
     cout << varName << "\t" << NUTRITIONAL_VALUE <<endl;

     getline(input_file, line);
     istringstream gg(line);
     gg >> varName >> UPTAKE_RATE;
     cout << varName << "\t" << UPTAKE_RATE <<endl;

     getline(input_file, line);
     istringstream hh(line);
     hh >> varName >> INCORPORATION_RATE;
     cout << varName << "\t" << INCORPORATION_RATE <<endl;

     getline(input_file, line);
     istringstream ii(line);
     ii >> varName >> LYSIS;
     cout << varName << "\t" << LYSIS <<endl;

     getline(input_file, line);
     istringstream jj(line);
     jj >> varName >> SELECTION;
     cout << varName << "\t" << SELECTION <<endl;

     getline(input_file, line);
     istringstream kk(line);
     kk >> varName >> FITNESS;
     cout << varName << "\t" << FITNESS <<endl;

     getline(input_file, line);
     istringstream ll(line);
     ll >> varName >> BURST_SIZE;
     cout << varName << "\t" << BURST_SIZE <<endl;

     getline(input_file, line);
     istringstream mm(line);
     mm >> varName >> NUM_LOCI;
     cout << varName << "\t" << NUM_LOCI <<endl;

     getline(input_file, line);
     istringstream nn(line);
     nn >> varName >> CONSOLE;
     cout << varName << "\t" << CONSOLE <<endl;

     getline(input_file, line);
     istringstream oo(line);
     oo >> varName >> TIME_MAX;
     cout << varName << "\t" << TIME_MAX <<endl;

     getline(input_file, line);
     istringstream pp(line);
     pp >> varName >> NUM_STRAIN;
     cout << varName << "\t" << NUM_STRAIN <<endl;

     getline(input_file, line);
     istringstream qq(line);
     qq >> varName >> PROP_STRAIN1;
     cout << varName << "\t" << PROP_STRAIN1 <<endl;

     getline(input_file, line);
     istringstream rr(line);
     rr >> varName >> PROP_STRAIN2;
     cout << varName << "\t" << PROP_STRAIN2 <<endl;

     getline(input_file, line);
     istringstream ss(line);
     ss >> varName >> PROP_STRAIN3;
     cout << varName << "\t" << PROP_STRAIN3 <<endl;

     getline(input_file, line);
     istringstream tt(line);
     tt >> varName >> STRAIN_NAME;
     cout << varName << "\t" << STRAIN_NAME <<endl;

     getline(input_file, line);
     istringstream uu(line);
     uu >> varName >> TOTAL_NUM_REACTIONS;
     cout << varName << "\t" << TOTAL_NUM_REACTIONS <<endl;

     getline(input_file, line);
     istringstream vv(line);
     vv >> varName >> TOTAL_NUM_POPULATION;
     cout << varName << "\t" << TOTAL_NUM_POPULATION <<endl;

     getline(input_file, line);
     istringstream ww(line);
     ww >> varName >> NUM_LYSIS;
     cout << varName << "\t" << NUM_LYSIS <<endl;

     getline(input_file, line);
     istringstream xx(line);
     xx >> varName >> NUM_G2;
     cout << varName << "\t" << NUM_G2 <<endl;

     getline(input_file, line);
     istringstream yy(line);
     yy >> varName >> NUM_RECOMB;
     cout << varName << "\t" << NUM_RECOMB <<endl;

     getline(input_file, line);
     istringstream zz(line);
     zz >> varName >> SAMPLING_FREQ;
     cout << varName << "\t" << SAMPLING_FREQ <<endl;

     getline(input_file, line);
     istringstream output(line);
     output >> varName >> OUTPUT;
     cout << varName << "\t" << OUTPUT.c_str() <<endl;

    // define the variables
    // loop index
	int a;
	int b;
	int c;
    int d;
    int e;
	int l;
	int k;
	int t;
	int u;
	int status;
    long *n = (long*)malloc(sizeof(long));
    bool condition1;
    bool condition2;
    bool condition3;

    // dynamic of number of GTA
	long *GTAIndex = (long *)malloc(sizeof(long));
	long *numOfGTA = (long *)malloc(sizeof(long));
	long *numOfGTA1 = (long *)malloc(sizeof(long));
	long *numOfGTA0 = (long *)malloc(sizeof(long));
	long *temp = (long *)malloc(sizeof(long));

	int *lineNum = (int *)malloc(sizeof(int));
    int *randomArray = (int *)malloc(TOTAL_NUM_POPULATION*sizeof(int));
    long *g = (long *)malloc(sizeof(long));
    double *delta_tau = (double *)malloc(sizeof(double));
    long *rnum1 = (long *)malloc(sizeof(long));
    double *rnum2 = (double *)malloc(sizeof(double));
    int *reactionsum = (int *)malloc(sizeof(int));
    double time_0;
    double *gtime = (double *)malloc(sizeof(double));
	int *totalPopSize = (int *)malloc(sizeof(int));
	int *popSizeArray = (int *)malloc(TOTAL_NUM_POPULATION*sizeof(int));
	int *popSizeIndex = (int *)malloc(sizeof(int));
	int *reactionIndex = (int *)malloc(sizeof(int));
	int *reactions =  (int *)malloc(TOTAL_NUM_REACTIONS*sizeof(int));

	//initialize population size and time
	*totalPopSize = INIT_TOTALPOPSIZE;

	popSizeIndex = initPopSize(INIT_TOTALPOPSIZE);

	for(e=0; e<TOTAL_NUM_POPULATION; e++){
		popSizeArray[e] = *(popSizeIndex+e);
		//printf("pop size [%d] \t %d \n", e, popSizeArray[e]);
	}

	time_0 = TIME_INIT;
    *gtime = time_0;
    *numOfGTA = 0;
    *numOfGTA1 = 0;
    *numOfGTA0 = 0;
    // setup reaction array, all the reaction functions are in the reaction file

	reactionIndex = reactArray(popSizeArray, totalPopSize, LYSIS, FITNESS, numOfGTA1, numOfGTA0, numOfGTA, STRAIN_NAME.c_str());
	for (a=0; a<TOTAL_NUM_REACTIONS; a++){
		reactions[a] = *(reactionIndex+a);
	    //printf("raection[%d]: \t %d\n", a, reactions[a]);
	}

    // generate the reaction status matrix
    statusMatrix reactionStatus = getArray(STRAIN_NAME.c_str());
    /*
    for(a=0; a<TOTAL_NUM_REACTIONS; a++){
    		for(d=0; d<TOTAL_NUM_POPULATION; d++){
    			status = reactionStatus.reactStat[a][d];
    			printf("status is %d\n", status);
    		}
    }
	*/

    printf("Starting Stochastic simulation...\n");
    printf("The system includes %d strain(s) and %d genotypes \n", NUM_STRAIN, TOTAL_NUM_POPULATION);
    printf("Generation\tTime\tG1S1\tG0S1\tG1S2\tG0S2\tG1S3\tG0S3\tG1S4\tG0S4\n");
    printf("[0]\t%e\t", time_0);
    for(b=0; b<TOTAL_NUM_POPULATION; b++){
    		printf("%d\t\t", popSizeArray[b]);
    }
    printf("\n");

    //write the results in files
    FILE *fp;  // pointer to a file type
    fp = fopen(OUTPUT.c_str(), "w"); // Change to match your path

    fprintf(fp, "[0]\t%e\t", time_0);
    for(b=0; b<TOTAL_NUM_POPULATION; b++){
    		fprintf(fp, "%d\t", popSizeArray[b]);
    }
    fprintf(fp, "\n");

    fclose(fp);

    srand(time_t(NULL));
    //init_sprng(SEED,SPRNG_DEFAULT, DEFAULT_RNG_TYPE);
    //print_sprng();
    fp = fopen(OUTPUT.c_str(), "a"); // Change to match your pathq
    *n=1;

    while((*gtime < TIME_MAX)&& isNegativePopSize(popSizeArray) && isNegativeReaction(reactions)){
    		// random sampling the reaction status matrix row, weighted by reaction array
    		*reactionsum=0;
    		for (l = 0; l < TOTAL_NUM_REACTIONS; l++){

    			*reactionsum = *reactionsum + reactions[l];
    		}
    		*rnum1 = rand();
    		//*rnum1 = isprng();

    		*rnum1 = *rnum1 % *reactionsum +1;
    		//printf("random number 1 is %li \n", *rnum1);
    		*lineNum = weightedSampling(rnum1, reactions);

    		//printf("random line is %d \n", *lineNum);


    		for (k=0; k < TOTAL_NUM_POPULATION; k++){
    			randomArray[k]=reactionStatus.reactStat[*lineNum][k];
    		}

    		// random sampling for tau
    		*rnum2 = (double)rand()/(double)RAND_MAX;
    		//printf("random number 2 is %f \n", *rnum2);

    		// update population size and tau

    		for (t=0; t<TOTAL_NUM_POPULATION; t++){
    			popSizeArray[t] = popSizeArray[t] +randomArray[t];
    			if (popSizeArray[t] < 0){
    				popSizeArray[t] = 0;
    			}
    			//printf("popsize[%d]: %d \t %d\n", t, popSizeArray[t], randomArray[t]);
    		}

    		*temp = 0;
    		for (u=0; u<TOTAL_NUM_POPULATION; u++){
    			*temp = *temp + popSizeArray[u];
    		}
    		*totalPopSize = *temp;
    		//printf("pop size and total are %d \t %d \t %d \t %d \t %d \n", popSizeArray[0], popSizeArray[1], popSizeArray[2], popSizeArray[3], *totalPopSize);

    		*delta_tau = calcTau(reactions, rnum2);
    		//printf("changes of tau is: %e\n", *delta_tau);
    		*gtime = *gtime + *delta_tau;


    		// turn on/off CONSOLE
    		if (CONSOLE){

    			if (*n % SAMPLING_FREQ == 0){
    				*g = (long)(*n/SAMPLING_FREQ);

    			    printf("[%li]\t%e\t", *g, *gtime);
    			    for(b=0; b<TOTAL_NUM_POPULATION; b++){
    			    		printf("%d\t\t", popSizeArray[b]);
    			    }
    			    printf("\n");
    			}
    		}

    		//write the results in files
    		if (*n % SAMPLING_FREQ == 0){
    			*g = (long)(*n/SAMPLING_FREQ);

			fprintf(fp, "[%li]\t%e\t", *g, *gtime);
			for(b=0; b<TOTAL_NUM_POPULATION; b++){
			     fprintf(fp, "%d\t\t", popSizeArray[b]);
			}
			fprintf(fp, "\n");
    		}

    		// update number of GTA
    		GTAIndex = GTAdynamic(lineNum, popSizeArray, numOfGTA1, numOfGTA0);
    		*numOfGTA1 = *GTAIndex;
    		*numOfGTA0 = *(GTAIndex + 1);
    		*numOfGTA = *numOfGTA1 + *numOfGTA0;

    		// update reaction array
    		reactionIndex = reactArray(popSizeArray, totalPopSize, LYSIS, FITNESS, numOfGTA1, numOfGTA0, numOfGTA, STRAIN_NAME.c_str());

    		for (d=0; d<TOTAL_NUM_REACTIONS; d++){
    			reactions[d] = *(reactionIndex+d);
    			//printf("reaction[%d]: %d\n",d, *(reactionIndex+d));
    		}
    		*n = *n+1;
    	    //condition1 = (*gtime < TIME_MAX);
    	    //condition2 = isNegativePopSize(popSizeArray, TOTAL_NUM_POPULATION);
    	    //condition3 = isNegativeReaction(reactions, TOTAL_NUM_REACTIONS);
    	    //printf("condition1: \t %d; condition2: \t %d; condition3: \t %d\n", condition1, condition2, condition3);
    	}

    fclose(fp);
    return(0);
}
