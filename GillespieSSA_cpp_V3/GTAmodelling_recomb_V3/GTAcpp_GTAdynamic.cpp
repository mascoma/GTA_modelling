/*GTAcpp_GTAdynamic.cpp
Jun 13, 2016
Xin
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <iostream>
#include "GTAcpp_recomb.h"

using namespace std;

long  *GTAdynamic(int *lineNum, int *popSizeArray, long *GTA1, long *GTA0){
	long temp;
	int *lysisIndex = 0;
	int *g2Index = 0;
	int *recombIndex = 0;
	int c;
	long numOfGTA1;
	long numOfGTA0;
	static long GTAarray[2];
    numOfGTA1 = *GTA1;
    numOfGTA0 = *GTA0;
	// update number of GTA
	// increase

	//printf("number of GTA, GTA0, GTA1 are %li \t %li \t %li\n", *numOfGTA, *numOfGTA0, *numOfGTA1);
	if (STRAIN_NAME.compare("Xa") != 0){
	    temp = 0;

	    	lysisIndex = lineNumOfLysisDeath(STRAIN_NAME.c_str());
	    	for (c=0; c<NUM_LYSIS; c++){
	    			if (*lineNum == *(lysisIndex+c)){
	    				//printf("found!!\n");
	    				if(STRAIN_NAME.find("Xa") != std::string::npos) {
	    					//printf("Xa included\n");
	    					temp = GTAIncrease(popSizeArray[c+2]);
	    					//printf("popSize is %d\n", popSizeArray[c]);
	    					if (c%2==0){
	    						numOfGTA1 = numOfGTA1 + temp;
	    					//printf("GTA1 increase %li\n", temp);
	    					}
	    					else{
	    						numOfGTA0 = numOfGTA0 + temp;
	    					//printf("GTA0 increase %li\n", temp);
	    					}
	    				}
	    				else{
	    					//printf("Xa not included\n");
	    					temp = GTAIncrease(popSizeArray[c]);
	    					//printf("popSize is %d\n", popSizeArray[c]);
	    					if (c%2==0){
	    						numOfGTA1 = numOfGTA1 + temp;
	    					   //	printf("GTA1 increase %li\n", *temp);
	    					}
	    					else{
	    						numOfGTA0 = numOfGTA0 + temp;
	    						//printf("GTA0 increase %li\n", *temp);
	    					}
	    				}
	    				GTAarray[0] = numOfGTA1;
	    				GTAarray[1] = numOfGTA0;
	    	    			return(GTAarray);
	    			}
	    	}
	}
	// decrease from growth2
	if (STRAIN_NAME.compare("Xa")!=0 ||STRAIN_NAME.compare("XaXr")!=0){
	    temp = 0;
	    	g2Index = lineNumOfGrowth2(STRAIN_NAME.c_str());
	    	for (c=0; c<NUM_G2; c++){
	    	if (*lineNum == *(g2Index+c)){
	    		//printf("found!!\n");
	    		temp = GTADecrease(GTA1);
	    		//printf("popSize is %d\n", popSizeArray[c]);
	    		numOfGTA1 = numOfGTA1 - temp;
	    		//printf("GTA1 decrease %li\n", *temp);
	    		temp = GTADecrease(GTA0);
	    		numOfGTA0 = numOfGTA0 - temp;
	    		//printf("GTA0 decrease %li\n", *temp);
			GTAarray[0] = numOfGTA1;
			GTAarray[1] = numOfGTA0;
			return(GTAarray);
	    		}
	    	}
	}

	// decrease from recomb
	if (STRAIN_NAME.compare("Xa")!=0||STRAIN_NAME.compare("Xn")!=0||STRAIN_NAME.compare("XaXn")!=0){
	    	temp = 0;
	    recombIndex = lineNumOfRecomb(STRAIN_NAME.c_str());
	    		for (c=0; c<NUM_RECOMB; c++){
	    		    	if (*lineNum == *(recombIndex+c)){
	        		//printf("found!!\n");
	        			if(STRAIN_NAME.find("XaXn") != std::string::npos) {
	        				//printf("XaXn included\n");
	        				//printf("popSize is %d\n", popSizeArray[c]);
	        				if (c%2==0){
	        					temp = GTADecrease(GTA1);
	        					numOfGTA1 = numOfGTA1 - temp;
	        						//printf("GTA1 decrease %li\n", *temp);
	        				}
	        				else{
	        					temp = GTADecrease(GTA0);
	        					numOfGTA0 = numOfGTA0 - temp;
	        						//printf("GTA0 decrease %li\n", *temp);
	        				}
	        			}
	        			else if(STRAIN_NAME.find("Xa") != std::string::npos||STRAIN_NAME.find("Xn") != std::string::npos){
	   					//printf("popSize is %d\n", popSizeArray[c]);
	   					if (c%2==0){
	   						temp = GTADecrease(GTA1);
	   						numOfGTA1 = numOfGTA1 - temp;
	   						// printf("GTA1 decrease %li\n", *temp);
	   					}
	   					else{
	   						temp = GTADecrease(GTA0);
	   						numOfGTA0 = numOfGTA0 - temp;
	   						//printf("GTA0 decrease %li\n", *temp);
	   					}
	        			}
	        			else{
	        					//printf("Xa or Xn not included\n");
	        					//printf("popSize is %d\n", popSizeArray[c]);
	        					if (c%2==0){
	               				temp = GTADecrease(GTA1);
	        						numOfGTA1 = numOfGTA1 - temp;
	        						//printf("GTA1 decrease %li\n", *temp);
	        					}
	        					else{
	        						temp = GTADecrease(GTA0);
	        						numOfGTA0 = numOfGTA0 - temp;
	        						//printf("GTA0 decrease %li\n", *temp);
	        					}
	        			}
	        			GTAarray[0] = numOfGTA1;
	        			GTAarray[1] = numOfGTA0;
	        			return(GTAarray);
	    		    	}
	    		}
	}
	//printf("number of GTA, GTA0, GTA1 are %li \t %li \t %li\n", *numOfGTA, *numOfGTA0, *numOfGTA1);
	return(GTAarray);
}
