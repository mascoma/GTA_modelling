/*reactstatus.cpp
Jun 2, 2016
Xin
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <iostream>
#include <sstream>
#include "GTAcpp_recomb.h"
using namespace std;

statusMatrix getArray(string index){ // get reactioin status matrix for one strain
    int a=0;
    int b=0;
    int c=0;
    int d=0;
    int e=0;
    int f=0;

    int x1 = 0;
    int x2 = 0;
    int x3 = 0;
    int x4 = 0;

    int reactionsysIndex = 0;
    statusMatrix r;

    reactionsysIndex = reactionSystem(index);

    switch (reactionsysIndex) {
        case 1: // Xa
            for (a=0; a < 8; a++){
                for (b=0; b < 2; b++){
                     r.reactStat[a][b]=0;
                }
            }
            x1 = 4; // num of processes in each strain

            	for (c = 0; c < 2; c++) {
            		r.reactStat[(x1*c)][c] = 1 ; // growth
            		r.reactStat[(x1*c+1)][c] = -1; // death
            		r.reactStat[(x1*c+2)][c] = 1; // positive mutation
            		r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
            		r.reactStat[(x1*c+3)][c] = -1; // negative mutation
            		r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            	}
            break;
            /*
             // status changes of gene fragments influenced by lysis
             for(colNum = 0; colNum < genotypeNum; colNum++){
             for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
             y = matrixIndex(rowNum, colNum)-1;
             //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
             reactionsStat[((colNum*NUM_OF_PROCESSES_Xr)+4)].reactionStatus[genotypeNum +y] = 1;
             }
             }

             for(d=genotypeNum; d < (genotypeNum+genetypeNum); d++){
             reactionsStat[((genotypeNum*NUM_OF_PROCESSES)+(d-genotypeNum))].reactionStatus[d]=0; // gene fragement taken up
             }
             */

        case 2: // Xn
            for (a=0; a < 10; a++){
                for (b=0; b < 2; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 5; // num of processes in each strain

            	for (c = 0; c < 2; c++) {
            		r.reactStat[(x1*c)][c] = 1 ; // growth2
            		r.reactStat[(x1*c+1)][c] = -1; // death
            		r.reactStat[(x1*c+2)][c] = 1; // positive mutation
            		r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
            		r.reactStat[(x1*c+3)][c] = -1; // negative mutation
            		r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            		r.reactStat[(x1*c+4)][c] = -1; // lysis
            	}
            break;

        case 3: // Xr
            for (a=0; a < 14; a++){
                for (b=0; b < 2; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 7; // num of processes in each strain

            	for (c = 0; c < 2; c++) {
            		r.reactStat[(x1*c)][c] = 1 ; // growth
            		r.reactStat[(x1*c+1)][c] = -1; // death
            		r.reactStat[(x1*c+2)][c] = 1; // positive mutation
            		r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
            		r.reactStat[(x1*c+3)][c] = -1; // negative mutation
            		r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+4)][c] = 1; // positive recomb
        			r.reactStat[(x1*c+4)][(1-c)] = -1; //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+5)][c] = -1; // negative recomb
        			r.reactStat[(x1*c+5)][(1-c)] = 1; // negative recomb has opposite effect to the opposite genotype
            		r.reactStat[(x1*c+6)][c] = -1; // lysis
            	}
            break;


        case 4: // Xnr
            for (a=0; a < 14; a++){
                for (b=0; b < 2; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 7; // num of processes in each strain

            for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth2
        			r.reactStat[(x1*c+1)][c] = -1; // death
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+4)][c] = 1; // positive recomb
        			r.reactStat[(x1*c+4)][(1-c)] = -1; //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+5)][c] = -1; // negative recomb
        			r.reactStat[(x1*c+5)][(1-c)] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+6)][c] = -1; // lysis
            }
            break;


        case 5: // XaXn

            for (a=0; a < 18; a++){
                 for (b=0; b < 4; b++){
                      r.reactStat[a][b]=0;
                 }
             }

             x1 = 4;
             x2 = 5;

             for (c = 0; c < 2; c++) {
         			r.reactStat[(x1*c)][c] = 1 ; // growth
         			r.reactStat[(x1*c+1)][c] = -1; // death
         			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
         			r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
         			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
         			r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
             }

             for (d = 0; d < 2; d++) {
         			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
         			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
         			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
         			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
         			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
         			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
         			r.reactStat[(x1*2+x2*d+4)][(d+2)] = -1; // lysis
         		}

             break;



        case 6: //XaXr

            for (a=0; a < 22; a++){
                for (b=0; b < 4; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 4;
            x2 = 7;

            for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth
        			r.reactStat[(x1*c+1)][c] = -1; // death
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            }

            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}

            break;


        case 7: //XaXnr

            for (a=0; a < 22; a++){
                for (b=0; b < 4; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 4;
            x2 = 7;

            for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth
        			r.reactStat[(x1*c+1)][c] = -1; // death
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            }

            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}

            break;


        case 8: //XnXr

            for (a=0; a < 24; a++){
                for (b=0; b < 4; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 5;
            x2 = 7;

            for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth2
        			r.reactStat[(x1*c+1)][c] = -1; // ceath
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][1-c] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][1-c] = 1; // negative mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+4)][c] = -1; // lysis
        		}


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}

            break;


        case 9: //XnXnr

            for (a=0; a < 24; a++){
                for (b=0; b < 4; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 5;
            x2 = 7;

            for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth2
        			r.reactStat[(x1*c+1)][c] = -1; // ceath
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][1-c] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][1-c] = 1; // negative mutation has opposite effect to the opposite genotype

        			r.reactStat[(x1*c+4)][c] = -1; // lysis
        		}


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype

        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}

            break;


        case 10: //XrXnr

            for (a=0; a < 28; a++){
                for (b=0; b < 4; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 7;
            x2 = 7;

        		for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth
        			r.reactStat[(x1*c+1)][c] = -1; // death
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+4)][c] = 1; // positive recomb
        			r.reactStat[(x1*c+4)][(1-c)] = -1; //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+5)][c] = -1; // negative recomb
        			r.reactStat[(x1*c+5)][(1-c)] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+6)][c] = -1; // lysis
        		}


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype

        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}

            break;


        case 11: //XaXnXr

            for (a=0; a < 32; a++){
                for (b=0; b < 6; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 4;
            x2 = 5;
            x3 = 7;

            for (c = 0; c < 2; c++) {
                	r.reactStat[(x1*c)][c] = 1 ; // growth
                	r.reactStat[(x1*c+1)][c] = -1; // death
                	r.reactStat[(x1*c+2)][c] = 1; // positive mutation
                	r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
                	r.reactStat[(x1*c+3)][c] = -1; // negative mutation
                	r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            }


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype

        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = -1; // lysis
        		}


        		for (e = 0; e < 2; e++) {
        			r.reactStat[(x1*2+x2*2+x3*e)][(e+4)] = 1 ; // growth
        			r.reactStat[(x1*2+x2*2+x3*e+1)][(e+4)] = -1; // death
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(e+4)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(9-(e+4))] = -1; // positive mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(e+4)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(9-(e+4))] = 1; // negative mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(e+4)] = 1; // positive reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(9-(e+4))] = -1; //  positive reeomb has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(e+4)] = -1; // negative reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(9-(e+4))] = 1; // negative reeomb has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+6)][(e+4)] = -1; // lysis
        		}
            break;



        case 12: //XaXnXnr

            for (a=0; a < 32; a++){
                for (b=0; b < 6; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 4;
            x2 = 5;
            x3 = 7;

            for (c = 0; c < 2; c++) {
                	r.reactStat[(x1*c)][c] = 1 ; // growth
                	r.reactStat[(x1*c+1)][c] = -1; // death
                	r.reactStat[(x1*c+2)][c] = 1; // positive mutation
                	r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
                	r.reactStat[(x1*c+3)][c] = -1; // negative mutation
                	r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            }


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype

        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = -1; // lysis
        		}


        		for (e = 0; e < 2; e++) {
        			r.reactStat[(x1*2+x2*2+x3*e)][(e+4)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*2+x3*e+1)][(e+4)] = -1; // death
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(e+4)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(9-(e+4))] = -1; // positive mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(e+4)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(9-(e+4))] = 1; // negative mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(e+4)] = 1; // positive reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(9-(e+4))] = -1; //  positive reeomb has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(e+4)] = -1; // negative reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(9-(e+4))] = 1; // negative reeomb has opposite effeet to the opposite genotype

        			r.reactStat[(x1*2+x2*2+x3*e+6)][(e+4)] = -1; // lysis
        		}
            break;


        case 13: //XaXrXnr

            for (a=0; a < 36; a++){
                for (b=0; b < 6; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 4;
            x2 = 7;
            x3 = 7;

            for (c = 0; c < 2; c++) {
                	r.reactStat[(x1*c)][c] = 1 ; // growth
                	r.reactStat[(x1*c+1)][c] = -1; // death
                	r.reactStat[(x1*c+2)][c] = 1; // positive mutation
                	r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
                	r.reactStat[(x1*c+3)][c] = -1; // negative mutation
                	r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
            }


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}


        		for (e = 0; e < 2; e++) {
        			r.reactStat[(x1*2+x2*2+x3*e)][(e+4)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*2+x3*e+1)][(e+4)] = -1; // death
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(e+4)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(9-(e+4))] = -1; // positive mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(e+4)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(9-(e+4))] = 1; // negative mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(e+4)] = 1; // positive reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(9-(e+4))] = -1; //  positive reeomb has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(e+4)] = -1; // negative reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(9-(e+4))] = 1; // negative reeomb has opposite effeet to the opposite genotype

        			r.reactStat[(x1*2+x2*2+x3*e+6)][(e+4)] = -1; // lysis
        		}
            break;



        case 14: //XnXrXnr

            for (a=0; a < 38; a++){
                for (b=0; b < 6; b++){
                     r.reactStat[a][b]=0;
                }
            }

            x1 = 5;
            x2 = 7;
            x3 = 7;

            for (c = 0; c < 2; c++) {
        			r.reactStat[(x1*c)][c] = 1 ; // growth2
        			r.reactStat[(x1*c+1)][c] = -1; // ceath
        			r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        			r.reactStat[(x1*c+2)][1-c] = -1; // positive mutation has opposite effect to the opposite genotype
        			r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        			r.reactStat[(x1*c+3)][1-c] = 1; // negative mutation has opposite effect to the opposite genotype

        			r.reactStat[(x1*c+5)][c] = -1; // lysis
        		}


            for (d = 0; d < 2; d++) {
        			r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth
        			r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
        			r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype
        			r.reactStat[(x1*2+x2*d+4)][(d+2)] = 1; // positive recomb
        			r.reactStat[(x1*2+x2*d+4)][(5-(d+2))] = -1;  //  positive recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+5)][(d+2)] = -1; // negative recomb
        			r.reactStat[(x1*2+x2*d+5)][(5-(d+2))] = 1; // negative recomb has opposite effect to the opposite genotype
        			r.reactStat[(x1*2+x2*d+6)][(d+2)] = -1; // lysis
        		}


        		for (e = 0; e < 2; e++) {
        			r.reactStat[(x1*2+x2*2+x3*e)][(e+4)] = 1 ; // growth2
        			r.reactStat[(x1*2+x2*2+x3*e+1)][(e+4)] = -1; // death
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(e+4)] = 1; // positive mutation
        			r.reactStat[(x1*2+x2*2+x3*e+2)][(9-(e+4))] = -1; // positive mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(e+4)] = -1; // negative mutation
        			r.reactStat[(x1*2+x2*2+x3*e+3)][(9-(e+4))] = 1; // negative mutation has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(e+4)] = 1; // positive reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+4)][(9-(e+4))] = -1; //  positive reeomb has opposite effeet to the opposite genotype
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(e+4)] = -1; // negative reeomb
        			r.reactStat[(x1*2+x2*2+x3*e+5)][(9-(e+4))] = 1; // negative reeomb has opposite effeet to the opposite genotype

        			r.reactStat[(x1*2+x2*2+x3*e+6)][(e+4)] = -1; // lysis
        		}
            break;


        case 15: //XaXnXrXnr
        		for (a=0; a < 46; a++){
        	        for (b=0; b < 8; b++){
        	            r.reactStat[a][b]=0;
        	        }
        	    }

        	    x1 = 4;
        	    x2 = 5;
        	    x3 = 7;
        	    x4 = 7;

        	    for (c = 0; c < 2; c++) {
        	        r.reactStat[(x1*c)][c] = 1 ; // growth
        	        r.reactStat[(x1*c+1)][c] = -1; // death
        	        r.reactStat[(x1*c+2)][c] = 1; // positive mutation
        	        r.reactStat[(x1*c+2)][(1-c)] = -1; // positive mutation has opposite effect to the opposite genotype
        	        r.reactStat[(x1*c+3)][c] = -1; // negative mutation
        	        r.reactStat[(x1*c+3)][(1-c)] = 1; // negative mutation has opposite effect to the opposite genotype
        	    }


            for (d = 0; d < 2; d++) {
            		r.reactStat[(x1*2+x2*d)][(d+2)] = 1 ; // growth2
            		r.reactStat[(x1*2+x2*d+1)][(d+2)] = -1; // death
            		r.reactStat[(x1*2+x2*d+2)][(d+2)] = 1; // positive mutation
            		r.reactStat[(x1*2+x2*d+2)][(5-(d+2))] = -1; // positive mutation has opposite effedt to the opposite genotype
            		r.reactStat[(x1*2+x2*d+3)][(d+2)] = -1; // negative mutation
            		r.reactStat[(x1*2+x2*d+3)][(5-(d+2))] = 1; // negative mutation has opposite effedt to the opposite genotype

            		r.reactStat[(x1*2+x2*d+4)][(d+2)] = -1; // lysis
            }


            for (e = 0; e < 2; e++) {
            		r.reactStat[(x1*2+x2*2+x3*e)][(e+4)] = 1 ; // growth
            		r.reactStat[(x1*2+x2*2+x3*e+1)][(e+4)] = -1; // death
            		r.reactStat[(x1*2+x2*2+x3*e+2)][(e+4)] = 1; // positive mutation
            		r.reactStat[(x1*2+x2*2+x3*e+2)][(9-(e+4))] = -1; // positive mutation has opposite effeet to the opposite genotype
            		r.reactStat[(x1*2+x2*2+x3*e+3)][(e+4)] = -1; // negative mutation
            		r.reactStat[(x1*2+x2*2+x3*e+3)][(9-(e+4))] = 1; // negative mutation has opposite effeet to the opposite genotype
            		r.reactStat[(x1*2+x2*2+x3*e+4)][(e+4)] = 1; // positive reeomb
            		r.reactStat[(x1*2+x2*2+x3*e+4)][(9-(e+4))] = -1; //  positive reeomb has opposite effeet to the opposite genotype
            		r.reactStat[(x1*2+x2*2+x3*e+5)][(e+4)] = -1; // negative reeomb
            		r.reactStat[(x1*2+x2*2+x3*e+5)][(9-(e+4))] = 1; // negative reeomb has opposite effeet to the opposite genotype
            		r.reactStat[(x1*2+x2*2+x3*e+6)][(e+4)] = -1; // lysis
            }


    			for (f = 0; f < 2; f++) {
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f)][(f+6)] = 1 ; // growth2
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+1)][(f+6)] = -1; // dfath
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+2)][(f+6)] = 1; // positivf mutation
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+2)][(13-(f+6))] = -1; // positivf mutation has oppositf ffffft to thf oppositf gfnotypf
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+3)][(f+6)] = -1; // nfgativf mutation
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+3)][(13-(f+6))] = 1; // nfgativf mutation has oppositf ffffft to thf oppositf gfnotypf
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+4)][(f+6)] = 1; // positivf rffomb
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+4)][(13-(f+6))] = -1; //  positivf rffomb has oppositf ffffft to thf oppositf gfnotypf
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+5)][(f+6)] = -1; // nfgativf rffomb
    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+5)][(13-(f+6))] = 1; // nfgativf rffomb has oppositf ffffft to thf oppositf gfnotypf

    				r.reactStat[(x1*2+x2*2+x3*2+x4*f+6)][(f+6)] = -1; // lysis
    			}

            break;




        default:
            printf("reaction system argument is missing\n");
            break;
    }

    return r;

}
