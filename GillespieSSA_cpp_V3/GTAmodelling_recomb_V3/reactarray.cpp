/*reactarray.cpp
Jun 3, 2016
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

int *reactArray(int *popSize, int *totalPopSize, double lysis, double fitness, long *numOfGTA1, long *numOfGTA0, long *totalGTA, string s){

	int systemIndex = 0;
	int a = 0;
	int b;
	double c = 0;
	int v = 0;
	int w = 0;
	int x = 0;
	int y = 0;
	int z = 0;

	static int reactions[MAX_NUM_ARRAY];
	for(b=0; b<MAX_NUM_ARRAY; b++){
		reactions[b] = 0;
	}

	systemIndex = reactionSystem(s);

	switch(systemIndex){
		case 1: //Xa

			x = 1; // num of strain
			y = 4; // num of processes

			for (a = 0; a < 2*x; a++) {
				reactions[a*y] = growth(popSize[a], totalPopSize, 0);
				if (a%2 == 0){
					c = fitness;
				}
				else {
					c = fitness-1;
				}
				reactions[(a*y+1)] =  death(c, 0, popSize[a]);
				reactions[(a*y+2)] =  mut(0, popSize[(1-a)]);
				reactions[(a*y+3)] =  mut(0, popSize[a]);
			}
			break;

		case 2: //Xn

			x = 1; // num of strain
			y = 5; // num of processes


			for (a = 0; a < 2*x; a++) {
				reactions[a*y] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
				if (a%2 == 0){
					c = fitness;
				}
				else {
					c = fitness-1;
				}
				reactions[(a*y+1)] = death(c, lysis, popSize[a]);
				reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
				reactions[(a*y+3)] = mut(lysis, popSize[a]);
				reactions[(a*y+4)] = lysisDeath(popSize[a]);
			}
			break;

		case 3: //Xr

			x = 1; // num of strain
			y = 7; // num of processes

			for (a = 0; a < 2*x; a++) {
				reactions[a*y] = growth(popSize[a], totalPopSize, lysis);
				if (a%2 == 0){
					c = fitness;
					reactions[(a*y+4)] = recomb(numOfGTA1, popSize[(1-a)]);
					reactions[(a*y+5)] = recomb(numOfGTA0, popSize[a]);
				}
				else {
					c = fitness-1;
					reactions[(a*y+4)] = recomb(numOfGTA0, popSize[(1-a)]);
					reactions[(a*y+5)] = recomb(numOfGTA1, popSize[a]);

				}
				reactions[(a*y+1)] = death(c, lysis, popSize[a]);
				reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
				reactions[(a*y+3)] = mut(lysis, popSize[a]);
				reactions[(a*y+6)] = lysisDeath(popSize[a]);
			}
			break;

		case 4: //Xnr

			x = 1; // num of strain
			y = 7; // num of processes

			for (a = 0; a < 2*x; a++) {
				reactions[a*y] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
				if (a%2 == 0){
					c = fitness;
					reactions[(a*y+4)] = recomb(numOfGTA1, popSize[(1-a)]);
					reactions[(a*y+5)] = recomb(numOfGTA0, popSize[a]);
				}
				else {
					c = fitness-1;
					reactions[(a*y+4)] = recomb(numOfGTA0, popSize[(1-a)]);
					reactions[(a*y+5)] = recomb(numOfGTA1, popSize[a]);
				}
				reactions[(a*y+1)] = death(c, lysis, popSize[a]);
				reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
				reactions[(a*y+3)] = mut(lysis, popSize[a]);
				reactions[(a*y+6)] = lysisDeath(popSize[a]);

			}
			break;

		case 5: //XaXn

			x = 2; // num of strain
			y = 4; // num of processes of strain 1
			z = 5; // num of processes of strain 2

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+4)] = lysisDeath(popSize[a]);
				}
			}
			break;

		case 6: //XaXr

			x = 2; // num of strain
			y = 4; // num of processes of strain 1
			z = 7; // num of processes of strain 2

			for (a = 0; a < 2*x; a++) {

				if (a < 2){
					reactions[a*y] = growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;


		case 7: //XaXnr

			x = 2; // num of strain
			y = 4; // num of processes of strain 1
			z = 7; // num of processes of strain 2

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y]= growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;


		case 8: //XnXr
			x = 2; // num of strain
			y = 5; // num of processes strain1
			z = 7; // num of processes strain2

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, lysis, popSize[a]);
					reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(lysis, popSize[a]);
					reactions[(a*y+4)] = lysisDeath(popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;

		case 9: //XnXnr

			x = 2; // num of strain
			y = 5; // num of processes strain1
			z = 7; // num of processes strain2

			for (a = 0; a < 2*x; a++) {

				if (a < 2){
					reactions[a*y] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, lysis, popSize[a]);
					reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(lysis, popSize[a]);
					reactions[(a*y+4)] = lysisDeath(popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;


		case 10: //XrXnr
			x = 2; // num of strain
			y = 7; // num of processes strain1
			z = 7; // num of processes strain2

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[(a*y)] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(a*y+4)] = recomb(numOfGTA1, popSize[(1-a)]);
						reactions[(a*y+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(a*y+4)] = recomb(numOfGTA0, popSize[(1-a)]);
						reactions[(a*y+5)] = recomb(numOfGTA1, popSize[a]);

					}
					reactions[(a*y+1)] = death(c, lysis, popSize[a]);
					reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(lysis, popSize[a]);
					reactions[(a*y+6)] = lysisDeath(popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(a*y)] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(a*y+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(a*y+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(a*y+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(a*y+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(a*y+1)] = death(c, lysis, popSize[a]);
					reactions[(a*y+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(a*y+3)] = mut(lysis, popSize[a]);
					reactions[(a*y+6)] = lysisDeath(popSize[a]);
				}

			}
			break;

		case 11: //XaXnXr

			x = 3; // num of strain
			y = 4; // num of processes of strain 1
			z = 5; // num of processes of strain 2
			w = 7;

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] =growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(2*y+z*(a-2)+1)] =death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+4)] = lysisDeath(popSize[a]);
				}
				else if(3<a && a<6){
					reactions[(2*y+z*2+w*(a-4))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA1, popSize[(9-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA0, popSize[(9-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*2+w*(a-4)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+2)] = mut(lysis, popSize[(9-a)]);
					reactions[(2*y+z*2+w*(a-4)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;

		case 12: //XaXnXnr

			x = 3; // num of strain
			y = 4; // num of processes of strain 1
			z = 5; // num of processes of strain 2
			w = 7;

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+4)] = lysisDeath(popSize[a]);
				}
				else if (3<a && a<6){
					reactions[(2*y+z*2+w*(a-4))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA1, popSize[(9-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA0, popSize[(9-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*2+w*(a-4)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+2)] = mut(lysis, popSize[(9-a)]);
					reactions[(2*y+z*2+w*(a-4)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;


		case 13: //XaXrXnr

			x = 3; // num of strain
			y = 4; // num of processes of strain 1
			z = 7; // num of processes of strain 2
			w = 7;

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
				else if (3<a && a<6){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(9-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(9-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(9-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);

				}
 			}
			break;

		case 14: //XnXrXnr

			x = 3; // num of strain
			y = 5; // num of processes of strain 1
			z = 7; // num of processes of strain 2
			w = 7;

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, lysis, popSize[a]);
					reactions[(a*y+2)] = mut(lysis, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(lysis, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(5-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
				else if (3<a && a<6){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA1, popSize[(9-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*(a-2)+4)] = recomb(numOfGTA0, popSize[(9-a)]);
						reactions[(2*y+z*(a-2)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(9-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;

		case 15: //XaXnXrXnr

			x = 4; // num of strain
			y = 4; // num of processes of strain 1
			z = 5; // num of processes of strain 2
			w = 7;
			v = 7;

			for (a = 0; a < 2*x; a++) {
				if (a<2){
					reactions[a*y] = growth(popSize[a], totalPopSize, 0);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(a*y+1)] = death(c, 0, popSize[a]);
					reactions[(a*y+2)] = mut(0, popSize[(1-a)]);
					reactions[(a*y+3)] = mut(0, popSize[a]);
				}
				else if (1<a && a<4){
					reactions[(2*y+z*(a-2))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
					}
					else {
						c = fitness-1;
					}
					reactions[(2*y+z*(a-2)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+2)] = mut(lysis, popSize[(5-a)]);
					reactions[(2*y+z*(a-2)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*(a-2)+4)] = lysisDeath(popSize[a]);
				}
				else if (3<a && a<6){
					reactions[(2*y+z*2+w*(a-4))] = growth(popSize[a], totalPopSize, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA1, popSize[(9-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA0, popSize[(9-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*2+w*(a-4)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+2)] = mut(lysis, popSize[(9-a)]);
					reactions[(2*y+z*2+w*(a-4)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+6)] = lysisDeath(popSize[a]);
				}
				else if (5<a && a<8){
					reactions[(2*y+z*2+w*(a-4))] = growth2(popSize[a], totalPopSize, totalGTA, lysis);
					if (a%2 == 0){
						c = fitness;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA1, popSize[(13-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA0, popSize[a]);
					}
					else {
						c = fitness-1;
						reactions[(2*y+z*2+w*(a-4)+4)] = recomb(numOfGTA0, popSize[(13-a)]);
						reactions[(2*y+z*2+w*(a-4)+5)] = recomb(numOfGTA1, popSize[a]);
					}
					reactions[(2*y+z*2+w*(a-4)+1)] = death(c, lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+2)] = mut(lysis, popSize[(13-a)]);
					reactions[(2*y+z*2+w*(a-4)+3)] = mut(lysis, popSize[a]);
					reactions[(2*y+z*2+w*(a-4)+6)] = lysisDeath(popSize[a]);
				}
			}
			break;

        default:
            printf("reaction system argument is missing\n");

			break;


	}
	return reactions;

}
