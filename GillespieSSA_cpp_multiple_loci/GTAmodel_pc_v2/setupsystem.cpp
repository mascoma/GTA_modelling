//
//  setupsystem.cpp
//  GTAcppmodelling
//
//  Created by Xin Chen on 6/22/15.
//  Copyright (c) 2015 Xin Chen. All rights reserved.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <iostream>
#include <sstream>
#include "GTAcpp_G.h"

using namespace std; 

statusMatrix getArray(){
    
    int a=0;
    int b=0;
    int c=0;
    int d=0;
    int colNum=0;
    int rowNum=0;
    int y;
    statusMatrix r;
    for (a=0; a <(numOfgeonReact+numofgeneReact); a++){
        for (b=0; b < (genotypeNum*NUM_OF_POPTYPE+genetypeNum); b++){
             r.reactStat[a][b]=0;
        }
    }
    int reactionsysIndex=0;
    reactionsysIndex = reactionSystem(REACTIONSYSTEM);
    
    switch (reactionsysIndex) {
        case 1:
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
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
        case 2:
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+4)][c] = -1; // lysis
            }
            // status changes of gene fragments influenced by lysis
            
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                    y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES)+4)][genotypeNum +y] = 1;
                }
            }
            
            for(d=genotypeNum; d < (genotypeNum+genetypeNum); d++){
                r.reactStat[((genotypeNum*NUM_OF_PROCESSES)+(d-genotypeNum))][d]=-1;// gene fragement taken up
            }
            
            break;
            
        case 3:
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+4)][c] = -1; // lysis
                r.reactStat[((c*NUM_OF_PROCESSES)+5)][c] = 1; // positive recomb
                r.reactStat[((c*NUM_OF_PROCESSES)+6)][c] = -1; // negative recomb
            }
            // status changes of gene fragments influenced by lysis
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                    y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES)+4)][genotypeNum +y] = 1;
                }
            }
            
            for(d=(genotypeNum*NUM_OF_POPTYPE); d < (genotypeNum+genetypeNum); d++){
                r.reactStat[((genotypeNum*NUM_OF_PROCESSES)+(d-genotypeNum))][d]=-1;// gene fragement taken up
            }
            
            break;
           
        case 4: // XaXn
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
            }
            
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)][(c+genotypeNum)] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+1)][(c+genotypeNum)] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+2)][(c+genotypeNum)] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+3)][(c+genotypeNum)] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][(c+genotypeNum)] = -1; // lysis
            }
            
            // status changes of gene fragments influenced by lysis
            
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                    y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][genotypeNum*2 +y] = 1;
                }
            }
            
            for(d=genotypeNum*2; d < (genotypeNum*2+genetypeNum); d++){
                r.reactStat[((genotypeNum*NUM_OF_PROCESSES*2)+(d-genotypeNum*2))][d]=-1;// gene fragement taken up
            }

            
            break;
    
        case 5: //XaXr
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
            }
            
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)][c+genotypeNum] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+1)][c+genotypeNum] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+2)][c+genotypeNum] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+3)][c+genotypeNum] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][c+genotypeNum] = -1; // lysis
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+5)][c+genotypeNum] = 1; // positive recomb
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+6)][c+genotypeNum] = -1; // negative recomb
            }

            // status changes of gene fragments influenced by lysis
            
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                    y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][genotypeNum*2 +y] = 1;
                }
            }
            
            for(d=genotypeNum*2; d < (genotypeNum*2+genetypeNum); d++){
                r.reactStat[((genotypeNum*NUM_OF_PROCESSES*2)+(d-genotypeNum*2))][d]=-1;// gene fragement taken up
            }

            break;
            
        case 6: //XnXr
            
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+4)][c] = -1; // lysis
            }
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                		y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                		r.reactStat[((colNum*NUM_OF_PROCESSES)+4)][genotypeNum*2 +y] = 1;
                }
            }
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)][c+genotypeNum] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+1)][c+genotypeNum] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+2)][c+genotypeNum] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+3)][c+genotypeNum] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][c+genotypeNum] = -1; // lysis
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+5)][c+genotypeNum] = 1; // positive recomb
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+6)][c+genotypeNum] = -1; // negative recomb
            }
            
            // status changes of gene fragments influenced by lysis
            
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                    y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][genotypeNum*2 +y] = 1;
                }
            }
            
            for(d=genotypeNum*2; d < (genotypeNum*2+genetypeNum); d++){
                r.reactStat[((genotypeNum*NUM_OF_PROCESSES*2)+(d-genotypeNum*2))][d]=-1;// gene fragement taken up
            }

            break;
            
        case 7: //XaXnXr
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES)][c] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES)+1)][c] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES)+2)][c] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES)+3)][c] = -1; // negative mutation
            }
            
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)][(c+genotypeNum)] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+1)][(c+genotypeNum)] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+2)][(c+genotypeNum)] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+3)][(c+genotypeNum)] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][(c+genotypeNum)] = -1; // lysis
            }
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                		y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES)+4)][genotypeNum*3 +y] = 1;
                    }
            }
            for (c=0; c < genotypeNum; c++){
                r.reactStat[(c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)][c+genotypeNum*2] = 1; // growth
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+1)][c+genotypeNum*2] = -1; // death
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+2)][c+genotypeNum*2] = 1; // positive mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+3)][c+genotypeNum*2] = -1; // negative mutation
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+4)][c+genotypeNum*2] = -1; // lysis
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+5)][c+genotypeNum*2] = 1; // positive recomb
                r.reactStat[((c*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+6)][c+genotypeNum*2] = -1; // negative recomb
            }
            
            // status changes of gene fragments influenced by lysis
            
            for(colNum = 0; colNum < genotypeNum; colNum++){
                for(rowNum=0; rowNum < NUM_LOCI; rowNum++){
                    y = matrixIndex(rowNum, colNum)-1;
                    //printf("Index number [%d][%d] is %d\n",colNum, rowNum, y);
                    r.reactStat[((colNum*NUM_OF_PROCESSES+genotypeNum*NUM_OF_PROCESSES*2)+4)][genotypeNum*3 +y] = 1;
                }
            }
            
            for(d=genotypeNum*3; d < (genotypeNum*3+genetypeNum); d++){
                r.reactStat[((genotypeNum*NUM_OF_PROCESSES*3)+(d-genotypeNum*3))][d]=-1;// gene fragement taken up
            }

            
            break;
        
        default:
            printf("reaction system argument is missing\n");
            break;
    }

    return r;
}

string headerName(string s){
	string header;
	string temp1;
	string header1;
	string temp2;
	string header2;
	string temp3;
	string header3;
	string temp4;
	string header4;
	long long j=0;
	long long k=0;
	long long l=0;
	long long m=0;
	int headerIndex;
	headerIndex = reactionSystem(REACTIONSYSTEM.c_str());
	switch (headerIndex) {
	        case 1:
	            for (j=0; j < genotypeNum; j++){
	            		std::string a = std::to_string(j);
	            		temp1 = "Xa" + a + "\t";
	            		header1 = header1 + temp1;
	            }
	            for (k=0; k < genetypeNum; k++){
	            		std::string b = std::to_string(k);
	            		temp2 = "f" + b + "\t";
	            		header2 = header2 + temp2;
	            }
	            header = header1 + header2;
	            break;

	        case 2:
	        		for (j=0; j < genotypeNum; j++){
	        		    std::string a = std::to_string(j);
	        		    temp1 = "Xn" + a + "\t";
	        		    header1 = header1 + temp1;
	        		}
	        		for (k=0; k < genetypeNum; k++){
	        		    std::string b = std::to_string(k);
	        		    temp2 = "f" + b + "\t";
	        		    header2 = header2 + temp2;
	        		}
	        		header = header1 + header2;
	        		break;


	        case 3:
	        		for (j=0; j < genotypeNum; j++){
	        		     std::string a = std::to_string(j);
	        		     temp1 = "Xr" + a + "\t";
	        		     header1 = header1 + temp1;
	        		}
	        		for (k=0; k < genetypeNum; k++){
	        		     std::string b = std::to_string(k);
	        		     temp2 = "f" + b + "\t";
	        		     header2 = header2 + temp2;
	        		}
	        		header = header1 + header2;
	        		break;

	        case 4: // XaXn
	        		for (j=0; j < genotypeNum; j++){
	        			std::string a = std::to_string(j);
	        		    temp1 = "Xa" + a + "\t";
	        		    header1 = header1 + temp1;
	        		}
	        		for (k=0; k < genotypeNum; k++){
	        			std::string b = std::to_string(k);
	        		    temp2 = "Xn" + b + "\t";
	        		    header2 = header2 + temp2;
	        		}
	        		for (l=0; l < genetypeNum; l++){
	        			std::string c = std::to_string(l);
	        			temp3 = "f" + c + "\t";
	        			header3 = header3 + temp3;
	        		}
	        		header = header1 + header2 + header3;
	        		break;



	        case 5: //XaXr
	     		for (j=0; j < genotypeNum; j++){
	     			std::string a = std::to_string(j);
	    	        		temp1 = "Xa" + a + "\t";
	    	        		header1 = header1 + temp1;
	    	        }
	    	        for (k=0; k < genotypeNum; k++){
	    	        		std::string b = std::to_string(k);
	    	        		temp2 = "Xr" + b + "\t";
	    	        		header2 = header2 + temp2;
	    	        }
	    	        for (l=0; l < genetypeNum; l++){
	    	        		std::string c = std::to_string(l);
	    	        		temp3 = "f" + c + "\t";
	    	        		header3 = header3 + temp3;
	    	        	}
	    	        header = header1 + header2 + header3;
	    	        	break;


	        case 6: //XnXr

	        		for (j=0; j < genotypeNum; j++){
	        			std::string a = std::to_string(j);
	        		    	temp1 = "Xn" + a + "\t";
	        		    	header1 = header1 + temp1;
	        		}
	        		for (k=0; k < genotypeNum; k++){
	        			std::string b = std::to_string(k);
	        		    	temp2 = "Xr" + b + "\t";
	        		    	header2 = header2 + temp2;
	        		}
	        		for (l=0; l < genetypeNum; l++){
	        			std::string c = std::to_string(l);
	        		    	temp3 = "f" + c + "\t";
	        		    	header3 = header3 + temp3;
	        		}
	        		header = header1 + header2 + header3;
	        		break;

	        case 7: //XaXnXr
	        		for (j=0; j < genotypeNum; j++){
	        			std::string a = std::to_string(j);
	        		   	temp1 = "Xa" + a + "\t";
	        		    header1 = header1 + temp1;
	        		}
	        		for (k=0; k < genotypeNum; k++){
	        			std::string b = std::to_string(k);
	        		    temp2 = "Xn" + b + "\t";
	        		   	header2 = header2 + temp2;
	        		}
	        		for (l=0; l < genotypeNum; l++){
	        			std::string c = std::to_string(l);
	        		    	temp3 = "Xr" + c + "\t";
	        		    	header3 = header3 + temp3;
	        		}
	        		for (m=0; m < genetypeNum; m++){
	        			std::string c = std::to_string(m);
	        			temp4 = "f" + c + "\t";
	        			header4 = header4 + temp4;
	        		}
	        		header = header1 + header2 + header3 + header4;
	            break;

	        default:
	            printf("header error!!");
	            break;
	    }


	return header;
}
