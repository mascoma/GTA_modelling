### GTA_modelling

This program is for Gillespie stochastic simulation to model the dynamic of host bacterial populations under the influence of GTA. 
* The selfish model is to test the gene level selection hypothesis without considering the benefits from GTA activities. 
* The recombination model is to test the group selection hypothesis considering that the GTA mediated HGT as nutritioanl supply and/or enhancing diversity. 

 

Xa: bacteria population without GTA and won't take GTA
Xn: bacteria population take GTA as "food"
Xr: bacteria population integrate infected DNA fragments into their own genome
Xnr: bacteria population take GTA as "food" and integrate infected DNA fragments into their own genome
 
 
PC version:
compile: g++ -o out *.cpp 
run program: ./out input.txt

modifying the parameters in "input.txt" to define the reaction system and processes

Cluster version:
#module load intel-compilers/15.0
using SPRNG for random number generator
compile: icc -O3 -mkl -I/opt/sprng/4.4/include -std=c++11 -o model_XaXnXr *.cpp /opt/sprng/4.4/lib/libsprng.a
