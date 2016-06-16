### GTA_modelling

This program (GillespieSSA_cpp_V3) is for Gillespie stochastic simulation to model the dynamic of host bacterial populations under the influence of GTA. 
* The selfish model is to test the gene level selection hypothesis without considering the benefits from GTA activities. 
* The recombination model is to test the group selection hypothesis considering that the GTA mediated HGT as nutritioanl supply and/or enhancing diversity. 


Xa: bacteria population without GTA and won't take GTA
Xn: bacteria population take GTA as "food"
Xr: bacteria population integrate infected DNA fragments into their own genome
Xnr: bacteria population take GTA as "food" and integrate infected DNA fragments into their own genome
 
 
**compile:** 

* PC version:
```bash
g++ -o out *.cpp 
```
* Cluster version:

using SPRNG for random number generator 

```bash
module load intel-compilers/15.0
```
 
```bash
icc -O3 -mkl -I/opt/sprng/4.4/include -std=c++11 -o model_XaXnXr *.cpp /opt/sprng/4.4/lib/libsprng.a
```

**run program:** 

modifying the parameters in "input.txt" to define the reaction system and processes
```bash
./out input.txt
```

