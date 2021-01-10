//Parameters for the coalescence simulation program : fsimcoal2.exe
3 samples to simulate :
//Population effective sizes (number of genes)
Np
Ng
Nb
//Samples sizes and samples age
30
30
30
//Growth rates  : negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 Ma Mb
Mc 0 Md
Me Mf 0
//Migration matrix 1
0 Ma 0
Mc 0 0
0 0 0
//Migration matrix 2
0 0 0
0 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2  historical events
Tb 2 0 1 RSb 0 1
Ta 1 0 1 RSa 0 2
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   1.2e-8 OUTEXP
