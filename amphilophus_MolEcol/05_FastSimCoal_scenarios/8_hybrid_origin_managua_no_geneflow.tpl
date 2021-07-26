//Parameters for the coalescence simulation program : fsimcoal2.exe
4 samples to simulate :
//Population effective sizes (number of genes)
Nnc
Nnl
Nmc
Nml
//Samples sizes and samples age
30
30
30
30
//Growth rates	: negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
4  historical events
Th 3 1 ALA RS1 0 0
Th 3 2 1 RS2 0 0
Tn 1 0 1 RSm 0 0
Ta 2 0 1 RSa 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   1.2e-8 OUTEXP
