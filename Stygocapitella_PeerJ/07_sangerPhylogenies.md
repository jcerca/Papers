I did include some sanger phylogenies in this work. You can download the COI, 18s, 16s and ITS1 datasets from NCBI following the NCBI assession numbers on Supplementary Table 2.

If you have a fasta file with sequences from a single gene, you only need to align, and make a tree:
```
# Alignement
mafft --auto 16s.fasta > 16s.aligned.fasta
# Tree
iqtree -s 16s.aligned.fasta -bb 1000
```

If you want to concatenate several genes, so you can analyse them, here's an example:

```
# First, you align
mafft --maxiterate 1000 --localpair --reorder 16s.fas > 16s.aligned.fas
mafft --maxiterate 1000 --localpair --reorder co1.fas > co1.aligned.fas
mafft --maxiterate 1000 --globalpair --reorder its1.fas > its1.aligned.fas
mafft --maxiterate 1000 --localpair --reorder 18s.fas > 18s.aligned.fas

# Then, we concatenate. To concatenate we use FasConCat. Notice, for it to work, fasta files need to be called ".fas"

~/perl5/perlbrew/perls/perl-5.20.1/bin/perl ~/00_manually_installed_stuff/FASconCAT-master/FASconCAT_v1.1.pl -s

# We need to make a partitions.nex file
nano partitions.nex
# It looks like this (without the first #).
# #nexus
# begin sets;
#         charset 16s_part1 = 1-944;
#         charset 18s_part2 = 945-2745;
#         charset co1_part3 = 2746-3441;
#         charset its1_part4 = 3442-4501;
#  end;

# And we make the tree:
iqtree -s FcC_smatrix.fas -spp partitions.nex -bb 1000 -nt AUTO
```
