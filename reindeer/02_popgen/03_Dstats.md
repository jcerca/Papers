We used Dsuite to calculate ABBA-BABA, ...

```
ml Dsuite/20210718-GCC-10.3.0
Dsuite Dtrios reindeer_outgroup_concatenated_noMT_noScaffolds_noSexChr_qual20_mindp15_maxdp75.vcf.gz ind_list2.tsv -t tree.nwk -o reindeer_trios
```

... and f_dM!

```
Dsuite Dinvestigate -w 50,25 $vcf ind_list2.tsv test_trios.tsv
```
