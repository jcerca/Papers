### Final analysis! 01 - getting the vcf and scaffold lists.
```
00_vcf
cp ../../00_vcf/02_renaming/spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf .
grep -v "^#" spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf  | awk '{print $1}' | sort -u > scaffold_list.tsv

# We make a genome coordinates file, so we can remove stuff.
awk '{print $1, $2}' ../../../00_theGenome/T_kauaiensis_ref.fasta.fai | tr " " "\t" > genome.coords

bgzip spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf
bcftools index spider.MinCoverage7.MaxCoverage30.repeatsRemoved.renamed.vcf.gz

```
### 02 - We split the vcf.

```
01_vcfSplitting
# We extract regions of 20kb along the genome.
# Loop over windows in genome file
while read chrom start end; do
  # Define output file name
  echo ${chrom}_${start}_${end}
  output_file="${chrom}_${start}_${end}.vcf"

  # Extract variants in window from VCF
  bcftools view -r "${chrom}:${start}-${end}" $input_vcf > $output_file

done < <(bedtools makewindows -g ../00_vcf/genome.coords -w 20000)

# I tried removing remove regions where there are <50 SNPs, and got 16.000 parts of the genome but tons of missing data in some loci. So I ended up going with 100 snps.


for file in *; do

  # Count non-comment lines in file
  num_lines=$(grep -v "^#" $file | wc -l)

  # If file has fewer than 100 non-comment lines, remove it
  if [[ $num_lines -lt 100 ]]; then
    rm $file
  fi

done


```

# 03 - Converting and making trees
```
for i in *; do echo $i; python ../vcf2phylip.py --nexus -i $i; done
mv *nexus ../02_vcf2phylip/
mv *phy ../02_vcf2phylip/

conda activate phyl
for i in *.nexus; do iqtree -s ${i} -B 1000 --wbt; done
```

# 04 - monophyletic greens?
```
cd ../03*
cat *treefile > allTrees.tsv

# Cleaning trees
cat allTrees.tsv  | sed "s/[0-9]*//g; s/\.//g; s/_://g; s/://g" > cleanedTrees.tsv

# all trees are different
sort cleanedTrees.tsv | uniq -c | wc -l
4275

# Got the mainTree from the major analysis
cat mainTree.tsv  | sed "s/[0-9]*//g; s/\.//g; s/_://g; s/://g" > maintree.cleaned

# MonoPhyly check
I used phykit on saga to get whether the greens or the browns were monophyletic in any of the 20k trees
conda activate phykit; for i in ScpjasG_1*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.1.tsv; done
conda activate phykit; for i in ScpjasG_2*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.2.tsv; done
conda activate phykit; for i in ScpjasG_3*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.3.tsv; done
conda activate phykit; for i in ScpjasG_4*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.4.tsv; done
conda activate phykit; for i in ScpjasG_5*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.5.tsv; done
conda activate phykit; for i in ScpjasG_6*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.6.tsv; done
conda activate phykit; for i in ScpjasG_7*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.7.tsv; done
conda activate phykit; for i in ScpjasG_8*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.8.tsv; done
conda activate phykit; for i in ScpjasG_9*treefile; do echo $i; phykit monophyly_check $i greens >> tree_greens.9.tsv; done


## Basically on tree_greens there is one monophyletic green group!!
It is this tree
(T_anuenue_046:0.0000010000,(T_anuenue_047:0.0000010000,(T_anuenue_049:0.0206833023,T_anuenue_066:0.0000010000)90:0.0000022535)60:0.0000023484,(T_anuenue_048:0.0000010000,(((((((T_anuenue_060:0.0000010000,(((((((((T_brevignatha_007:0.0000022653,(T_brevignatha_050:0.0000010000,T_brevignatha_067:0.0000010000)85:0.0000022883)97:0.0435700646,T_tantalus_039:0.0000023599)53:0.0000025316,((T_waikamoi_062:0.0000010000,((((T_waikamoi_073:0.0000010000,T_waikamoi_076:0.0000010000)46:0.0000010000,(T_waikamoi_043:0.0000010000,T_waikamoi_045:0.0000010000)36:0.0000010000)62:0.0000010000,T_waikamoi_041:0.0000010000)79:0.0000021214,T_waikamoi_044:0.0000010000)82:0.0000021422)58:0.0001126005,T_waikamoi_042:0.0199672425)65:0.0259839229)36:0.0000021039,((T_polychromata_053:0.0000010000,T_tantalus_040:0.0000010000)86:0.0000010000,T_brevignatha_059:0.0000010000)95:0.0203040953)98:0.0624248376,(((T_kauaiensis_016:0.0000010000,(T_kauaiensis_019:0.0000022994,T_kauaiensis_020:0.0000010000)59:0.0000010000)68:0.0000027800,T_kauaiensis_018:0.0000010000)87:0.0172049200,T_kauaiensis_017:0.0000029648)99:0.1578832410)55:0.0078547295,((((((T_filiciphila_074:0.0000010000,T_filiciphila_075:0.0000010000)96:0.0112817743,T_maka_023:0.0000010000)73:0.0000028010,T_maka_052:0.0000010000)59:0.0000026471,T_stelarobusta_001:0.0000010000)89:0.0086591081,T_acuta_071:0.0402063903)100:0.2214238540,(T_pilosa_028:0.0000010000,(((T_pilosa_029:0.0000010000,T_pilosa_032:0.0000010000)73:0.0000010000,T_pilosa_031:0.0000010000)64:0.0000010000,T_pilosa_030:0.0000010000)66:0.0000010000)100:0.1760814149)48:0.0037245095)31:0.0065655009,((T_mohihi_021:0.0000010000,T_mohihi_051:0.0000010000)63:0.0000010000,T_mohihi_022:0.0000010000)100:0.1068035598)79:0.0709793929,((T_restricta_037:0.0000010000,(T_restricta_069:0.0000010000,T_restricta_070:0.0000010000)77:0.0000010000)67:0.0000010000,T_restricta_038:0.0000010000)99:0.0897438974)64:0.0606388913,(((((T_kamakou_008:0.0000010000,((((T_kamakou_010:0.0000010000,T_kamakou_012:0.0000025190)39:0.0000023855,T_kamakou_013:0.0000010000)26:0.0000010000,T_kamakou_014:0.0000010000)24:0.0000010000,T_kamakou_015:0.0000010000)45:0.0000010000)37:0.0000020526,T_kamakou_011:0.0000010000)28:0.0000020405,T_kamakou_009:0.0000010000)34:0.0085664144,T_kamakou_068:0.0091981578)99:0.1603252996,((T_perreirai_024:0.0000026599,T_perreirai_025:0.0000010000)63:0.0000029445,T_perreirai_061:0.0000010000)91:0.0226507198)40:0.0294027364)100:0.1302249789)33:0.0000021032,(((T_kukuiki_026:0.0000010000,(T_quasimodo_002:0.0000010000,((T_quasimodo_004:0.0000010000,T_quasimodo_005:0.0000010000)93:0.0086366396,T_quasimodo_065:0.0000010000)45:0.0000021389)51:0.0085875909)55:0.0000027358,T_quasimodo_003:0.0000010000)47:0.0000020201,(((((((T_quasimodo_027:0.0000010000,(T_quasimodo_056:0.0000010000,T_quasimodo_033:0.0000010000)91:0.0000010000)86:0.0000021869,T_quasimodo_054:0.0000010000)81:0.0193577117,T_quasimodo_034:0.0000010000)82:0.0000029626,T_quasimodo_035:0.0000010000)79:0.0000010000,T_quasimodo_057:0.0000010000)79:0.0000023959,T_quasimodo_055:0.0000010000)93:0.0392501854,T_quasimodo_036:0.0000010000)52:0.0000028411)20:0.0000026661)30:0.0000022114,T_quasimodo_064:0.0000010000)46:0.0000027406,T_quasimodo_063:0.0000010000)75:0.0163882587,T_kikokiko_072:0.0000022016)69:0.0000024759,T_obscura_006:0.0252701758)86:0.0331430803,T_kukuiki_058:0.0000010000)64:0.0000020114)57:0.0000022091);

From the genomic region:
ScpjasG_1342_1420000_1440000.min4.nexus.treefile

# inside or nearby this region there are 5 genes
ScpjasG_1342    AUGUSTUS        gene    1416535 1416789 0.97    -       .       ID=g18403;
ScpjasG_1342    AUGUSTUS        gene    1420715 1421002 0.69    +       .       ID=g18404;
ScpjasG_1342    AUGUSTUS        gene    1425542 1431096 0.83    -       .       ID=g18405;
ScpjasG_1342    AUGUSTUS        gene    1432873 1435588 1       -       .       ID=g18406;
ScpjasG_1342    AUGUSTUS        gene    1441099 1479043 0.21    -       .       ID=g18407;

# NCBI blast says tat 18405 is this gene
# query cover 100%, e-value 1e-103, identity 99.32% to a gene in Nephila pilipes

```
## 05 - Astral
```
# 01_runWithout_uf
cat ../../03_*/*treefile > allTrees.tsv
java -jar /data/bigexpansion/jcerca/local_bin/Astral/astral.5.7.8.jar -i allTrees.tsv -o species_pp.tre  2>out.log



# 06_Astral/02_runWith_uf
ls ../../03*/*.ufboot > ml_boot.txt
cp ../06*/allTrees.tsv .

java -jar /data/bigexpansion/jcerca/local_bin/Astral/astral.5.7.8.jar -i allTrees.tsv -b ml_boot.txt -o species_boot.trees 2>out.log

# 06_Astral/03_geneConcordance
iqtree -t ../01_runWithout_uf/species_pp.tre --gcf ../01_runWithout_uf/allTrees.tsv --prefix gene_concordance
# the file gene_concordance.cf.tree has posterior probabilies from Astral, and the second one is gene concordance factors

# 06_Astral/04_SiteConcordanceFactor
ruby concatenate.rb -i ../../02_vcf2phylip/*nexus -o concatenated.nex -f nexus


iqtree -t ../01_runWithout_uf/species_pp.tre --gcf ../01_runWithout_uf/allTrees.tsv -s concatenated.nex --scf 100 --prefix gene_and_site_concordance

# gene_and_site_concordance.cf.tree has the things we need. We can visualize them as labels in Figtree
```
