AMOVA style.

Getting a list of individuals formatted:
`
# First, we do a list of the individuals:
bcftools query -l amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.IndsWith0.5MissingnessRemoved.IndsFromHAremoved.MigrantsRemoved.vcf > ind.list

# Then we make a list for the populations:
cat ind.list | awk -F "_" '{print $(NF-1)"_"$(NF)}' > pop.list

# and homogeneize names
sed -i "s/C_LV/C_Ni/g; s/C_Om/C_Ni/g; s/C_LV/C_Ni/g;s/C_PD/C_Ni/g; s/C_So/C_Ni/g;  s/L_LV/L_Ni/g; s/L_Om/L_Ni/g; s/L_LV/L_Ni/g;s/L_PD/L_Ni/g; s/L_So/L_Ni/g;" pop.list

# and make a file that reads: <ind><tab><population> by pasting:
paste ind.list pop.list > popmap.tsv

# Now, using nano, we add, so that Arlequin works:
#Individual<tab>Population
nano popmap.tsv

`

Now, converting the vcf using PGDspider:

`
ml PGDSpider/2.1.1.5-Java-1.8

#let's generate a template.. it'll give us an error!
java -Xmx1024m -Xms512m -jar $EBROOTPGDSPIDER/PGDSpider2-cli.jar -inputfile amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf.IndsWith0.5MissingnessRemoved.IndsFromHAremoved.MigrantsRemoved.vcf -inputformat VCF -outputfile final.arq -outputformat ARLEQUIN
`
This code generates a template but gives out a error. In template, and using nano, I changed:

`
#               VCF_PARSER_PLOIDY_QUESTION=DIPLOID
#               VCF_PARSER_MONOMORPHIC_QUESTION=false
#               VCF_PARSER_POP_FILE_QUESTION=./popmap.tsv
#               VCF_PARSER_POP_QUESTION=true
`

Ready to convert!!
`
java -Xmx1024m -Xms512m -jar $EBROOTPGDSPIDER/PGDSpider2-cli.jar -inputfile amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf.IndsWith0.5MissingnessRemoved.IndsFromHAremoved.MigrantsRemoved.vcf -inputformat VCF -outputfile final.arq -outputformat ARLEQUIN -spid template_VCF_ARLEQUIN.spid
`

Now, on Arquelin, we point and click:
1. Arquelin v3.5.2.2
  1.1 Click on "Open project" > select file created above (i.e. "final.FasSimCoal.onlyHAremoved.arq")
  1.2 Click on "Structure editor",  define groups of populations based on groups. Basically, all different populations are part of "group 1". Click "update project" down below
  1.3 Click on "Settings", & select
    1.3.1 AMOVA > "Standard AMOVA computations", "Locus by locus amova"
    1.3.2 Population comparisons > "Pairwise FST"
