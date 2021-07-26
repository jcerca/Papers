We want a VCF without LD and MAF cut-offs.
However, this VCF had too many SNPs, and dividing it in populations (code silenced just below) did not work well. I'm working with 9,000 snps.

`
# Getting 9,000 snps
grep -v "^#" final.FasSimCoal.onlyHAremoved.vcf | shuf | head -n 9000 > 9000.snps

# Separating the VCF
grep "^#" final.FasSimCoal.onlyHAremoved.vcf > header.tsv
cat header.tsv 9000.snps > filtered.vcf

# Getting a list of individuals
ml BCFtools/1.9-intel-2018b
bcftools query -l filtered.vcf > ind.list

# Separing the species onto different files:
grep "C_" ind.list > all_citrinellus.tsv
grep "L_" ind.list > all_labiatus.tsv

# Getting a vcf for each species:
vcftools --vcf filtered.vcf --keep all_citrinellus.tsv --out onlyCitrinellus.vcf --recode --recode-INFO-all
vcftools --vcf filtered.vcf --keep all_labiatus.tsv --out onlyLabiatus.vcf  --recode --recode-INFO-all

# Getting a population list
bcftools query -l onlyCitrinellus.vcf.recode.vcf > Citrinellus.ind.list
# Now, because we want to understand NE by  lake, we will convert C_LV, C_Om, C_PD, C_So to C_Ni
bcftools query -l onlyCitrinellus.vcf.recode.vcf | awk -F "_" '{print $(NF-1)"_"$(NF)}' | sed "s/C_LV/C_Ni/; s/C_Om/C_Ni/; s/C_PD/C_Ni/; s/C_So/C_Ni/" > Citrinellus.pop.list
bcftools query -l onlyLabiatus.vcf.recode.vcf > Labiatus.ind.list
bcftools query -l onlyLabiatus.vcf.recode.vcf | awk -F "_" '{print $(NF-1)"_"$(NF)}' | sed "s/L_LV/L_Ni/; s/L_Om/L_Ni/; s/L_PD/L_Ni/; s/L_So/L_Ni/" > Labiatus.pop.list

paste Citrinellus.ind.list Citrinellus.pop.list > Citrinellus.popmap.tsv
nano Citrinellus.popmap.tsv
        #as first row, one must add [so NeEstimator runs]:
        #Individual<tab>Population

paste Labiatus.ind.list Labiatus.pop.list > Labiatus.popmap.tsv
nano Labiatus.popmap.tsv
        #as first row, add:
        #Individual<tab>Population
        
`
Now, we will be using PGDSpider to convert the file.

`
ml PGDSpider/2.1.1.5-Java-1.8
java -Xmx1024m -Xms512m -jar $EBROOTPGDSPIDER/PGDSpider2-cli.jar -inputfile filtered.vcf -inputformat VCF -outputfile final.FasSimCoal.onlyHAremoved.genepop -outputformat GENEPOP

`
This gives us the following error:

ERROR, it needs a spid file. However, it also created a .spid file.
I saw how to correct the file, here: https://github.com/jpuritz/BIO_594_2018/blob/master/FinalAssignment/Ian%20Bishop/BSsnp.spid
So, using nano we add:

`
nano template_VCF_GENEPOP.spid
#       added:
#               VCF_PARSER_PLOIDY_QUESTION=DIPLOID
#               VCF_PARSER_MONOMORPHIC_QUESTION=false
#               VCF_PARSER_POP_FILE_QUESTION=./Citrinellus.popmap.tsv
#               VCF_PARSER_POP_FILE_QUESTION=./Labiatus.popmap.tsv
#               VCF_PARSER_POP_QUESTION=true
`

I made one for each species.

`
java -Xmx1024m -Xms512m -jar $EBROOTPGDSPIDER/PGDSpider2-cli.jar -inputfile onlyCitrinellus.vcf.recode.vcf -inputformat VCF -outputfile Citrinellus.genepop -outputformat GENEPOP -spid template_CITRI_VCF_GENEPOP.spid
java -Xmx1024m -Xms512m -jar $EBROOTPGDSPIDER/PGDSpider2-cli.jar -inputfile onlyLabiatus.vcf.recode.vcf -inputformat VCF -outputfile Labiatus.genepop -outputformat GENEPOP -spid template_LABI_VCF_GENEPOP.spid
`

Done!
