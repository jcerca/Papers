First we need to convert the vcf to eigenstrat using a modified version of JMeier's script I retrieved somewhere online.

```
# Script to convert vcf to eigenstrat format for ADMIXTOOLS
# Originally written by Joana Meier - modified by José Cerca
# It takes a single argument: the vcf file (can be gzipped)
# It requires vcftools and admixtools

loaded_vcf=Stygocapitella.subterranea.josemariobrancoi.westheidei.r50.p8.stacks.maf0.05.maxMeanDP100.minMeanDP10.indswith90missingnessRemoved.randomSNP.vcf
file=${loaded_vcf%.vcf}

vcftools --vcf ${file}.vcf \
         --plink --mac 1.0 --remove-indels --max-alleles 2 --out $file



# Change the .map file to match the requirements of ADMIXTOOLS
awk -F"\t" 'BEGIN{scaff="";add=0}{
        split($2,newScaff,":")
        if(!match(newScaff[1],scaff)){
                scaff=newScaff[1]
                add=lastPos
        }
        count+=0.0001
        pos=add+$4
        print newScaff[1]"\t"$2"\t"count"\t"pos
        lastPos=pos
}' ${file}.map  | sed 's/^chr//' > better.map
mv better.map ${file}.map


# Change the .ped file to match the ADMIXTOOLS requirements
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t";
 for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' ${file}.ped > tmp.ped
mv tmp.ped ${file}.ped


# create an inputfile for convertf
echo "genotypename:    ${file}.ped" > par.PED.EIGENSTRAT.${file}
echo "snpname:         ${file}.map" >> par.PED.EIGENSTRAT.${file}
echo "indivname:       ${file}.ped" >> par.PED.EIGENSTRAT.${file}
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT.${file}
echo "genotypeoutname: ${file}.eigenstratgeno" >> par.PED.EIGENSTRAT.${file}
echo "snpoutname:      ${file}.snp" >> par.PED.EIGENSTRAT.${file}
echo "indivoutname:    ${file}.ind" >> par.PED.EIGENSTRAT.${file}
echo "familynames:     NO" >> par.PED.EIGENSTRAT.${file}

# MY CODE
# OK. The .map file needs to have "clean names" on the chromossomes.
# My chromossomes here have the following names: 87701; 58325
# we need to format it so it has "chr1" on the first row; and "chr1:<snp>" on the second row.

awk '{print $2}' ${file}.map | awk -F ":" '{print $2}' > snp.pos
# then, awk it again and get just the snp position

paste snp.pos ${file}.map | awk '{print "chr1\tchr1:"$1"\t"$4"\t"$5}' > new.map

# Now, we paste the snp.position, and the map file. ...
# THEN WE BUILD A MAP FILE WITH THE FOLLOWING STRUCTURE:
# chr1<tab>chr1:<snp position><tab><distance><distance>
# This file looks something like:
# chr1    chr1:52 0.0001  52
# chr1    chr1:225        0.0002  225
# chr1    chr1:51 0.0003  51
# chr1    chr1:130        0.0004  130
mv new.map ${file}.map

# Use CONVERTF to parse PED to eigenstrat:
convertf -p par.PED.EIGENSTRAT.${file}
```

Now we need a list of populations:

```
awk -F "_" '{print $(NF-1)}' ${file}.ind > pop
# pop is a single column file with population ID

paste ${file}.ind pop | awk '{print $1 "\t" $2 "\t" $4}' > new.${file}
mv new.${file} ${file}.ind

# change the snp file for ADMIXTOOLS:
awk 'BEGIN{i=0}{i=i+1; print $1"\t"$2"\t"$3"\t"i"\t"$5"\t"$6}' $file.snp > $file.snp.tmp
mv $file.snp.tmp $file.snp

#It needs a pop file which is the same as the .ind file so...
cp $file.ind $file.pop

#from your population list, just make a file called Dstat_quadruples with <taxon><taxon><taxon><taxon> to test
echo "st canoe glenancross reid" > Dstat_quadruples
echo "st canoe reid glenancross" >> Dstat_quadruples

echo "genotypename: $file.eigenstratgeno" > par.$file
echo "snpname:      $file.snp" >> par.$file
echo "indivname:    $file.pop" >> par.$file
echo "popfilename:  Dstat_quadruples" >> par.$file
```
Now, we get our results for the F3 test.

```
qp3Pop -p par.$file
```
