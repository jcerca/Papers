This one is easy. First, we make a list of individuals and assign them to populations:

```
VCF=amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.IndsWith0.5MissingnessRemoved.vcf

bcftools query -l $VCF | grep "C_Ma" > C_Ma # getting C_Ma individuals
bcftools query -l $VCF | grep "L_Ma" > L_Ma # etc etc

bcftools query -l $VCF | grep "_L_" | grep -v "_Ma" > L_Nic
bcftools query -l $VCF | grep "_C_" | grep -v "_Ma" > C_Nic
```


Now, we use fst. One example:

```
vcftools --vcf amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.IndsWith0.5MissingnessRemoved.vcf \
--fst-window-size 50000 --fst-window-step 10000 --weir-fst-pop C_Ma \
--weir-fst-pop C_Nic \
--out Cman_VS_Cnic
```

Finally, we plot:

```
rm(list = ls())
library(tidyverse)

setwd("C:/Users/josece_adm/Desktop/AMPHILOPHUS/FST/run_with_new_data/")

# We need to index and get a *fai file and then get a genome index.
genomeindex <- read.table("./CCOE01000001-CCOE01006637.fasta.fai", sep="\t") %>%
  rename(contig=V1, contigSize=V2) %>%
  select(-V3, -V4, -V5)


genomeindex$scaffoldindex<-as.factor(unlist(as.data.frame(read.table("./scaffold_ID.txt"))))


filtered_genomeindex <- filter(genomeindex, contigSize > 500000) %>% # here we filter for scaffolds bigger than 500kb
  arrange(desc(contigSize)) %>% #here we arrange the column in a descending order (biggest scaffold > 2nd biggest > ... > shortest scaffold)
  mutate(cumulativeSize=(cumsum(contigSize))) %>% # Here we create a column called cumulativesize by cumulative summing of the contig size.
  mutate(cumSize_from0 = lag(cumulativeSize)) %>% # Here we create a column called cumSize_from0 ... which is the same as the previous code, YET, the first entry is a NA.
  replace(is.na(.),0) #We replace the NA for 0

#10 biggest scaffolds
filtered_genomeindex <- filter(genomeindex, contigSize > 4500000) %>% # here we filter for scaffolds bigger than 500kb
  arrange(desc(contigSize)) %>% #here we arrange the column in a descending order (biggest scaffold > 2nd biggest > ... > shortest scaffold)
  mutate(cumulativeSize=(cumsum(contigSize))) %>% # Here we create a column called cumulativesize by cumulative summing of the contig size.
  mutate(cumSize_from0 = lag(cumulativeSize)) %>% # Here we create a column called cumSize_from0 ... which is the same as the previous code, YET, the first entry is a NA.
  replace(is.na(.),0) #We replace the NA for 0


# read in point data
fst <- read_delim("./Cman_VS_Cnic.windowed.weir.fst", delim = '\t') %>%
  rename(contig=CHROM, fst = MEAN_FST, snp=BIN_START, bp=BIN_END)

fst <- read_delim("./Lman_VS_Lnic.windowed.weir.fst", delim = '\t') %>%
  rename(contig=CHROM, fst = MEAN_FST, snp=BIN_START, bp=BIN_END)

fst <- read_delim("./Cman_VS_Lman.windowed.weir.fst", delim = '\t') %>%
  rename(contig=CHROM, fst = MEAN_FST, snp=BIN_START, bp=BIN_END)

fst <- read_delim("./Cnic_VS_Lnic.windowed.weir.fst", delim = '\t') %>%
  rename(contig=CHROM, fst = MEAN_FST, snp=BIN_START, bp=BIN_END)

#Joining the dataset based on matches. We add columns to filtered, based on fst
#inner_join :: Return all rows from x where there are matching values in y, and all columns from x and y. If there are multiple matches between x and y, all combination of the matches are returned.
fused_dataset<-inner_join(fst, filtered_genomeindex, by="contig")


a <- ggplot(fused_dataset, aes(bp/10^6, fst), colour=contig)
a <- a + geom_point(alpha = 0.1)
a <- a + xlab("Position (Mb)") + ylab(expression(italic("F")[ST]))
a <- a + facet_grid(~scaffoldindex, scales = "free_x")
#a <- a + facet_wrap(~scaffoldindex)
a + theme_light()
a + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line())
```
