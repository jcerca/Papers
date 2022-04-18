I used Jellyfish.

Folder structure 00_genome, 01_kmerCounts, 02_dumping, 03_plotting

I ran Jellyfish separately for each chromosome. Remember to read (& cite ;) ) the paper to understand all the details.

01_kmerCounts
```
# then, we do kmer counting with jellyfish
cd ../00_genome
for i in *fa; do echo $i; jellyfish count -m 13 -s 100M -t 70 $i; mv mer_counts.jf ../02_kmer13/01_kmerCount/${i%fa.}.mer_counts.jf; done
```

02_dumping
```
cd 01*
for i in *jf; do jellyfish dump -c -L 100 $i > ${i%.jf}.dumps.larger100.col; mv *.col ../02_dumping/; done
# -L 100 means: Don't output k-mer with count < lower-count
# Adding "kmer frequency" to the top of the files"
for i in *col; do sed -i '1 i\kmer\tfrequency' $i; done
# Adding tabs so it is read-able in R
for i in *col; do sed -i "s/ /\t/" $i; done
```

03_plotting
Enter R!

```
getwd()
setwd("/Users/josepc/Desktop/tmp")

library(tidyverse)

##### PAIR 01 - we read in the individual chromosome dumps, larger than 100x.
s12<-read.table("./ScDrF4C_12.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s12=frequency)
s25<-read.table("./ScDrF4C_25.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s25=frequency)

# We join the pairs
pair01<-full_join(s12, s25, by="kmer")
head(pair01)

# filter based on the 2x data - again, read the paper
filtered_pair01<- pair01 %>% filter(s12 > 2*s25 | s25 > 2*s12)
nrow(filtered_pair01)

##### PAIR 02 - we do what we dd with pair 1 in every chr pair.
s1<-read.table("./ScDrF4C_1.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s1=frequency)
s1634<-read.table("./ScDrF4C_1634.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s1634=frequency)

pair02<-full_join(s1, s1634, by="kmer")
head(pair02)

filtered_pair02<- pair02 %>% filter(s1 > 2*s1634 | s1634 > 2*s1)
nrow(filtered_pair02)

##### PAIR 03
s116<-read.table("./ScDrF4C_116.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s116=frequency)
s10<-read.table("./ScDrF4C_10.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s10=frequency)

pair03<-full_join(s10, s116, by="kmer")
head(pair03)

filtered_pair03<- pair03 %>% filter(s10 > 2*s116 | s116 > 2*s10)
nrow(filtered_pair03)


##### PAIR 04
s11<-read.table("./ScDrF4C_11.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s11=frequency)
s30<-read.table("./ScDrF4C_30.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s30=frequency)

pair04<-full_join(s11, s30, by="kmer")
head(pair04)

filtered_pair04<- pair04 %>% filter(s11 > 2*s30 | s30 > 2*s11)
nrow(filtered_pair04)

##### PAIR 05
s15<-read.table("./ScDrF4C_15.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s15=frequency)
s13<-read.table("./ScDrF4C_13.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s13=frequency)

pair05<-full_join(s15, s13, by="kmer")
head(pair05)

filtered_pair05<- pair05 %>% filter(s15 > 2*s13 | s13 > 2*s15)
nrow(filtered_pair05)

##### PAIR 06
s14<-read.table("./ScDrF4C_14.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s14=frequency)
s23<-read.table("./ScDrF4C_23.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s23=frequency)

pair06<-full_join(s14, s23, by="kmer")
head(pair06)

filtered_pair06<- pair06 %>% filter(s14 > 2*s23 | s23 > 2*s14)
nrow(filtered_pair06)

##### PAIR 07
s9<-read.table("./ScDrF4C_9.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s9=frequency)
s1632<-read.table("./ScDrF4C_1632.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s1632=frequency)

pair07<-full_join(s9, s1632, by="kmer")
head(pair07)

filtered_pair07<- pair07 %>% filter(s9 > 2*s1632 | s1632 > 2*s9)
nrow(filtered_pair07)

##### PAIR 08
s20<-read.table("./ScDrF4C_20.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s20=frequency)
s1633<-read.table("./ScDrF4C_1633.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s1633=frequency)

pair08<-full_join(s20, s1633, by="kmer")
head(pair08)

filtered_pair08<- pair08 %>% filter(s20 > 2*s1633 | s1633 > 2*s20)
nrow(filtered_pair08)


##### PAIR 09
s28<-read.table("./ScDrF4C_28.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s28=frequency)
s16<-read.table("./ScDrF4C_16.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s16=frequency)

pair09<-full_join(s28, s16, by="kmer")
head(pair09)

filtered_pair09<- pair09 %>% filter(s28 > 2*s16 | s16 > 2*s28)
nrow(filtered_pair09)

##### PAIR 10
s26<-read.table("./ScDrF4C_26.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s26=frequency)
s17<-read.table("./ScDrF4C_17.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s17=frequency)

pair10<-full_join(s26, s17, by="kmer")
head(pair10)

filtered_pair10<- pair10 %>% filter(s26 > 2*s17 | s17 > 2*s26)
nrow(filtered_pair10)

##### PAIR 11
s8<-read.table("./ScDrF4C_8.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s8=frequency)
s18<-read.table("./ScDrF4C_18.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s18=frequency)

pair11<-full_join(s8, s18, by="kmer")
head(pair11)

filtered_pair11<- pair11 %>% filter(s8 > 2*s18 | s18 > 2*s8)
nrow(filtered_pair11)

##### PAIR 12
s27<-read.table("./ScDrF4C_27.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s27=frequency)
s19<-read.table("./ScDrF4C_19.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s19=frequency)

pair12<-full_join(s27, s19, by="kmer")
head(pair12)

filtered_pair12<- pair12 %>% filter(s27 > 2*s19 | s19 > 2*s27)
nrow(filtered_pair12)


##### PAIR 13
s4<-read.table("./ScDrF4C_4.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s4=frequency)
s2<-read.table("./ScDrF4C_2.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s2=frequency)

pair13<-full_join(s4, s2, by="kmer")
head(pair13)

filtered_pair13<- pair13 %>% filter(s4 > 2*s2 | s2 > 2*s4)
nrow(filtered_pair13)

##### PAIR 14
s22<-read.table("./ScDrF4C_22.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s22=frequency)
s21<-read.table("./ScDrF4C_21.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s21=frequency)

pair14<-full_join(s22, s21, by="kmer")
head(pair14)

filtered_pair14<- pair14 %>% filter(s21 > 2*s22 | s22 > 2*s21)
nrow(filtered_pair14)

##### PAIR 15
s29<-read.table("./ScDrF4C_29.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s29=frequency)
s5<-read.table("./ScDrF4C_5.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s5=frequency)

pair15<-full_join(s29, s5, by="kmer")
head(pair15)

filtered_pair15<- pair15 %>% filter(s5 > 2*s29 | s29 > 2*s5)
nrow(filtered_pair15)

##### PAIR 16
s6<-read.table("./ScDrF4C_6.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s6=frequency)
s7<-read.table("./ScDrF4C_7.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s7=frequency)

pair16<-full_join(s6, s7, by="kmer")
head(pair16)

filtered_pair16<- pair16 %>% filter(s7 > 2*s6 | s6 > 2*s7)
nrow(filtered_pair16)


##### PAIR 17
s3<-read.table("./ScDrF4C_3.fa.mer_counts.dumps.larger100.col",sep="\t", header=T) %>%
  rename(s3=frequency)
s24<-read.table("./ScDrF4C_24.fa.mer_counts.dumps.larger100.col",sep="\t",header=T) %>%
  rename(s24=frequency)

pair17<-full_join(s24, s3, by="kmer")
head(pair17)

filtered_pair17<- pair17 %>% filter(s24 > 2*s3 | s3 > 2*s24)
nrow(filtered_pair17)



### Finally!
# Now, we join the data
df<-full_join(filtered_pair01, filtered_pair02, by="kmer") %>%
  full_join(filtered_pair03, by="kmer") %>%
  full_join(filtered_pair04, by="kmer") %>%
  full_join(filtered_pair05, by="kmer") %>%
  full_join(filtered_pair06, by="kmer") %>%
  full_join(filtered_pair07, by="kmer") %>%
  full_join(filtered_pair08, by="kmer") %>%
  full_join(filtered_pair09, by="kmer") %>%
  full_join(filtered_pair10, by="kmer") %>%
  full_join(filtered_pair11, by="kmer") %>%
  full_join(filtered_pair12, by="kmer") %>%
  full_join(filtered_pair13, by="kmer") %>%
  full_join(filtered_pair14, by="kmer") %>%
  full_join(filtered_pair15, by="kmer") %>%
  full_join(filtered_pair16, by="kmer") %>%
  full_join(filtered_pair17, by="kmer")

# Make sure ew remove missing data and make a hierarchical clustering.
par(mfrow=c(1,1))
cleaned_df<-df[complete.cases(df), ]
rownames(cleaned_df)<-cleaned_df$kmer
squeaky_cleaned_df<-cleaned_df[,-1]
transversed_squeaky_cleaned_df<-t(squeaky_cleaned_df)
dist_transversed_squeaky_cleaned_df <- dist(transversed_squeaky_cleaned_df)
hclust_avg <- hclust(dist_transversed_squeaky_cleaned_df)
plot(hclust_avg, main="13bp kmers at least 100x per chromossome, no NA; scalesia")


###
#rownames(df)<-df$kmer

#df<-df[,-1]

#transversed_df<-t(df)

#dist_transversed_df <- dist(transversed_df)
#hclust_avg <- hclust(dist_transversed_df)
#plot(hclust_avg, main="13bp kmers at least 100x per chromossome, NAs included; scalesia")

```
