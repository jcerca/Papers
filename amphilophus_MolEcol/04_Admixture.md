Stop. [Hammer] Admixture time.


First, we will convert the data using plink v1.9. Remember we wantthe LD-filtered dataset.

```
VCF=amphilophus_snp_filter_hardpass.PRUNNED_MAF_COVERAGE.LDprunned.vcf.IndsWith0.5MissingnessRemoved.vcf
plink --vcf $VCF --recode12 --out admixture.input --allow-extra-chr --double-id --make-bed
```

Now; admixture wants chromossomes to be integers and not CCOE01 (like we have on this work) ... so we do some black magic, which will convert the chromosome Ids to numbers. This is dataset specific.

```
cat admixture.input.bim | sed -s "s/ENA|CCOE........|CCOE//" | sed -s "s/\.1//" | less > newinput.bim


head -n 5 newinput.bim # see the first column
01000001        ENA|CCOE01000001|CCOE01000001.1:83630   0       83630   A       C
01000001        ENA|CCOE01000001|CCOE01000001.1:128953  0       128953  C       G
01000001        ENA|CCOE01000001|CCOE01000001.1:130030  0       130030  T       G
01000001        ENA|CCOE01000001|CCOE01000001.1:154460  0       154460  C       T
01000001        ENA|CCOE01000001|CCOE01000001.1:154578  0       154578  G       T
mv newinput.bim admixture.input.bim 
```

Now, we create a new file for plotting:

```
awk -F "_" '{print $1"\t" $3 "\t" $4}' admixture.input.fam > ind.species.pop.tsv

head -n 6 ind.species.pop.tsv
13NIei01        C       Om
13NIei16        C       Om
13NIei17        C       Om
13NIei18        C       Om
13NIei19        C       Om
13NIei20        C       Om
```

Now the actual ADMIXTURE run:
````
INPUT=/usit/abel/u1/josece/003_amphilophus_radseq/2_20190319_newdata/3_admixture/1_fileconvertion/admixture.input.bed
MIX=/usit/abel/u1/josece/000_locally_installed_stuff/miniconda2/bin/admixture

for i in 2 3 4 5 6; do K=${i}; echo "This is my K for today = $i";
cd K$i ;
mkdir run1 run2 run3 run4 run5 ; 
cd run1 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ; 

cd run2 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ;

cd run3; 
$MIX --cv $INPUT $K -j4 | tee log; 
cd .. ;

cd run4 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ;

cd run5 ;
$MIX --cv $INPUT $K -j4 | tee log ;
cd .. ;

cd .. ;
done
````


... and plotting:

````
## examining Admixture data ##
rm(list = ls())

library(tidyverse)
library(gridExtra)

setwd("C:/Users/josece_adm/Desktop/AMPHILOPHUS/admixture/")

#we need the ind.species.pop.tsv file which is a "awk-modified" file from plink.fam
#the .Q files (they contain admixture proportions)
#the cv_error file :-)

# first read in and examine cross-validation error
cv_error<- tbl_df(read.table("./cv_error.log", sep="\t", header=T))

# read in individual data
#ind<-as.character(read.table("./admixture.input.fam")[,1])
myData <- tbl_df(read.table("./ind.species.pop.tsv", col.names = c("ind","species","site")))
#The ind.pop.species was created with a custom awk script based on the .fam file

# plot
a <- ggplot(cv_error, aes(K, error))+ geom_point() + geom_line()
a + ylab("Cross-valudation error") + theme_light()

## K2
# set up the dataset
K2 <- tbl_df(read.table("./admixture.input.2.Q",
                        col.names = c("one", "two")))
# add individuals
K2 <- mutate(K2, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K2 <- gather(K2, key = cluster, value = q, -ind, -site, -species)


## K3
# set up the dataset
K3 <- tbl_df(read.table("./admixture.input.3.Q",
                        col.names = c("one", "two", "three")))
# add individuals
K3 <- mutate(K3, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K3 <- gather(K3, key = cluster, value = q, -ind, -site, -species)


## K4
# set up the dataset
K4 <- tbl_df(read.table("./admixture.input.4.Q",
                        col.names = c("one", "two", "three", "four")))

# add individuals
K4 <- mutate(K4, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K4 <- gather(K4, key = cluster, value = q, -ind, -site, -species)


## K5
# set up the dataset
K5 <- tbl_df(read.table("./admixture.input.5.Q",
                        col.names = c("one", "two", "three", "four", "five")))
# add individuals
K5 <- mutate(K5, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K5 <- gather(K5, key = cluster, value = q, -ind, -site, -species)


##K6
# set up the dataset
K6 <- tbl_df(read.table("./admixture.input.6.Q",
                        col.names = c("one", "two", "three", "four", "five","six")))
# add individuals
K6 <- mutate(K6, ind = myData$ind, species = myData$species, site = myData$site)
# gather to plot
K6 <- gather(K6, key = cluster, value = q, -ind, -site, -species)


# set colour palettes
k2_cols <- c("dodgerblue", "firebrick1")
k3_cols <- c("dodgerblue", "gold",  "firebrick1")
k4_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen")
k5_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid")
k6_cols <- c("dodgerblue", "gold",  "firebrick1", "forestgreen", "darkorchid", "grey66")


#adding a lake column..
K2$lake<-K2$site
K2<-mutate(K2, lake = fct_recode(lake, "Nic"="Om", "Nic"="So", "Nic"="LV", "Nic"="PD"))

K3$lake<-K3$site
K3<-mutate(K3, lake = fct_recode(lake, "Nic"="Om", "Nic"="So", "Nic"="LV", "Nic"="PD"))


#by sites!
# make plot
a <- ggplot(K2, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = k2_cols)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
a <- a + facet_wrap(~site, scales = "free_x")
a <- a + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
a

# make plot
b <- ggplot(K3, aes(ind, q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = k3_cols)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
b <- b + facet_wrap(~site, scales = "free_x")
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
b

grid.arrange(a + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()),
             b + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()))



#by species!
a <- ggplot(K2, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = k2_cols)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
a <- a + facet_wrap(~species, scales = "free_x")
a <- a + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
a

# make plot
b <- ggplot(K3, aes(ind, q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = k3_cols)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
b <- b + facet_wrap(~species, scales = "free_x")
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
b

grid.arrange(a + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()),
             b + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()))



#by lake!
a <- ggplot(K2, aes(as.character(ind), q, fill = cluster)) + geom_bar(stat = "identity")
a <- a + scale_fill_manual(values = k2_cols)
a <- a + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
a <- a + facet_wrap(~lake, scales = "free_x")
a <- a + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
a

# make plot
b <- ggplot(K3, aes(ind, q, fill = cluster)) + geom_bar(stat = "identity")
b <- b + scale_fill_manual(values = k3_cols)
b <- b + xlab(NULL) + ylab("Ancestry") + scale_y_continuous(limits = c(0,1.01), expand = c(0, 0))
b <- b + facet_wrap(~lake, scales = "free_x")
b <- b + theme_light() + theme(legend.position = "none",
                               panel.border = element_blank(),
                               panel.grid = element_blank(),
                               axis.text.y = element_text(size = 10),
                               axis.text.x = element_text(size = 6, angle = 90),
                               axis.ticks.x = element_blank(),
                               strip.text.x = element_text(colour = "black", face = "bold"),
                               strip.background = element_rect(colour = "black", fill = "white"))
b

grid.arrange(a + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()),
             b + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()))
````
