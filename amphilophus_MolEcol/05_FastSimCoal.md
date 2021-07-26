FastSimCoal, a really powerful beast.


First order of business, we need the Site Frequency Spectrum (SFS). For this we need the vcf file (I included depth and missing data filters).
Second, we need a file formated as:
````
ind<tab>pop
ind<tab>pop
...
ind<tab>pop
````

To format this file, I did:
````
bcftools query -l final.FasSimCoal.onlyHAremoved.vcf > inds.ind.tsv 			# getting individuals as a columns
cat inds.ind.tsv | awk -F "_" '{print $(NF-1)"_"$(NF)}' > pops.pop.tsv			# getting populations as a columns

# Cleaning up some names - so it's easier to read.
cat pops.pop.tsv | sed 's/C_Ma/Citro_LakeMAN/; s/L_Ma/Labi_LakeMAN/; s/C_So/Citro_LakeNIC/; s/L_LV/Labi_LakeNIC/; s/L_So/Labi_LakeNIC/; s/C_Om/Citro_LakeNIC/; s/L_Om/Labi_LakeNIC/; s/C_LV/Citro_LakeN IC/; s/C_PD/Citro_LakeNIC/; s/L_PD/Labi_LakeNIC/' > populations.by.lakes.tsv

# Pasting files so it looks like I mentioned above: ind<tab>pop
paste inds.ind.tsv populations.by.lakes.tsv > easySFS.input.populations.tsv

````

Let's do this. We need python2 and the easySFS script:

````
module load python2

easySFS.py -i final.FasSimCoal.vcf -p easySFS.input.populations.tsv --preview
easySFS.py -i final.FasSimCoal.vcf -p easySFS.input.populations.tsv --proj 6,6,6,6
````

Now! FSC needs scenarios and models. See the folder 05_FastSimCoal_scenarios included as part of this folder.
We run each scenarios 10,000 times. Example of one scenario, on a array (SLURM_ARRAY_TASK_ID):

```
MODEL=no_geneflow


TEMPLATE=/cluster/projects/nn9408k/cerca/003_amphilophus_radseq/2_20190319_newdata/4_FasSimCoal/04_UpdatedModels/01_fsc_scenarios/$MODEL.tpl
ESTIMATES=/cluster/projects/nn9408k/cerca/003_amphilophus_radseq/2_20190319_newdata/4_FasSimCoal/04_UpdatedModels/01_fsc_scenarios/$MODEL.est
obs_file=/cluster/projects/nn9408k/cerca/003_amphilophus_radseq/2_20190319_newdata/4_FasSimCoal/01_inputfiles/2_SFS_convertion/2_without_HA/output/fastsimcoal2/final_MSFS.obs

# set up output
OUTPUT_DIR=/cluster/work/users/josece/fsc/$MODEL/run_log

# set up workdir
WORKDIR=/cluster/work/users/josece/fsc/$MODEL/runs/myarray_${SLURM_ARRAY_TASK_ID}
mkdir $WORKDIR

# We move to the working directory and copy the files there
cd $WORKDIR

echo "Running array number ${SLURM_ARRAY_TASK_ID}"

cp $obs_file ./${MODEL}_MSFS.obs
cp $TEMPLATE .
cp $ESTIMATES .

~/local_bin/fsc26_linux64/fsc26 -t $MODEL.tpl -e $MODEL.est -n 1 -C 10 -M -L 40 -c 20 -B 20 -m --multiSFS | tee run_${SLURM_ARRAY_TASK_ID}.log
```

Now, we need to see which of the 10,000 runs has the best likelihoods. For this we run a loop for each scenario and:
a) move in to the directory (cd $i/runs/)
b) open a second loop, to get every one of the 10,000 runs
c) go inside each of the 10,000 runs (cd myarray_$j/$i) ... copy the .bestlhoodsfile to a folder where we will have *all* the .bestlhood files (cp ${i}.bestlhoods ../../../run_log/$i.$j.bestlhoods)
d) return where we began (cd ../../)
```
for i in Managua_geneflow all_geneflow hybrid_origin_managua_geneflow hybrid_origin_nicaragua_geneflow modern_lineages_geneflow sympatric_geneflow Nicaragua_geneflow ancestral_geneflow hybrid_origin_managua_no_geneflow hybrid_origin_nicaragua_no_geneflow no_geneflow dispersal_geneflow dispersal_no_geneflow pulse_of_geneflow; do
	cd $i/runs/; for j in $(seq 1 10000); do cd myarray_$j/$i; cp ${i}.bestlhoods ../../../run_log/$i.$j.bestlhoods; cd ../..; done;
	cd ../../; done
````

Tired of loops? Here's another one for you. Now we have all .bestlhood files on a folder, we need to:
a) make a loop for each model, and move in a folder called run_log (cd $k/run_log)
b) inside we run another look for the 10,000 bestlhoods files.
c) we open the bestlhoods file, grep out the header and get the column we need, and estimate the "difference between the observed and estimated likelihoods" for each run (do cat $k.${i}.bestlhoods | grep -v MaxObsLhood | awk -v var=$i '{print var,$NF-$(NF-1)}' | cut -f 1 -d "."; done ), thus formating a "bestlhood.tsv" file
d) we go back to where we were (cd ../../)

````

for k in Managua_geneflow all_geneflow hybrid_origin_managua_geneflow hybrid_origin_nicaragua_geneflow modern_lineages_geneflow sympatric_geneflow Nicaragua_geneflow ancestral_geneflow hybrid_origin_managua_no_geneflow hybrid_origin_nicaragua_no_geneflow no_geneflow dispersal_geneflow dispersal_no_geneflow pulse_of_geneflow; do
	cd $k/run_log; for i in $(seq 1 10000); do cat $k.${i}.bestlhoods | grep -v MaxObsLhood | awk -v var=$i '{print var,$NF-$(NF-1)}' | cut -f 1 -d "."; done > ../bestlhood.$MODEL.tsv;
	cd ../../; done
v
````

Now, we identify the best run so we do some likelihood modeling:

`
sort -k2 bestlhood.$MODEL.tsv | head
`

In this case, I formatted something like:
Run_id is the ID from the 1-10,000 runs, and then the smallest difference between observed and estimated likelihoods.
`
<model>							<run_ID>	<smallest difference between observed and estimated likelihoods>
Managua_geneflow 					272			2538
Nicaragua_geneflow 					925			2802
all_geneflow 						976			2794
ancestral_geneflow 					756			2987
hybrid_origin_managua_geneflow 				631			2651
hybrid_origin_managua_no_geneflow 			60			2630
hybrid_origin_nicaragua_geneflow 			602			2598
hybrid_origin_nicaragua_no_geneflow 			695			2993
modern_lineages_geneflow 				503			2698
no_geneflow 						100			2688
sympatric_geneflow 					246			2739
pulse_of_geneflow					890			2777
dispersal_geneflow 					166			3072
dispersal_no_geneflow					875			2983
`

Now, we the .par file (generated by the best run) and the .obs file (generated in the beginning) and we run fsc 100 times.

Example:

``
cp /cluster/work/users/josece/fsc/dispersal_no_geneflow/runs/myarray_875/sympatric_geneflow/sympatric_geneflow_maxL.par .

MODEL=sympatric_geneflow
obs_file=/cluster/projects/nn9408k/cerca/003_amphilophus_radseq/2_20190319_newdata/4_FasSimCoal/01_inputfiles/2_SFS_convertion/2_without_HA/output/fastsimcoal2/final_MSFS.obs

OUTPUT_DIR=/cluster/work/users/josece/fsc/$MODEL/best_likelihood_distributions
bestpar=/cluster/projects/nn9408k/cerca/003_amphilophus_radseq/2_20190319_newdata/4_FasSimCoal/04_UpdatedModels/13_bestPARscenarios/${MODEL}_maxL.par

cd $OUTPUT_DIR

# set up workdir
WORKDIR=/cluster/work/users/josece/fsc/$MODEL/best_likelihood_distributions/myarray_${SLURM_ARRAY_TASK_ID}

mkdir $WORKDIR
cd $WORKDIR

#move to output dir and get things running
echo "Running array number ${TASK_ID}"

cp $bestpar .
cp $obs_file ./${MODEL}_maxL_jointMAFpop1_0.obs

~/local_bin/fsc26_linux64/fsc26 -i ${MODEL}_maxL.par -n1000000 -m -q -0
``

We get the best .lhoods file from each of the 100 runs, and format a file with the 100 ML distributions for each model:

`
for i in dispersal_no_geneflow dispersal_geneflow pulse_of_geneflow; do cd $i/; touch ${i}.lhoods; cd best_likelihood_distributions; for j in `ls`; do sed -n '2,3p' $j/*_maxL/*lhoods; done >> ../${i}.lhoods; cd ../..; done
`

Now plot this created file over R.
`
#fastSIMCOAL

setwd("C:/Users/Cerca/Desktop/ongoing_stuff_ambientedetrabalho/AMPHILOPHUS/amphi_fasSimCoal/04_UpdatedModels/bestLhood/")

Managua_geneflow<-read.table("Managua_geneflow.lhoods")
Nicaragua_geneflow<-read.table("Nicaragua_geneflow.lhoods")
all_geneflow<-read.table("all_geneflow.lhoods")
ancestral_geneflow<-read.table("ancestral_geneflow.lhoods")
dispersal_geneflow<-read.table("dispersal_geneflow.lhoods")
dispersal_no_geneflow<-read.table("dispersal_no_geneflow.lhoods")
hybrid_origin_managua_geneflow<-read.table("hybrid_origin_managua_geneflow.lhoods")
hybrid_origin_managua_no_geneflow<-read.table("hybrid_origin_managua_no_geneflow.lhoods")
hybrid_origin_nicaragua_geneflow<-read.table("hybrid_origin_nicaragua_geneflow.lhoods")
hybrid_origin_nicaragua_no_geneflow<-read.table("hybrid_origin_nicaragua_no_geneflow.lhoods")
modern_lineages_geneflow<-read.table("modern_lineages_geneflow.lhoods")
no_geneflow<-read.table("no_geneflow.lhoods")
pulse_of_geneflow<-read.table("pulse_of_geneflow.lhoods")
sympatric_geneflow<-read.table("sympatric_geneflow.lhoods")


par(mfrow=c(1,1))


boxplot(range = 0,
        Managua_geneflow$V1,
        Nicaragua_geneflow$V1,
        all_geneflow$V1,
        ancestral_geneflow$V1,
        dispersal_geneflow$V1,
        dispersal_no_geneflow$V1,
        hybrid_origin_managua_geneflow$V1,
        hybrid_origin_managua_no_geneflow$V1,
        hybrid_origin_nicaragua_geneflow$V1,
        hybrid_origin_nicaragua_no_geneflow$V1,
        modern_lineages_geneflow$V1,
        no_geneflow$V1,
        pulse_of_geneflow$V1,
        sympatric_geneflow$V1,
        ylab="Likelihood",xaxt="n")

lablist.x<-as.vector(c("Managua geneflow",
                       "Nicaragua geneflow",
                       "all geneflow",
                       "ancestral geneflow",
                       "dispersal geneflow",
                       "dispersal no geneflow",
                       "hybrid origin managua geneflow",
                       "hybrid origin managua no geneflow",
                       "hybrid origin nicaragua geneflow",
                       "hybrid origin nicaragua no geneflow",
                       "modern lineages geneflow",
                       "no geneflow",
                       "pulse of geneflow",
                       "sympatric geneflow"))
axis(1, at=seq(1, 14, by=1), labels = FALSE)
text(x = seq(1, 14, by=1), par("usr")[3] - 0.2, labels = lablist.x, srt = 45, pos = 1, xpd = TRUE)



library(tidyverse)
#install.packages("ggthemes")
library(ggthemes)
Managua_geneflow<-read.table("Managua_geneflow.lhoods") %>%
        rename(flow=V1)
Managua_geneflow$ID<-rep("Managua_geneflow.lhoods",nrow(Managua_geneflow))
Nicaragua_geneflow<-read.table("Nicaragua_geneflow.lhoods") %>%
        rename(flow=V1)
Nicaragua_geneflow$ID<-rep("Nicaragua_geneflow.lhoods",nrow(Nicaragua_geneflow))
all_geneflow<-read.table("all_geneflow.lhoods") %>%
        rename(flow=V1)
all_geneflow$ID<-rep("all_geneflow.lhoods",nrow(all_geneflow))
ancestral_geneflow<-read.table("ancestral_geneflow.lhoods") %>%
        rename(flow=V1)
ancestral_geneflow$ID<-rep("ancestral_geneflow.lhoods",nrow(ancestral_geneflow))
dispersal_geneflow<-read.table("dispersal_geneflow.lhoods") %>%
        rename(flow=V1)
dispersal_geneflow$ID<-rep("dispersal_geneflow.lhoods",nrow(dispersal_geneflow))
dispersal_no_geneflow<-read.table("dispersal_no_geneflow.lhoods") %>%
        rename(flow=V1)
dispersal_no_geneflow$ID<-rep("dispersal_no_geneflow.lhoods",nrow(dispersal_no_geneflow))
hybrid_origin_managua_geneflow<-read.table("hybrid_origin_managua_geneflow.lhoods") %>%
        rename(flow=V1)
hybrid_origin_managua_geneflow$ID<-rep("hybrid_origin_managua_geneflow.lhoods",nrow(hybrid_origin_managua_geneflow))
hybrid_origin_managua_no_geneflow<-read.table("hybrid_origin_managua_no_geneflow.lhoods") %>%
        rename(flow=V1)
hybrid_origin_managua_no_geneflow$ID<-rep("hybrid_origin_managua_no_geneflow.lhoods",nrow(hybrid_origin_managua_no_geneflow))
hybrid_origin_nicaragua_geneflow<-read.table("hybrid_origin_nicaragua_geneflow.lhoods") %>%
        rename(flow=V1)
hybrid_origin_nicaragua_geneflow$ID<-rep("hybrid_origin_nicaragua_geneflow.lhoods",nrow(hybrid_origin_nicaragua_geneflow))
hybrid_origin_nicaragua_no_geneflow<-read.table("hybrid_origin_nicaragua_no_geneflow.lhoods") %>%
        rename(flow=V1)
hybrid_origin_nicaragua_no_geneflow$ID<-rep("hybrid_origin_nicaragua_no_geneflow.lhoods",nrow(hybrid_origin_nicaragua_no_geneflow))
modern_lineages_geneflow<-read.table("modern_lineages_geneflow.lhoods") %>%
        rename(flow=V1)
modern_lineages_geneflow$ID<-rep("modern_lineages_geneflow.lhoods",nrow(modern_lineages_geneflow))
no_geneflow<-read.table("no_geneflow.lhoods") %>%
        rename(flow=V1)
no_geneflow$ID<-rep("no_geneflow.lhoods",nrow(no_geneflow))
pulse_of_geneflow<-read.table("pulse_of_geneflow.lhoods") %>%
        rename(flow=V1)
pulse_of_geneflow$ID<-rep("pulse_of_geneflow.lhoods",nrow(pulse_of_geneflow))
sympatric_geneflow<-read.table("sympatric_geneflow.lhoods") %>%
        rename(flow=V1)
sympatric_geneflow$ID<-rep("sympatric_geneflow.lhoods",nrow(sympatric_geneflow))

datacombined <- rbind(Managua_geneflow,
                      Nicaragua_geneflow,
                      all_geneflow,
                      ancestral_geneflow,
                      dispersal_geneflow,
                      dispersal_no_geneflow,
                      hybrid_origin_managua_geneflow,
                      hybrid_origin_managua_no_geneflow,
                      hybrid_origin_nicaragua_geneflow,
                      hybrid_origin_nicaragua_no_geneflow,
                      modern_lineages_geneflow,
                      no_geneflow,
                      pulse_of_geneflow,
                      sympatric_geneflow)

qplot( x=ID , y=flow , data=datacombined , geom=c("boxplot","jitter") , colour=ID) +
        theme_minimal()

qplot( x=ID , y=flow , data=datacombined , geom=c("boxplot","jitter") , colour=ID) +
        theme_hc()

`

The reviewer asked us for a AIC evaluation, so we did it to confirm our results!
For this we need the .est files from the begining.

`
cp ../01_fsc_scenarios/*est .

`

`
# And the bestlhood, from the best run out of 1000; e.g.:
cp /cluster/work/users/josece/fsc/dispersal_no_geneflow/runs/myarray_875/dispersal_no_geneflow/dispersal_no_geneflow.bestlhoods .
`

I got a script from Joana Meier to get the AIC assessment. It looks like this:

`
#!/usr/bin/env Rscript
# Joana Meier

# Usage: calculateAIC.sh modelprefix


# This script calculates AIC from fsc modeling results
# Run in the folder with the highest likelihood

# Read model name
args=commandArgs(TRUE)

# Checks if model name was given
if(length(args)<1){
  stop("ERROR: No input / model name given\nUsage: fsc-calculateAIC.R modelname")
}

# Check if model.bestlhoods file exists
if(file.exists(paste(args[1],".bestlhoods",sep=""))){
  bestlhoods<-read.delim(paste(args[1],".bestlhoods",sep=""))
}else{
  stop(paste("ERROR: Aborted. No file ",args[1],".bestlhoods file exists",sep=""))
}

# Check if model.est file exists
if(file.exists(paste(args[1],".est",sep=""))){
  est<-readLines(paste(args[1],".est",sep=""))
}else{
  stop(paste("ERROR: Aborted. No file ",args[1],".est file exists in this directory!\nUsage: fsc-calculateAIC.R modelname",sep=""))
}

# Count number of parameters
k<-(grep("RULES",est))-(grep("//all Ns are",est)+1)

# Calculate AIC
AIC<-2*k-2*(bestlhoods$MaxEstLhood/log10(exp(1)))

# Calculate delta-likelihood
deltaL<-bestlhoods$MaxObsLhood-bestlhoods$MaxEstLhood

# Output model.AIC file in simulation folder
write.table(cbind(deltaL,AIC),paste(args[1],".AIC",sep=""),row.names = F,col.names = T,sep = "\t",quote = F)
`


OK. So I did:
`
ml R/4.0.0-foss-2020a
./AIC.R pulse_of_geneflow

`
Now, we plot the results:

`
Managua_geneflow<-("16778.8511764356")
all_geneflow<-("17963.6260224633")
hybrid_origin_managua_geneflow<-("17307.327510856")
hybrid_origin_nicaragua_geneflow<-("17064.5567541613")
modern_lineages_geneflow<-("17519.096798611")
sympatric_geneflow<-("17709.6693073868")
Nicaragua_geneflow<-("17995.0720173849")
ancestral_geneflow<-("18848.1199271268")
hybrid_origin_managua_no_geneflow<-("17205.6011943391")
hybrid_origin_nicaragua_no_geneflow<-("18880.5568530252")
no_geneflow<-("17467.651086324")
dispersal_geneflow<-("19245.5363670848")
dispersal_no_geneflow<-("18829.7176670635")
pulse_of_geneflow<-("17877.4836018559")

modelnames<-c("Managua_geneflow",
              "all_geneflow",
              "hybrid_origin_managua_geneflow",
              "hybrid_origin_nicaragua_geneflow",
              "modern_lineages_geneflow",
              "sympatric_geneflow",
              "Nicaragua_geneflow",
              "ancestral_geneflow",
              "hybrid_origin_managua_no_geneflow",
              "hybrid_origin_nicaragua_no_geneflow",
              "no_geneflow",
              "dispersal_geneflow",
              "dispersal_no_geneflow",
              "pulse_of_geneflow")

modelAIC<-c("16778.8511764356",
            "17963.6260224633",
            "17307.327510856",
            "17064.5567541613",
            "17519.096798611",
            "17709.6693073868",
            "17995.0720173849",
            "18848.1199271268",
            "17205.6011943391",
            "18880.5568530252",
            "17467.651086324",
            "19245.5363670848",
            "18829.7176670635",
            "17877.4836018559")

library(ggplot2)
# Basic dot plot
p<-ggplot(df, aes(x=modelnames, y=modelAIC, fill = factor(modelnames))) +
  geom_dotplot(binaxis='y', stackdir='center')
p
`
