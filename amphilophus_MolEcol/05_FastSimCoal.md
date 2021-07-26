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

`
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
`

We get the best .lhoods file from each of the 100 runs, and format a file with the 100 ML distributions for each model:

`
for i in `ls`; do cd $i/; touch ${i}.lhoods; cd best_likelihood_distributions; for j in `ls`; do sed -n '2,3p' $j/*_maxL/*lhoods; done >> ../${i}.lhoods; cd ../..; done
`

for i in dispersal_no_geneflow dispersal_geneflow pulse_of_geneflow; do cd $i/; touch ${i}.lhoods; cd best_likelihood_distributions; for j in `ls`; do sed -n '2,3p' $j/*_maxL/*lhoods; done >> ../${i}.lhoods; cd ../..; done

Now plot this over R.


#### AIC evaluation
#We need the est files
cp ../01_fsc_scenarios/*est .

# And the bestlhood, from the best run out of 1000
cp /cluster/work/users/josece/fsc/Managua_geneflow/runs/myarray_272/Managua_geneflow/Managua_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/Nicaragua_geneflow/runs/myarray_925/Nicaragua_geneflow/Nicaragua_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/all_geneflow/runs/myarray_976/all_geneflow/all_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/ancestral_geneflow/runs/myarray_756/ancestral_geneflow/ancestral_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/hybrid_origin_managua_geneflow/runs/myarray_631/hybrid_origin_managua_geneflow/hybrid_origin_managua_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/hybrid_origin_managua_no_geneflow/runs/myarray_60/hybrid_origin_managua_no_geneflow/hybrid_origin_managua_no_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/hybrid_origin_nicaragua_geneflow/runs/myarray_602/hybrid_origin_nicaragua_geneflow/hybrid_origin_nicaragua_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/hybrid_origin_nicaragua_no_geneflow/runs/myarray_695/hybrid_origin_nicaragua_no_geneflow/hybrid_origin_nicaragua_no_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/modern_lineages_geneflow/runs/myarray_503/modern_lineages_geneflow/modern_lineages_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/no_geneflow/runs/myarray_100/no_geneflow/no_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/sympatric_geneflow/runs/myarray_246/sympatric_geneflow/sympatric_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/pulse_of_geneflow/runs/myarray_890/pulse_of_geneflow/pulse_of_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/dispersal_geneflow/runs/myarray_166/dispersal_geneflow/dispersal_geneflow.bestlhoods .
cp /cluster/work/users/josece/fsc/dispersal_no_geneflow/runs/myarray_875/dispersal_no_geneflow/dispersal_no_geneflow.bestlhoods .

ml R/4.0.0-foss-2020a
# Note, I had to modify Joana's script. I basically changed the shebang from "#! /usr/bin/Rscript" to "#!/usr/bin/env Rscript"


#Notice, we need the .est file and the .bestlhood file together.
./AIC.R Managua_geneflow
./AIC.R all_geneflow
./AIC.R hybrid_origin_managua_geneflow
./AIC.R hybrid_origin_nicaragua_geneflow
./AIC.R modern_lineages_geneflow
./AIC.R sympatric_geneflow
./AIC.R Nicaragua_geneflow
./AIC.R ancestral_geneflow
./AIC.R hybrid_origin_managua_no_geneflow
./AIC.R hybrid_origin_nicaragua_no_geneflow
./AIC.R no_geneflow
./AIC.R dispersal_geneflow
./AIC.R dispersal_no_geneflow
./AIC.R pulse_of_geneflow

#Use my plotting script :)
