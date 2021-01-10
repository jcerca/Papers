This folder contains 4 population maps:

1. popmap_complete.tsv - including all individuals (3 species).
2. popmap_subterranea.tsv - including only Stygocapitella subterranea
3. popmap_josemariobrancoi.tsv - including only Stygocapitella josemariobrancoi
4. popmap_westheidei.tsv - including only Stygocapitella westheidei
5. popmap_cleaned.tsv - including individuals from the three species, but after cleaning. See the document ../01_stacks for cleaning.

Notice that the structure of population maps may have been corrupted as I copied them to github. Population maps should consist of two columns, separated by a tab:
<sample_01><tab><population_or_species_A>
<sample_02><tab><population_or_species_A>
<sample_03><tab><population_or_species_B>
<sample_04><tab><population_or_species_B>
...


It also contains the barcodes for the two sequencing lanes, barcode01, which corresponds to lane 01, and barcode02 which corresponds to lane02.
barcodes_lane01.tsv
barcodes_lane02.tsv

The structure of these consists of two columns, separated by a tab:
<Barcode_01><tab><sample_01>
<Barcode_02><tab><sample_02>
<Barcode_03><tab><sample_03>
<Barcode_04><tab><sample_04>
...
