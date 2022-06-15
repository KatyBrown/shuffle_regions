# Shuffle Regions

Description
-----------
This script is used to compare sets of genomic intervals and determine if an annotation of interest overlaps more with an annotated set of conserved codons than would be expected by chance.

* Compare a set of test intervals and a set of conserved codons (both in bed format) and calculate:
  * The number of conserved intervals which overlap with at least one test interval
  * The number of test intervals which overlap at least one conserved interval
  * The total number of nucleotides of overlap between the conserved and test intervals
* Create a specified number of random sets of intervals within the coding seqeuence with the same size distribution as the test intervals
* Calculate the above statistics again for these shuffled intervals
* Calculate and tabulate the percentile of the original result vs the shuffled results
* Plot the original and shuffled results as a histogram

Usage
-----
```shuffle_regions.py --cds_file cds.bed --genome_file genome.bed --conserved_codon_file conserved.bed --test_file myintervals.bed --outfile_stem out --n_shuffles 1000```


This statement would create 1000 shuffled datasets, each with a set of intervals where the length distribution and number of intervals is equal to those in myintervals.bed, all of which are fully within the coding sequence defined in cds.bed.
The number of 1) intervals in conserved.bed which overlap with at least one interval in myintervals.bed, 2) intervals in myintervals.bed which overlap with at least one interval in conserved.bed and 3) nucleotides of overlap between conserved.bed and myintervals.bed would then be calculated.
These would be compared to the same three statistics comparing conserved.bed and the intervals in the shuffled bed files.
The results would be shown in a table in out_shuffled_results.tsv and plotted as histograms in out_distribution_total_overlap.png. 

Parameters
----------
`--cds_file`

Path to bed file containing co-ordinates of coding sequences

`--genome_file`

Path to bed file containing co-ordinates of the whole genome

`--conserved_codon_file`

Path to bed file containing the conserved codons of interest
                        
`--test_file`

Path to bed file containing the regions you want to intersect with the conserved regions

`--plot_color`

Colour for histograms

`--outfile_stem`

Output file prefix

`--n_shuffles`

Number of times to shuffle the intervals

`--min_overlap`

Minimum proportion of interval to be in the overlap to classify a region as overlapping


Output
------
`xxx_shuffled_results.tsv`
This table has 6 columns


* `ID` - row ID, either 'orig' for the original test file or shuffle_x for each shuffled dataset
* `N_Conserved_in_Test` - the number of conserved intervals which overlap with at least one test interval (or shuffled interval) by at least fraction `min_overlap`.
* `Prop_Conserved_in_Test`- the proportion of conserved intervals which overlap with at least one test interval (or shuffled interval) by at least fraction `min_overlap`.
* `N_Test_in_Conserved` - the number of test intervals (or shuffled intervals) which overlap with at least one conserved interval  by at least fraction `min_overlap`.
* `Prop_Test_in_Conserved` - the proportion of test intervals (or shuffled intervals) which overlap with at least one conserved interval by at least fraction `min_overlap`.
* `Total_Overlap` - the total number of nucleotides covered by both the test (or shuffled) intervals and the conserved intervals

`xxx_distribution_total_overlap.png`

Three histograms showing the comparison between the test intervals and the shuffled intervals for (left to right) 1) the proportion of conserved intervals which overlap with at least one test interval, 2) the proportion of test intervals which overlap with at least one conserved interval, 3) the total number of nucleotides which overlap between the conserved and test intervals.
