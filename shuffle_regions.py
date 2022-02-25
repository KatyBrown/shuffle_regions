#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import argparse
import pybedtools
import logging
import os
import scipy.stats as stats


def getCDS(cds_file, log):
    '''
    Read the CDS regions from the CDS bed file --cds_file
    '''
    log.info("Reading CDS regions from %s" % cds_file)
    print ("Reading CDS regions from %s" % cds_file)
    # Read the bed file of coding regions
    cds = pd.read_csv(cds_file, sep="\t", header=None)
    
    # Sort the coding regions and merge adjacent windows
    cds = pybedtools.BedTool.from_dataframe(cds).sort().merge()
    
    log.info("Read %i regions in CDS" % len(cds))
    print ("Read %i regions in CDS" % len(cds))

    return (cds)


def getConservedCodingRegions(conserved_file, cds_bed, outfile_stem, log):
    '''
    - Read the bed file of conserved codons from --conserved_codon_file
    - Remove anything without a valid chromosome ID
    - Remove any codon not of length 3
    - Sort by start position
    - Merge adjacent codons to give conserved intervals
    - Ensure all the intervals are within the CDS defined in the CDS bed file
    - Save to outfile_stem_conserved_cds.bed
    '''

    log.info("Reading conserved codons from %s" % conserved_file)
    print ("Reading conserved codons from %s" % conserved_file)

    # Read the bed file of conserved codons
    conserved_codons_df = pd.read_csv(conserved_file,
                                      sep="\t", header=None,
                                      usecols=[0, 1, 2],
                                      names=['chrom', 'start', 'end'])
    
    # Remove anything without a valid chromosome ID
    conserved_codons_df = conserved_codons_df[
        ~conserved_codons_df['chrom'].isnull()]
    
    # Remove anything where codon length != 3
    conserved_codons_df = conserved_codons_df[
        conserved_codons_df['end'] - conserved_codons_df['start'] == 2]
    
    # Make sure start and end positions are integers
    conserved_codons_df['start'] = conserved_codons_df['start'].astype(int)
    conserved_codons_df['end'] = conserved_codons_df['end'].astype(int)

    # Convert to a BedTools object
    conserved_codons = pybedtools.BedTool.from_dataframe(conserved_codons_df)

    log.info("Sorting and merging conserved regions")
    print ("Sorting and merging conserved regions")

    # Sort the BedTools object and merge adjacent codons to give conserved
    # windows
    conserved_intervals = conserved_codons.sort().merge(d=1)
    log.info("Read %i conserved intervals" % len(conserved_intervals))
    print ("Read %i conserved intervals" % len(conserved_intervals))

    log.info("Saving conserved intervals")
    print ("Saving conserved intervals")    
    
    log.info("Intersecting coding and conserved regions")
    print("Intersecting coding and conserved regions")

    # Take only the parts of the conserved windows which are within the CDS
    # This should be all of them - it's just a safeguard
    conserved_cds = conserved_intervals.intersect(cds_bed).sort().merge()

    log.info("Found %i conserved coding regions" % len(conserved_cds))
    print ("Found %i conserved coding regions" % len(conserved_cds))   
    
    conserved_cds.saveas("%s_conserved_cds.bed" % outfile_stem)
    # Return the conserved window and cds BedTools objects
    return (conserved_cds)


def getTestIntervals(test_file, cds_bed, log):
    '''
    Get the test intervals from --test_file and intersect them with the CDS
    '''
    log.info("Reading intervals to test from %s" % test_file)
    print("Reading intervals to test from %s" % test_file)
    # Read the test bed file
    test_bed = pybedtools.BedTool(test_file)
    
    log.info("Finding test intervals within the CDS")
    print("Finding test intervals within the CDS")
    
    # Get only the parts of the test bed file which are within the CDS
    # Sort and merge just for safety
    test_bed_cds = test_bed.intersect(cds_bed).sort().merge()

    log.info("Found %i test intervals within the CDS" % len(test_bed_cds))
    print("Found %i test intervals within the CDS" % len(test_bed_cds))
    return (test_bed_cds)    


def makeNonCodingBed(genome_file, cds_bed, outfile_stem, log):
    '''
    Make a bed file of only non-coding regions of the genome to use
    with bedtools shuffle -excl to get only coding intervals.
    Used this rather than -incl becuase -incl doesn't work with -f
    which I want to use to limit the shuffle to only regions
    fully within the CDS.
    '''

    log.info("Reading genome from %s" % genome_file)
    print("Reading genome from %s" % genome_file)
    genome = pybedtools.BedTool(genome_file)

    log.info("Finding non-coding regions")
    print("Finding non-coding regions")  
    non_cds = genome.subtract(cds_bed)
    
    log.info("Saving non-coding regions as %s_non_coding.bed" % outfile_stem)
    print("Saving non-coding regions as %s_non_coding.bed" % outfile_stem)

    non_cds.saveas('%s_non_coding.bed' % outfile_stem)

    log.info("Saving chromosome sizes as %s_csizes.txt" % outfile_stem)
    print("Saving chromosome sizes as %s_csizes.txt" % outfile_stem)
    
    genome_df = genome.to_dataframe()
    genome_df = genome_df.drop('start', 1)
    genome_df.to_csv('%s_csizes.txt' % outfile_stem, sep="\t", index=None)


def calcStats(test_bed, conserved_bed, min_overlap):
    '''
    Intersect the conserved intervals with the test intervals and
    calculate:
        * n_conserved_intervals - the number of conserved intervals which
          overlap with at least one test interval
        * prop_conserved_intervals - the proportion of conserved intervals of
          which overlap with at least one conserved interval
        * n_test_intervals - the number of test intervals which overlap with
          at least one conserved interval.
        * prop_test_intervals - the proportion of test intervals which
          overlap with at least one conserved interval
        * total_overlap - the total number of nucleotides covered by
          both the test intervals and the conserved intervals
        * min_overlap - this proportion of at least one interval or the other
          must be in the intersection for n_conserved_intervals,
          prop_conserved_intervals,
          n_test_intervals and prop_test_intervals.
    '''

    # wb - return the intervals in b (conserved_bed) which overlap
    # the intervals in a (test_bed)
    # f - overlap must be at least this proportion of the interval in test_bed
    # F - overlap must be at least this proportion of the interval in conserved_bed
    # e - only F or f has to be True
    by_conserved_interval = conserved_bed.intersect(
        test_bed, wa=True, f=min_overlap, F=min_overlap,
        e=True).sort().merge().to_dataframe()
    
    # wa - return the intervals in a (test_bed) which overlap the intervals
    # in b (conserved_bed)
    # f - overlap must be at least this proportion of the interval in test_bed
    # F - overlap must be at least this proportion of the interval in conserved_bed
    # e - only F or f has to be True
    by_test_interval = test_bed.intersect(
        conserved_bed, wa=True, f=min_overlap, F=min_overlap,
        e=True).sort().merge().to_dataframe()
    
    
    # Return only regions which are in both files
    total_overlap = test_bed.intersect(
        conserved_bed).sort().merge().to_dataframe()

    # Remove duplicate rows from all three bed files
    by_conserved_interval = by_conserved_interval.drop_duplicates()
    by_test_interval = by_test_interval.drop_duplicates()
    total_overlap = total_overlap.drop_duplicates()
    
    # Calculate the statistics
    # Number of intervals in conserved_bed which contain at least
    # proportion min_overlap of an interval in test_bed
    n_conserved_intervals = len(by_conserved_interval)
    

    # Proportion of the total conserved intervals this applies to
    prop_conserved_intervals = n_conserved_intervals / len(conserved_bed)
    
    # Number of intervals in test_bed which overlap by at least proportion
    # min_overlap with an interval in conserved_bed
    n_test_intervals = len(by_test_interval)
    # Proportion of the test intervals this applies to
    
    if len(test_bed) != 0:
        prop_test_intervals = n_test_intervals / len(test_bed)
    else:
        prop_test_intervals = 0
    
    
    if len(total_overlap) != 0:
        # Total number of nucleotides covered by both the test intervals and
        # conserved intervals where at least proportion min_overlap of the
        # test interval overlaps
        total_overlap = sum(total_overlap['end'] - total_overlap['start'])
    else:
        total_overlap = 0

    return ([n_conserved_intervals, prop_conserved_intervals,
             n_test_intervals, prop_test_intervals, total_overlap])
    

def addPercentiles(df):
    '''
    For each of the count statistics in the results table, calculate the
    percentile of the data point compared to other points in the same column.
    '''
    # Iterate through all the columns
    for string in ['N_Conserved_in_Test',
                   'N_Test_in_Conserved',
                   'Total_Overlap']:
        # Iterate through all the rows
        for ind in df.index.values:
            # Calculate the percentile and put it in the dataframe
            df.loc[ind, 'Percentile_%s' % string] = stats.percentileofscore(
                df.loc[df.index != ind, string],
                df.loc[ind, string])
    return (df.round(5))


def shuffle(test_bed, conserved_bed, cds_bed, outfile_stem, n_shuffles,
            min_overlap, log):
    # List to store the results
    rows = []
    
    # Calculate overlap statistics for the input test file
    log.info("Calculating overlap between the test intervals and conserved intervals")
    print ("Calculating overlap between the test intervals and conserved intervals")

    stat = calcStats(test_bed, conserved_bed, min_overlap)
    
    log.info("""
Conserved intervals overlapping test, original: %i
Proportion of conserved intervals overlapping test, original: %.5f
Test intervals overlapping conserved intervals, original: %i
Proportion of test intervals overlapping conserved intervals, original: %.5f
Total number of overlapping nucleotides, original: %i
""" % tuple(stat))

    print("""
Conserved intervals overlapping test, original: %i
Proportion of conserved intervals overlapping test, original: %.5f
Test intervals overlapping conserved intervals, original: %i
Proportion of test intervals overlapping conserved intervals, original: %.5f
Total number of overlapping nucleotides, original: %i
""" % tuple(stat))
                
    # Store this result
    rows.append(['orig'] + stat)


    log.info("Shuffling test file %i times" % n_shuffles)
    print ("Shuffling test file %i times" % n_shuffles)
    
    # How often to log
    c = int(n_shuffles / 10)
    for i in np.arange(0, n_shuffles):
        if i % c == 0:
            log.info("Shuffling %i of %i" % (i, n_shuffles))
            print ("Shuffling %i of %i" % (i, n_shuffles))
        
        # Shuffle regions of the genome which are not in the non-coding
        # file (i.e. coding regions) to give regions of the same size as
        # those in the test file which are fully inside the CDS
        shuf = test_bed.shuffle(g='%s_csizes.txt' % outfile_stem,
                                excl='%s_non_coding.bed' % outfile_stem,
                                f=0)
        
        # Calculate how much the shuffled file overlaps with the conserved
        # intervals
        stat = calcStats(shuf, conserved_bed, min_overlap)
        
        # Store the results
        rows.append(['shuffle_%i' % i] + stat)
    
    log.info("Calculating overlap statistics for shuffled data")
    print ("Calculating overlap statistics for shuffled data")
    # Combine the results into a dataframe
    df = pd.DataFrame(rows, columns=['ID',
                                     'N_Conserved_in_Test',
                                     'Prop_Conserved_in_Test',
                                     'N_Test_in_Conserved',
                                     'Prop_Test_in_Conserved',
                                     'Total_Overlap'])
    log.info("""
Conserved intervals overlapping test, shuffled mean: %i
Proportion of conserved intervals overlapping test, shuffled mean: %.5f
Test intervals overlapping conserved intervals, shuffled mean: %i
Proportion of test intervals overlapping conserved intervals, shuffled mean: %.5f
Total number of overlapping nucleotides, shuffled mean: %i
""" % tuple([np.mean(df.loc[1:, x]) for x in df.columns[1:]]))

    print("""
Conserved intervals overlapping test, shuffled mean: %i
Proportion of conserved intervals overlapping test, shuffled mean: %.5f
Test intervals overlapping conserved intervals, shuffled mean: %i
Proportion of test intervals overlapping conserved intervals, shuffled mean: %.5f
Total number of overlapping nucleotides, shuffled mean: %i
""" % tuple([np.mean(df.loc[1:, x]) for x in df.columns[1:]]))

    log.info("Calculating percentiles")
    print ("Calculating percentiles")   
    # Add a percentile for each count statistic
    df = addPercentiles(df)
    log.info("""
Conserved intervals overlapping test, percentile: %.3f
Test intervals overlapping conserved intervals, percentile: %.3f
Total number of overlapping nucleotides, percentile: %.3f
""" % (df.loc[0, 'Percentile_N_Conserved_in_Test'],
       df.loc[0, 'Percentile_N_Test_in_Conserved'],
       df.loc[0, 'Percentile_Total_Overlap']))
    print("""
Conserved intervals overlapping test, percentile: %.3f
Test intervals overlapping conserved intervals, percentile: %.3f
Total number of overlapping nucleotides, percentile: %.3f
""" % (df.loc[0, 'Percentile_N_Conserved_in_Test'],
       df.loc[0, 'Percentile_N_Test_in_Conserved'],
       df.loc[0, 'Percentile_Total_Overlap']))

    df.to_csv("%s_shuffle_results.tsv" % outfile_stem, sep="\t", index=None)
    return (df)


    
def plotResults(df, outfile_stem, log, color='#ea0c31'):
    '''
    Plot density plots of the results
    '''
    # Create the figure
    log.info("Plotting results as density plots")
    print ("Plotting results as density plots")
    f = plt.figure(figsize=(20, 10))
    for i, column in enumerate(['N_Conserved_in_Test',
                                'N_Test_in_Conserved',
                                'Total_Overlap']):
        # Add a subplot
        a = f.add_subplot(1, 3, i+1)
        # Take the shuffled values
        shufs = df[1:]
        prop_column = column.replace("N", "Prop")
        shuf_totals = shufs[prop_column]
        # Take the original unshuffled value
        orig_total = df.loc[0, prop_column]
        rgb = matplotlib.colors.to_rgb(color)
        rgb_new = [r * 0.8 for r in rgb]
        #c = int(max(shuf_totals) / 20)
        # Plot the curve
        sns.histplot(shuf_totals, ax=a, color=rgb, edgecolor=None,
                     kde=True, kde_kws={'clip': None})
        #curve = a.lines[0].get_data()
        # Add a vertical line for the unshuffled point
        a.vlines(orig_total, 0, a.get_ylim()[1], color='#5b2c6f', lw=2)
        perc = df['Percentile_%s' % column].values[0]
        #a.text(orig_total*0.98, a.get_ylim()[1]*0.8,
        #       "%s\npercentile" % perc, ha='right', fontsize=8)
        #a.fill_between(curve[0], curve[1], color=color, alpha=0.3)
        #a.fill_between(curve[0],
        #               curve[1],
        #               where=curve[0] >= orig_total, color=color)
        # Make sure all the points are within the axis limits
        xmin = min([orig_total, a.get_xlim()[0]]) * 0.9
        xmax = max([orig_total, a.get_xlim()[1]]) * 1.1
        a.set_xlim(xmin, xmax)
        # Remove the right and top borders
        a.set_xlabel(prop_column.replace(
            "_", " ").replace("Prop", "Proportion"))
        a.set_ylabel("Frequency")
        for tick in a.get_xticklabels():
            tick.set_rotation('vertical')
        sns.despine()
    f.tight_layout()
    # Save and close the figure
    plt.savefig("%s_distribution_total_overlap.png" % outfile_stem,
                dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="""Calculate and plot the overlap between a set of test
                       intervals and a set of conserved codons, then shuffle
                       the test intervals and find and plot the percentile
                       of the original result vs the shuffled result""")

    parser.add_argument('--cds_file', dest='cds_file', type=str,
                        help='''Path to bed file containing co-ordinates of
                        coding sequences''')

    parser.add_argument("--genome_file", dest="genome_file", type=str,
                        help='''Path to bed file containing co-ordinates of
                        non-coding sequences''')

    parser.add_argument("--conserved_codon_file",
                        dest="conserved_file", type=str,
                        help='''Path to bed file containing the conserved
                        codons of interest''')                  
                        
    parser.add_argument("--test_file", dest="test_file", type=str,
                        help='''Path to bed file containing the regions
                        you want to intersect with the conserved regions''')

    parser.add_argument("--plot_color", dest="plot_color", type=str,
                        default='#7a9675', help="Colour for histograms")

    parser.add_argument("--outfile_stem", dest="outfile_stem", type=str,
                        help="Output file prefix")

    parser.add_argument("--n_shuffles", dest="n_shuffles", type=int,
                        default=1000,
                        help="Number of times to shuffle the intervals")

    parser.add_argument("--min_overlap", dest="min_overlap", type=float,
                        default=0.5,
                        help='''Minimum proportion of interval to be in the overlap
                        to classify a region as overlapping, this fraction can be of
                        either interval''')

    
    args = parser.parse_args()
    # Set up logger
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    logfile = "%s_log.txt" % args.outfile_stem
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)

    # Create a logging format
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)

    # add the handlers to the logger
    log.addHandler(handler)
    
    # pybedtools requires a temporary directory for output files
    try:
        os.mkdir('%s_temp' % (args.outfile_stem))
    except:
        pass
    pybedtools.set_tempdir('%s_temp' % args.outfile_stem)
    
    # Make a BedTool object of the CDS
    cds = getCDS(args.cds_file, log)
    
    # Make a BedTool object of the conserved intervals
    conserved_cds = getConservedCodingRegions(args.conserved_file,
                                              cds,
                                              args.outfile_stem,
                                              log)
    
    # Make a BedTool object of the test intervals
    test_cds = getTestIntervals(args.test_file, cds, log)
    
    # Make a bed file containing non-coding regions only
    makeNonCodingBed(args.genome_file, cds, args.outfile_stem, log)
    
    # Shuffle the test intervals and calculate the overlap with the conserved
    # intervals
    results = shuffle(test_cds, conserved_cds, cds, args.outfile_stem,
                      args.n_shuffles, args.min_overlap, log)
    
    # Plot the results
    plotResults(results, args.outfile_stem, log, args.plot_color)
    
    # Clean up temporary files
    pybedtools.helpers.cleanup()
    os.rmdir("%s_temp" % args.outfile_stem)


if __name__ == "__main__":
    main()
