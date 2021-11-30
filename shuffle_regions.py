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
    log.info("Reading CDS regions from %s" % cds_file)
    print ("Reading CDS regions from %s" % cds_file)
    # Read the bed file of coding regions
    cds = pd.read_csv(cds_file, sep="\t", header=None)
    
    # Sort the coding regions and merge adjacent windows
    cds = pybedtools.BedTool.from_dataframe(cds).sort().merge()
    
    log.info("Read %i regions in CDS" % len(cds))
    print ("Read %i regions in CDS" % len(cds))

    return (cds)


def getConservedCodingRegions(conserved_file, cds_bed, outfile_stem,
                              log):

    log.info("Reading conserved regions from %s" % conserved_file)
    print ("Reading conserved regions from %s" % conserved_file)

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
    conserved_codons_df['start'] = conserved_codons_df['start'].astype(int)
    conserved_codons_df['end'] = conserved_codons_df['end'].astype(int)

    # Convert to a BedTools object
    conserved_codons = pybedtools.BedTool.from_dataframe(conserved_codons_df)

    log.info("Sorting and merging conserved regions")
    print ("Sorting and merging conserved regions")

    # Sort the BedTools object and merge adjacent codons to give conserved
    # windows
    conserved_windows = conserved_codons.sort().merge(d=1)
    log.info("Read %i conserved regions" % len(conserved_windows))
    print ("Read %i conserved regions" % len(conserved_windows))
    
    log.info("Intersecting coding and conserved regions")
    print("Intersecting coding and conserved regions")

    # Take only the parts of the conserved windows which are within the CDS
    conserved_cds = conserved_windows.intersect(cds_bed).sort().merge()
    log.info("Found %i conserved coding regions" % len(conserved_cds))
    print ("Found %i conserved coding regions" % len(conserved_cds))   
    
    conserved_cds.saveas("%s_conserved_cds.bed" % outfile_stem)
    # Return the conserved window and cds BedTools objects
    return (conserved_cds)


def getTestIntervals(test_file, cds_bed, log):
    log.info("Reading intervals to test from %s" % test_file)
    print("Reading intervals to test from %s" % test_file)
    # Read the test bed file
    test_bed = pybedtools.BedTool(test_file)
    
    log.info("Finding test intervals within the CDS")
    print("Finding test intervals within the CDS")
    
    # Get only the parts of the test bed file which are within the CDS
    test_bed_cds = test_bed.intersect(cds_bed).sort().merge()

    log.info("Found %i test intervals within the CDS" % len(test_bed_cds))
    print("Found %i test intervals within the CDS" % len(test_bed_cds))
    return (test_bed_cds)    


def makeNonCodingBed(genome_file, cds_bed, outfile_stem, log):
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
    
    non_cds_df = non_cds.to_dataframe()
    non_cds_df = non_cds_df.drop('start', 1)
    non_cds_df.to_csv('%s_csizes.txt' % outfile_stem, sep="\t", index=None)


def calcStats(bed, conserved_bed):
    
    by_conserved_peak = bed.intersect(
        conserved_bed, wb=True).sort().merge().to_dataframe()
    by_test_interval = bed.intersect(
        conserved_bed, wa=True).sort().merge().to_dataframe()
    total_overlap = bed.intersect(
        conserved_bed).sort().merge().to_dataframe()

    by_conserved_peak = by_conserved_peak.drop_duplicates()
    by_test_interval = by_test_interval.drop_duplicates()
    total_overlap = total_overlap.drop_duplicates()
    
    n_conserved_peaks = len(by_conserved_peak)
    n_test_intervals = len(by_test_interval)
    total_overlap = sum(total_overlap['end'] - total_overlap['start'])

    return ([n_conserved_peaks, n_test_intervals, total_overlap])
    

def addPercentiles(df):
    for string in ['N_Conserved_Overlapping_Test',
                   'N_Test_Overlapping_Conserved',
                   'Total_Overlap']:
        for ind in df.index.values:
            df.loc[ind, 'Percentile_%s' % string] = stats.percentileofscore(
                df.loc[df.index != ind, string],
                df.loc[ind, string])
    return (df)
    
def shuffle(test_bed, conserved_bed, cds_bed, outfile_stem, n_shuffles,
            log):

    # List to store the results
    rows = []
    log.info("Shuffling test file %i times" % n_shuffles)
    print ("Shuffling test file %i times" % n_shuffles)
    
    stat = calcStats(test_bed, conserved_bed)
    rows.append(['orig'] + stat)
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
        stat = calcStats(shuf, conserved_bed)

        rows.append(['shuffle_%i' % i] + stat)
    df = pd.DataFrame(rows, columns=['ID',
                                     'N_Conserved_Overlapping_Test',
                                     'N_Test_Overlapping_Conserved',
                                     'Total_Overlap'])
    df = addPercentiles(df)
    df.to_csv("%s_shuffle_results.tsv" % outfile_stem, sep="\t", index=None)
    return (df)


    
def plotResults(df, outfile_stem, log):
    f = plt.figure(figsize=(10, 20))
    a = f.add_subplot(131)
    shufs = df[1:]
    shuf_totals = shufs['Total_Overlap']
    orig_total = df.loc[0, 'Total_Overlap']
    sns.kdeplot(shuf_totals, ax=a)
    a.vlines(orig_total, 0, a.get_ylim()[1])
    xmin = min([orig_total, a.get_xlim()[0]]) * 0.9
    xmax = max([orig_total, a.get_xlim()[0]]) * 1.1
    a.set_xlim(xmin, xmax)
    sns.despine()
    plt.savefig("%s_distribution_total_overlap.png" % outfile_stem,
                dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="""Plot a sashimi plot of chimeric reads, e.g.\n\
                      plot_chimeric.py --annot_file EqTV_annotations.tsv \
                      --junction_file_1 junction_file_1.tsv \
                      --junction_file_2 junction_file_2.tsv \
                      --ref_fasta BEV.fasta \
                      --outfile test.png""")

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


    parser.add_argument("--outfile_stem", dest="outfile_stem", type=str,
                        help="Output file prefix")

    parser.add_argument("--n_shuffles", dest="n_shuffles", type=int,
                        default=1000,
                        help="Number of times to shuffle the intervals")

    
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
    try:
        os.mkdir('%s_temp' % (args.outfile_stem))
    except:
        pass

    pybedtools.set_tempdir('%s_temp' % args.outfile_stem)
    cds = getCDS(args.cds_file, log)

    if os.path.exists("%s_conserved_cds.bed" % args.outfile_stem):
        conserved_cds = pybedtools.BedTool(
            "%s_conserved_cds.bed" % args.outfile_stem)
    else:
        conserved_cds = getConservedCodingRegions(args.conserved_file,
                                                  cds,
                                                  args.outfile_stem,
                                                  log)
    test_cds = getTestIntervals(args.test_file, cds, log)
    
    makeNonCodingBed(args.genome_file, cds, args.outfile_stem, log)
    
    results = shuffle(test_cds, conserved_cds, cds, args.outfile_stem,
                      args.n_shuffles, log)
    plotResults(results, args.outfile_stem, log)
    pybedtools.helpers.cleanup()
    os.rmdir("%s_temp" % args.outfile_stem)


if __name__ == "__main__":
    main()
