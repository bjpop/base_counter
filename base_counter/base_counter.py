'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 19 Sep 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

Count frequency of DNA bases at genomic positions.
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pysam
import numpy as np
import pathlib
import csv
from intervaltree import Interval, IntervalTree
from collections import namedtuple


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "base_counter"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Count frequency of DNA bases at genomic positions'
    parser = ArgumentParser(description=description)
    parser.add_argument('--basequal', type=int, required=False,
        default=0, help='Minimum base quality. Bases below the minimum quality will be silently ignored')
    parser.add_argument('--mapqual', type=int, required=False,
        default=0, help='Minimum mapping quality. Reads below the minimum quality will be silently ignored')
    parser.add_argument('--alignlen', type=int, required=False,
        default=0, help='Minimum alignment length of reads. Reads with fewer bases aligned to the reference will be silently ignored')
    parser.add_argument('--ampliconprop', type=float, required=True,
        help='Minimum proportion of amplicon overlapped by a read')
    parser.add_argument('--regions', type=str, required=True,
        help='Bed file coordinates of genomic regions to consider')
    parser.add_argument('--amplicons', type=str, required=True,
        help='Coordinates of amplicons')
    parser.add_argument('--overlap', required=False, action='store_true', default=False,
        help='Only consider bases which are covered by forward and reverse reads and where both reads agree')
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('bam_files',
                        nargs='*',
                        metavar='BAM_FILE',
                        type=str,
                        help='Input BAM files')
    return parser.parse_args()


def get_regions(options):
    result = []
    with open(options.regions) as regions_file:
        for line in regions_file:
            fields = line.split()
            chrom, start, end = fields[:3]
            start = int(start)
            end = int(end)
            result.append((chrom, start, end))
    return result

'''
Amplicons are read in from BED file format.

Region coordinates are zero-based.

The first coordinate is within the region.
The second coordinate is one past the end of the region.
'''

Amplicon = namedtuple("Amplicon", ["start_outer", "end_outer", "name", "start_inner", "end_inner", "size_inner", "size_outer"])

def get_amplicons(options):
    amplicons = IntervalTree()
    with open(options.amplicons) as amplicons_file: 
        for line in amplicons_file:
            fields = line.split()
            chrom, start_outer, end_outer, name, start_inner, end_inner = fields[:6]
            start_outer = int(start_outer)
            end_outer = int(end_outer)
            start_inner = int(start_inner)
            end_inner = int(end_inner)
            size_inner = end_inner - start_inner
            size_outer = end_outer - start_outer
            this_amplicon = Amplicon(start_outer=start_outer, end_outer=end_outer,
                                     name=name, start_inner=start_inner, end_inner=end_inner,
                                     size_inner=size_inner, size_outer=size_outer)
            logging.info(f"Amplicon: {this_amplicon}")
            amplicons[start_outer:end_outer] = this_amplicon
    return amplicons 

'''
Check if a position in the genome is inside one of our amplicons.
A position may overlap more than one amplicon, so we have to choose
which amplicon to compare with. 
For each amplicon overlapping with the read alignment, we compute
the overlap between the read alignment and the inner segment of the
amplicon (the targetted bases of the amplicon), as a proportion
of the amplicon inner segment.
We choose the amplicon where this proportion is the largest.
If the proportion is greater than a threshold we check that the position
is within the inner segment of the amplicon.
'''

def is_pos_within_amplicon(min_amplicon_overlap_threshold, alignment, pos_zero_based, amplicons):
    # reference_end points to one past the last aligned residue
    # reference_start 0-based leftmost coordinate
    align_start = alignment.reference_start
    align_end = alignment.reference_end
    align_size = align_end - align_start
    max_overlap = None
    origin_amplicon = None
    for overlap in amplicons[align_start:align_end]:
        amplicon = overlap.data
        amplicon_start = amplicon.start_inner
        amplicon_end = amplicon.end_inner
        overlap_start = max(align_start, amplicon_start) 
        overlap_end = min(align_end, amplicon_end)
        overlap_size = overlap_end - overlap_start
        amplicon_size = amplicon.size_inner
        overlap_proportion = overlap_size / amplicon_size
        if max_overlap is None or overlap_proportion > max_overlap:
            max_overlap = overlap_proportion
            origin_amplicon = amplicon
    if origin_amplicon is None:
        # This position does not overlap any amplicon
        return False
    elif max_overlap < min_amplicon_overlap_threshold:
        # The best matching amplicon is not sufficiently overlapped by the read 
        return False
    elif (pos_zero_based < amplicon_start) or (pos_zero_based >= amplicon_end):
        # The position is outside of the amplicon's inner segment (e.g. could be in a primer)
        return False
    else:
        return True
    

class Counts(object):
    def __init__(self):
        self.A = 0
        self.T = 0
        self.G = 0
        self.C = 0

    def increment_base_count(self, base):
        if base == 'A':
            self.A += 1
        elif base == 'T':
            self.T += 1
        elif base == 'G':
            self.G += 1
        elif base == 'C':
            self.C += 1
        else:
            exit(f"Unrecognised base: {base}")

    def __str__(self):
        return f"A:{self.A}, T:{self.T}, G:{self.G}, C:{self.C}"


VALID_DNA_BASES = "ATGC"

def process_bam_files(options, regions, amplicons):
    total_reads = 0 
    total_failed_mapping_quality = 0 
    total_failed_base_quality = 0
    total_failed_align_length = 0 
    total_failed_overlapping_positions = 0
    total_failed_valid_dna_base = 0
    total_failed_amplicon = 0
    fieldnames = ['chrom', 'pos', 'sample', 'A', 'T', 'G', 'C', 'unfiltered coverage', 'filtered coverage', 'failed amplicon location', 'failed mapping quality', 'failed align length', 'failed base quality', 'failed valid DNA base', 'base absent', 'base not overlapped', 'failed overlapping positions']
    writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=fieldnames) 
    writer.writeheader()
    for bam_filename in options.bam_files:
        num_queries_common = 0
        num_queries_read1_only = 0
        num_queries_read2_only = 0
        samfile = pysam.AlignmentFile(bam_filename, "rb" )
        sample = pathlib.PurePath(bam_filename).stem
        logging.info(f"Processing BAM file from {bam_filename} for sample {sample}")
        for coord in regions:
            chrom, start, end = coord
            # ignore_orphans (bool) – ignore orphans (paired reads that are not in a proper pair). 
            # ignore_overlaps (bool) – If set to True, detect if read pairs overlap and only take the higher quality base. 
            for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True, stepper='samtools',
                                               ignore_overlaps=False, ignore_orphans=True,
                                               max_depth=1000000000):
                pos_failed_mapping_quality_reads = set()
                pos_failed_align_length_reads = set()
                pos_failed_base_quality = 0
                pos_failed_overlapping_positions = 0
                pos_failed_valid_dna_base = 0
                pos_failed_amplicon = 0
                pos_del_skip = 0
                num_pos_reads = 0 
                read1s = {}
                read2s = {}
                this_pos_zero_based = pileupcolumn.pos
                this_pos_one_based = this_pos_zero_based + 1
                counts = Counts()
                # the maximum number of reads aligning to this position before filtering
                unfiltered_coverage = pileupcolumn.nsegments
                # collect all bases in the current column, and assign them to either read1s or read2s
                for pileupread in pileupcolumn.pileups:
                    num_pos_reads += 1
                    count_this_read = True
                    this_alignment = pileupread.alignment
                    if not is_pos_within_amplicon(options.ampliconprop, this_alignment, this_pos_zero_based, amplicons):
                        pos_failed_amplicon += 1
                        count_this_read = False
                    query_name = this_alignment.query_name
                    mapping_quality = this_alignment.mapping_quality
                    alignment_length = this_alignment.query_alignment_length
                    # is_read1 and is_read2 are mutually exclusive
                    is_read1 = this_alignment.is_read1
                    is_read2 = this_alignment.is_read2
                    unique_query_name = query_name
                    if is_read1:
                        unique_query_name += "_R1"
                    if is_read2:
                        unique_query_name += "_R2"
                    # Skip reads that fail mapping quality
                    if mapping_quality < options.mapqual:
                        pos_failed_mapping_quality_reads.add(unique_query_name)
                        count_this_read = False 
                    # Skip reads that fail alignment length 
                    if alignment_length < options.alignlen:
                        pos_failed_align_length_reads.add(unique_query_name)
                        count_this_read = False 
                    if not pileupread.is_del and not pileupread.is_refskip:
                        this_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        this_base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                        if this_base_qual < options.basequal:
                            pos_failed_base_quality += 1
                            count_this_read = False 
                        if this_base not in VALID_DNA_BASES:
                            pos_failed_valid_dna_base += 1
                            count_this_read = False
                        if count_this_read:
                            if options.overlap:
                                if is_read1:
                                   read1s[query_name] = this_base
                                elif is_read2:
                                   read2s[query_name] = this_base
                            else:
                                counts.increment_base_count(this_base)
                    else:
                        pos_del_skip += 1
                if options.overlap:
                    read1_queries = set(read1s.keys())
                    read2_queries = set(read2s.keys())
                    queries_common = read1_queries.intersection(read2_queries)
                    num_queries_common += len(queries_common)
                    queries_read1_only = read1_queries - read2_queries
                    num_queries_read1_only += len(queries_read1_only) 
                    queries_read2_only = read2_queries - read1_queries
                    num_queries_read2_only += len(queries_read2_only) 
                    for this_query in queries_common:
                        read1_base = read1s[this_query]
                        read2_base = read2s[this_query]
                        if read1_base == read2_base:
                            # Increment the base twice because it appears on two reads
                            counts.increment_base_count(read1_base)
                            counts.increment_base_count(read1_base)
                        else:
                            pos_failed_overlapping_positions += 1
                num_pos_failed_mapping_quality = len(pos_failed_mapping_quality_reads)
                num_pos_failed_align_length = len(pos_failed_align_length_reads) 
                num_pos_not_overlapped = num_queries_read1_only + num_queries_read2_only
                output_row = {'chrom': 17,
                              'pos': this_pos_one_based,
                              'sample': sample,
                              'A': counts.A, 'T': counts.T, 'G': counts.G, 'C': counts.C,
                              #'unfiltered coverage': num_pos_reads,
                              'unfiltered coverage': unfiltered_coverage,
                              'filtered coverage': counts.A + counts.T + counts.G + counts.C,
                              'failed amplicon location': pos_failed_amplicon,
                              'failed mapping quality': num_pos_failed_mapping_quality,
                              'failed align length': num_pos_failed_align_length,
                              'failed base quality': pos_failed_base_quality,
                              'failed valid DNA base': pos_failed_valid_dna_base,
                              'base absent': pos_del_skip,
                              'base not overlapped': num_pos_not_overlapped,
                              'failed overlapping positions': pos_failed_overlapping_positions}
                writer.writerow(output_row) 
                total_failed_mapping_quality += num_pos_failed_mapping_quality
                total_failed_align_length += num_pos_failed_align_length
                total_failed_base_quality += pos_failed_base_quality
                total_failed_overlapping_positions += pos_failed_overlapping_positions
                total_failed_valid_dna_base += pos_failed_valid_dna_base
                total_failed_overlapping_positions += num_queries_read1_only + num_queries_read2_only
                total_failed_amplicon += pos_failed_amplicon
                total_reads += num_pos_reads
    samfile.close()
    logging.info(f"Total number of reads considered: {total_reads}")
    logging.info(f"Number of reads that failed amplicon overlap: {total_failed_amplicon}")
    logging.info(f"Number of reads that failed mapping quality threshold: {total_failed_mapping_quality}")
    logging.info(f"Number of reads that failed the alignment length threshold: {total_failed_align_length}")
    logging.info(f"Number of bases that failed base quality threshold: {total_failed_base_quality}")
    logging.info(f"Number of bases that failed valid DNA base: {total_failed_valid_dna_base}")
    if options.overlap:
        logging.info(f"Number of bases that failed the overlapping reads requirement: {total_failed_overlapping_positions}")


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def output_counts(options, counts):
    header = "pos\tA\tT\tG\tC"
    print(header)
    for pos in sorted(counts):
        print(f"{pos}\t{counts[pos].A}\t{counts[pos].T}\t{counts[pos].G}\t{counts[pos].C}")


def output_coverage(options, coverage):
    # pos -> [int]
    all_sample_coverage = {}
    for sample in coverage:
        this_sample_coverage = coverage[sample]
        with open(sample + ".tsv", "w") as file:
            print("\t".join(["pos","depth"]), file=file)
            for pos in sorted(this_sample_coverage):
                this_coverage = this_sample_coverage[pos]
                print(f"{pos}\t{this_coverage}", file=file)
                if pos not in all_sample_coverage:
                    all_sample_coverage[pos] = []
                all_sample_coverage[pos].append(this_coverage)
    with open(options.coverageout, "w") as file:
        print("\t".join(["pos","min","max","mean","std_dev"]), file=file)
        for pos in sorted(all_sample_coverage):
            this_counts = all_sample_coverage[pos]
            this_min = min(this_counts)
            this_max = max(this_counts)
            if len(this_counts) > 0:
                this_mean = sum(this_counts) / len(this_counts)
                this_dev = np.std(this_counts)
            else:
                this_mean = float('nan') 
                this_dev = float('nan') 
            print(f"{pos}\t{this_min}\t{this_max}\t{this_mean}\t{this_dev}", file=file)


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    regions = get_regions(options) 
    amplicons = get_amplicons(options)
    process_bam_files(options, regions, amplicons)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
