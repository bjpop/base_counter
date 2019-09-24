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
    parser.add_argument('--regions', type=str, required=True,
        help='Bed file coordinates of genomic regions to consider')
    parser.add_argument('--coverageout', type=str, required=False,
        help='Output file for coverage per position')
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

def process_bam_files(options, regions):
    # pos -> (base -> Int)
    counts = {}
    # sample -> (pos -> int)
    coverage = {}
    num_queries_common = 0
    num_queries_read1_only = 0
    num_queries_read2_only = 0
    num_matching_bases = 0
    num_mismatching_bases = 0
    for bam_filename in options.bam_files:
        sample = pathlib.PurePath(bam_filename).stem
        if sample not in coverage:
            coverage[sample] = {}
        sample_coverage = coverage[sample]
        logging.info("Processing BAM file from %s", bam_filename)
        samfile = pysam.AlignmentFile(bam_filename, "rb" )
        for coord in regions:
            chrom, start, end = coord
            # ignore_orphans (bool) – ignore orphans (paired reads that are not in a proper pair). 
            # ignore_overlaps (bool) – If set to True, detect if read pairs overlap and only take the higher quality base. 
            for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True, stepper='samtools',
                                               ignore_overlaps=False, ignore_orphans=True,
                                               max_depth=1000000000, min_base_quality=options.basequal):
                read1s = {}
                read2s = {}
                this_pos_zero_based = pileupcolumn.pos
                this_pos_one_based = this_pos_zero_based + 1
                if this_pos_one_based not in counts:
                    counts[this_pos_one_based] = Counts()
                sample_coverage[this_pos_one_based] = pileupcolumn.n 
                # collect all bases in the current column, and assign them to either read1s or read2s
                for pileupread in pileupcolumn.pileups:
                    this_alignment = pileupread.alignment
                    query_name = this_alignment.query_name
                    mapping_quality = this_alignment.mapping_quality
                    alignment_length = this_alignment.query_alignment_length
                    # is_read1 and is_read2 are mutually exclusive
                    is_read1 = this_alignment.is_read1
                    is_read2 = this_alignment.is_read2
                    if mapping_quality >= options.mapqual and alignment_length >= options.alignlen:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            this_base = pileupread.alignment.query_sequence[pileupread.query_position]
                            if this_base in VALID_DNA_BASES:
                                if options.overlap:
                                    if is_read1:
                                        read1s[query_name] = this_base
                                    elif is_read2:
                                        read2s[query_name] = this_base
                                else:
                                    counts[this_pos_one_based].increment_base_count(this_base)
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
                            counts[this_pos_one_based].increment_base_count(read1_base)
                            counts[this_pos_one_based].increment_base_count(read1_base)
                            num_matching_bases += 1
                        else:
                            num_mismatching_bases += 1
        samfile.close()
    if options.overlap:
        logging.info(f"Number of bases covered by two reads: {num_queries_common}")
        logging.info(f"Number of bases covered by read 1 only: {num_queries_read1_only}")
        logging.info(f"Number of bases covered by read 2 only: {num_queries_read2_only}")
        logging.info(f"Number of bases covered by two reads that agree: {num_matching_bases}")
        logging.info(f"Number of bases covered by two reads that disagree: {num_mismatching_bases}")
    return counts, coverage


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
    counts, coverage = process_bam_files(options, regions)
    output_counts(options, counts)
    if options.coverageout is not None:
        output_coverage(options, coverage)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
