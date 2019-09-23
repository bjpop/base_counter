'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 19 Sep 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

Detect variants
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pysam
import numpy as np
import pathlib
import csv
from scipy.stats import multinomial
import math


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
    description = 'Detect variants'
    parser = ArgumentParser(description=description)
    parser.add_argument('--counts', type=str, required=True,
        help='TSV file of base counts per genomic position')
    parser.add_argument('--regions', type=str, required=True,
        help='Bed file coordinates of genomic regions to consider')
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
            logging.warning(f"Unrecognised base: {base}")

    def __str__(self):
        return f"A: {self.A}, T:{self.T}, G:{self.G}, C:{self.C}"


def process_bam_file(bam_filename, regions):
    # pos -> (base -> Int)
    counts = {}
    sample = pathlib.PurePath(bam_filename).stem
    logging.info("Processing BAM file from %s", bam_filename)
    samfile = pysam.AlignmentFile(bam_filename, "rb" )
    for coord in regions:
        chrom, start, end = coord
        for pileupcolumn in samfile.pileup(chrom, start, end):
            for pileupread in pileupcolumn.pileups:
                this_pos = pileupcolumn.pos
                if start <= this_pos <= end:
                    if this_pos not in counts:
                        counts[this_pos] = Counts()
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # query position is None if is_del or is_refskip is set.
                        query_name = pileupread.alignment.query_name
                        this_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        counts[this_pos].increment_base_count(this_base)
    samfile.close()
    return counts


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
    with open(options.countsout, "w") as file:
        header = "pos\tA\tT\tG\tC"
        print(header, file=file)
        for pos in sorted(counts):
            print(f"{pos}\t{counts[pos].A}\t{counts[pos].T}\t{counts[pos].G}\t{counts[pos].C}", file=file)


def process_control_counts(options):
    # pos -> (base -> Int)
    counts = {}
    with open(options.counts) as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            this_pos = int(row['pos'])
            this_A = int(row['A'])
            this_T = int(row['T'])
            this_G = int(row['G'])
            this_C = int(row['C'])
            if this_pos not in counts:
                counts[this_pos] = {}
            counts[this_pos]['A'] = this_A
            counts[this_pos]['T'] = this_T
            counts[this_pos]['G'] = this_G
            counts[this_pos]['C'] = this_C
    return counts
            

def compare_case_control(case_counts, control_counts):
    results = []
    for pos in case_counts:
        this_case_counts = case_counts[pos]
        case_A = this_case_counts.A
        case_T = this_case_counts.T
        case_G = this_case_counts.G
        case_C = this_case_counts.C
        case_total = case_A + case_T + case_G + case_C
        if pos in control_counts:
            this_control_counts = control_counts[pos]
            control_A = this_control_counts['A']
            control_T = this_control_counts['T']
            control_G = this_control_counts['G']
            control_C = this_control_counts['C']
            control_total = control_A + control_T + control_G + control_C
            control_A_proportion = control_A / control_total
            control_T_proportion = control_T / control_total
            control_G_proportion = control_G / control_total
            control_C_proportion = control_C / control_total
            rv = multinomial(case_total, [control_A_proportion, control_T_proportion, control_G_proportion, control_C_proportion])
            probability = rv.pmf([case_A, case_T, case_G, case_C])
            if probability <= 0.001:
                results.append((probability, pos, case_A, case_T, case_G, case_C, control_A, control_T, control_G, control_C))
    return results 


def output_results(sample, results):
    for (prob, pos, case_A, case_T, case_G, case_C, control_A, control_T, control_G, control_C)  in sorted(results):
        case_total = case_A + case_T + case_G + case_C
        control_total = control_A + control_T + control_G + control_C
        case_A_prop = case_A / case_total
        case_T_prop = case_T / case_total
        case_G_prop = case_G / case_total
        case_C_prop = case_C / case_total
        control_A_prop = control_A / control_total
        control_T_prop = control_T / control_total
        control_G_prop = control_G / control_total
        control_C_prop = control_C / control_total
        delta_A = math.sqrt(case_A_prop ** 2 + control_A_prop)
        delta_T = math.sqrt(case_T_prop ** 2 + control_T_prop)
        delta_G = math.sqrt(case_G_prop ** 2 + control_G_prop)
        delta_C = math.sqrt(case_C_prop ** 2 + control_C_prop)
        out_str = "\t".join([sample] + [str(x) for x in [pos, prob, case_A, case_T, case_G, case_C, control_A, control_T, control_G, control_C, case_A_prop, case_T_prop, case_G_prop, case_C_prop, control_A_prop, control_T_prop, control_G_prop, control_C_prop, delta_A, delta_T, delta_G, delta_C]])
        print(out_str)

header = "\t".join(["sample", "pos", "prob", "case_A", "case_T", "case_G", "case_C", "control_A", "control_T", "control_G", "control_C", "case_A_prop", "case_T_prop", "case_G_prop", "case_C_prop", "control_A_prop", "control_T_prop", "control_G_prop", "control_C_prop", "delta_A", "delta_T", "delta_G", "delta_C"])

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    regions = get_regions(options) 
    print(header)
    for bam_filename in options.bam_files:
        sample = pathlib.PurePath(bam_filename).stem
        case_counts = process_bam_file(bam_filename, regions)
        control_counts = process_control_counts(options)
        results = compare_case_control(case_counts, control_counts)
        output_results(sample, results)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
