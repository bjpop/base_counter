'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 19 Sep 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

Get reference bases
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pysam


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "reference_bases"


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
    parser.add_argument('--regions', type=str, required=True,
        help='Bed file coordinates of genomic regions to consider')
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('reference',
                        metavar='REFERENCE_FILE',
                        type=str,
                        help='Input REFERANCE file')
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


def get_positions(reference_filename, regions):
    fasta_file = pysam.FastaFile(reference_filename)
    for chrom, start, end in regions:
        reference_sequence = fasta_file.fetch(reference=chrom, start=start, end=end) 
        for pos, base in zip(range(start, end), reference_sequence):
            # report the position in 1-based coordinates
            print(f"{chrom},{pos+1},{base}")
    fasta_file.close()


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    regions = get_regions(options) 
    get_positions(options.reference, regions)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
