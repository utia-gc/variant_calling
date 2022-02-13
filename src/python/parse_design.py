#!/usr/bin/env python3
"""
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-06
Purpose: Parse input design file for nextflow pipeline
"""

import argparse
import pytest
import sys
import os


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Parse input design file for nextflow pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('design',
                        metavar='DESIGN',
                        nargs='+',
                        type=argparse.FileType('rt'),
                        help='Input csv design file')


    args = parser.parse_args()

    return args


# --------------------------------------------------
def main():
    """main program"""

    args = get_args()

    for in_fh in args.design:

        # read design file into a list of lists
        dsgn_in = [line.rstrip().split(',') for line in in_fh]

        # pull header into its own object
        header_in = dsgn_in.pop(0)


# --------------------------------------------------
def print_error(error, context):
    """Print error message"""

    error_msg = f"ERROR: Samplesheet -> {error}\n\tLINE: {context}"
    sys.exit(error_msg)


def check_header(header):
    """
    Check design file has proper header format:

    lib_ID,sample_name,replicate,reads1,reads2

    Note: reads2 column is optional to handle single end reads
    """

    VALID_HEADERS = ['lib_ID,sample_name,replicate,reads1'.split(','), 'lib_ID,sample_name,replicate,reads1,reads2'.split(',')]

    if header not in VALID_HEADERS:
        print_error(error="Missing or invalid header.", context=','.join(header))
    elif header == VALID_HEADERS[0]:
        return "single_end", "lib_ID,sample_rep,fq1"
    else:
        return "paired_end", "lib_ID,sample_rep,fq1,fq2"


# --------------------------------------------------
def test_print_error():
    """Test print_error"""

    error_str = f"ERROR: Samplesheet -> Invalid number of columns.\n\tLINE: wt_L1,HSL-1.fastq.gz"

    with pytest.raises(SystemExit) as out:
        print_error(error="Invalid number of columns.", context="wt_L1,HSL-1.fastq.gz")

    assert out.type == SystemExit
    assert out.value.code == error_str


def test_check_header(design_file):
    """Test check_header"""

    # test no or incorrect header
    error_str = f"ERROR: Samplesheet -> Missing or invalid header.\n\tLINE: HSL-1,wt_control,1,data/HSL-1_R1.fastq.gz,data/HSL-1_R2.fastq.gz"
    with pytest.raises(SystemExit) as out:
        check_header('./test/PE_design_noHeader.csv')
    
    assert out.type == SystemExit
    assert out.value.code == error_str


    # test correct header output
    assert check_header('./test/PE_design_good.csv') == "lib_ID,sample_rep,fq1,fq2"
    assert check_header('./test/SE_design_good.csv') == "lib_ID,sample_rep,fq1"


    
# --------------------------------------------------
if __name__ == '__main__':
    main()
