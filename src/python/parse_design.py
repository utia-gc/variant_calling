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
                        help='Input csv design file')


    args = parser.parse_args()

    return args


# --------------------------------------------------
def main():
    """main program"""

    args = get_args()


# --------------------------------------------------
def print_error(error, context):
    """Print error message"""

    error_msg = f"ERROR: Samplesheet -> {error}\n\tLINE: {context}"
    sys.exit(error_msg)


def check_header(design_file):
    """
    Check design file has proper header format:

    lib_ID,sample_name,replicate,reads1,reads2

    Note: reads2 column is optional for single end reads
    """


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

    assert check_header('./test/test_PE_design.csv') == "lib_ID,sample_name,replicate,reads1,reads2"
    assert check_header('./test/test_SE_design.csv') == "lib_ID,sample_name,replicate,reads1"


    
# --------------------------------------------------
if __name__ == '__main__':
    main()
