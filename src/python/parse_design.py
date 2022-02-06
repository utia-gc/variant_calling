#!/usr/bin/env python3
"""
Author : Trevor F. Freeman <trvrfreeman@gmail.com>
Date   : 2022-02-06
Purpose: Parse input design file for nextflow pipeline
"""

import argparse


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
if __name__ == '__main__':
    main()
