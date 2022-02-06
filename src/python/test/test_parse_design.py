#!/usr/bin/env python3
"""
tests for parse_design.py
"""


import os
import re
from subprocess import getstatusoutput

prg = '../parse_design.py'

# --------------------------------------------------
def test_exists():
    """exists"""

    assert os.path.isfile(prg)


# --------------------------------------------------
def test_usage():
    """usage"""

    for flag in ['-h', '--help']:
        rv, out = getstatusoutput(f'{prg} {flag}')
        assert rv == 0
        assert re.match("usage", out, re.IGNORECASE)


# --------------------------------------------------
def test_single_end():
    """single end"""

    rv, out = getstatusoutput(f'{prg} test_SE_design.csv')
    assert rv == 0
    expected = ('wt_L1_rep1,data/library_123.fastq.gz\n'
                'wt_L1_rep2,data/library_456.fastq.gz')
    assert out.strip() == expected


# --------------------------------------------------
def test_paired_end():
    """paired end"""

    rv, out = getstatusoutput(f'{prg} test_PE_design.csv')
    assert rv == 0
    expected = ('wt_L1_rep1,data/library_123_R1.fastq.gz,data/library_123_R2.fastq.gz\n'
                'wt_L1_rep2,data/library_456_R1.fastq.gz,data/library_456_R2.fastq.gz')
    assert out.strip() == expected
