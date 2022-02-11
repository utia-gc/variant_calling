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
    expected = ('lib_ID,sample_rep,fq1\n'
                'HSL-1,wt_control_rep1,data/HSL-1_R1.fastq.gz\n'
                'HSL-2,wt_control_rep2,data/HSL-2_R1.fastq.gz\n'
                'HSL-3,wt_DMSO_rep1,data/HSL-3_R1.fastq.gz\n'
                'HSL-4,wt_DMSO_rep2,data/HSL-4_R1.fastq.gz')
    assert out.strip() == expected


# --------------------------------------------------
def test_paired_end():
    """paired end"""

    rv, out = getstatusoutput(f'{prg} test_PE_design.csv')
    assert rv == 0
    expected = ('lib_ID,sample_rep,fq1,fq2\n'
                'HSL-1,wt_control_rep1,data/HSL-1_R1.fastq.gz,data/HSL-1_R2.fastq.gz\n'
                'HSL-2,wt_control_rep2,data/HSL-2_R1.fastq.gz,data/HSL-2_R2.fastq.gz\n'
                'HSL-3,wt_DMSO_rep1,data/HSL-3_R1.fastq.gz,data/HSL-3_R2.fastq.gz\n'
                'HSL-4,wt_DMSO_rep2,data/HSL-4_R1.fastq.gz,data/HSL-4_R2.fastq.gz')
    assert out.strip() == expected
