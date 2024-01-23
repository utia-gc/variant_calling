---
title: Required params
layout: default
parent: Input / Output
---

# Required params
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
- TOC
{:toc}
</details>

## Problem

What params do I need to run the pipeline?

Bioinformatics pipelines require input of user data and settings.
We have designed `utia-gc/ngs` and pipelines built on `utia-gc/ngs` to require limited params that are simple to specify.

## Solution

We identified a few types of information that users must supply to be able to run the pipeline:

* Input sample files
  * `samplesheet` -- Sample files and metadata [formatted as a structured CSV file](https://github.com/utia-gc/ngs/wiki/Samplesheet-format).
* Reference assembly information
  * `genome` -- A URL or path to a reference genome fasta file. May be gzip compressed.
  * `annotations` -- A URL path to a reference annotations GTF file. May be gzip compressed.

## Usage
