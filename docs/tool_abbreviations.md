# Tool Abbreviations

## Purpose

One of the key purposes of `SRAlign` (and pipelines built on `SRAlign`) is to simplify analysis.
Part of this simplification is making it easy to keep track of how each file has been produced, i.e. which tools were used to produce a file. 
In order to keep track of the steps for producing a file, this pipeline uses abbreviations appended onto file names to encode the tools used to produce each file.
Each major step is encoded in a brief, unique identifier.

## Unique tool identifiers

### Trim reads

#### `fastp`

fsp - `fastp`
