# XTree (CrossTree) Documentation

## Overview

XTree (CrossTree) is a fast, k-mer-based sequence lookup and tallying tool designed for efficient analysis of metagenomic data. It offers a range of functionalities including genome coverage analysis, maximum-likelihood read redistribution, cross-tabulation of annotation hierarchies, and per-query annotation and statistics.

## Version

This documentation covers XTree version 2.00c.

## Features

- Total and unique genome coverage analysis
- Maximum-likelihood read redistribution
- Cross-tabulation of two annotation hierarchies (e.g., taxonomy and function)
- Per-query annotation and statistics
- Multiple granularity reporting: per-reference (genome-level), per-rank (redistributed), or per-rank (LCA [lowest common ancestor])
- High-speed processing (often 10-100x faster than alignment)
- Support for standard input, compressed files, and both FASTA and FASTQ formats
- Adamantium mode for robust k-mer analysis

## Usage

```
xtree {BUILD,ALIGN} [options]
```

### Options for both BUILD and ALIGN modes

- `--seqs <arg>`: Input sequences file
- `--log-out <arg>`: Log output file
- `--threads <arg>`: Number of threads to use
- `--db <arg>`: Database file (input for ALIGN, output for BUILD)

### BUILD Mode Options

- `--map <arg>`: Taxonomy mapping file
- `--comp <arg>`: Compression level (0-2)
- `--k <arg>`: K-mer size

### ALIGN Mode Options

- `--confidence <arg>`: Confidence threshold for LCA
- `--perq-out <arg>`: Per-query output file
- `--ref-out <arg>`: Reference output file
- `--tax-out <arg>`: Taxonomy output file
- `--cov-out <arg>`: Coverage output file
- `--orthog-out <arg>`: Orthogonal output file
- `--redistribute`: Enable redistribution algorithm
- `--fast-redistribute`: Enable single-pass redistribution
- `--shallow-lca`: Use shallow LCA algorithm
- `--copymem`: Copy database to local memory
- `--doforage`: Include all references in coverage table
- `--half-forage`: Include references with at least 2/3 of max k-mer matches
- `--no-adamantium`: Disable adamantium mode

## Input File Requirements

### Reference Sequences (BUILD mode)
- FASTA format
- Linearized: one line per header, one line per sequence
- Use `lingenome` tool to create compliant files

### Taxonomy Mapping (BUILD mode)
- Tab-delimited text file
- Columns: sequence name, taxonomy string, optional second taxonomy/function string
- Taxonomy strings are semicolon-delimited

### Query Sequences (ALIGN mode)
- FASTA or FASTQ format (gzipped or uncompressed)
- Maximum sequence length: 100MB

## Build Mode

Build mode creates an XTree database from reference sequences and taxonomic information.

### Key Options
- Compression levels (0-2): Higher levels recommended for large databases
- K-mer size: Affects sensitivity and specificity (default: 29)
- Log output: Provides detailed statistics on reference sequences

## Align Mode

Align mode processes query sequences against the XTree database.

### Key Features
- Confidence threshold for LCA assignment
- Redistribution algorithm for ambiguous mappings
- Adamantium mode for robust k-mer analysis

## Output Options

1. Per-query Output (`--perq-out`)
   - Detailed information for each query sequence

2. Reference Output (`--ref-out`)
   - Tabulated summary of matched reference sequences

3. Taxonomy Output (`--tax-out`)
   - Tabulated summary of matched taxonomies
   - LCA applied when `--redistribute` is not used

4. Coverage Output (`--cov-out`)
   - Detailed coverage statistics for reference sequences
   - Additional columns when adamantium mode is enabled

5. Orthogonal Output (`--orthog-out`)
   - Cross-tabulation of two hierarchies (e.g., taxonomy and function)

## Performance Considerations

- Multithreading support for improved performance
- Compression level 2 recommended for large databases
- `--copymem` option for improved database access speed
- `--shallow-lca` for faster LCA calculation in certain cases

## Limitations

- Maximum read length: 100MB
- Large databases may require significant memory
- Confidence threshold interpretation may change in future versions

## Conclusion

XTree provides a powerful and flexible tool for analyzing metagenomic data, offering a balance between speed and accuracy. Its various output options and optimization features make it suitable for a wide range of metagenomic studies and large-scale sequence analysis projects. Users should refer to the detailed addendum for in-depth information on database structure, input formats, output details, and advanced features.


# Addendum: Comprehensive XTree Documentation

## Database Structure

XTree uses a custom database format (`.xtr`) with the following structure:

1. Version byte [1] | compression level
2. Size of prefix [1]
3. Size of suffix [1]
4. Size of kmer_t [1]
5. Number of references [4]
6. Number of k-mers [8]
7. All prefix indices [1 << (2 * #2) x 8]
8. All k-mer data [(#2 + #4) * #6]
9. Number of adamantine k-mer targets [8]
10. All adamantine k-mer targets [(#2 + #4) * #9]
11. All reference lengths [#5 x 4]
12. Size of string data [8]
13. All strings [#9]
14. Number of Ix1 in map [4] [0 means skip rest of file]
15. String size for Ix1 [8]
16. Ix1 strings dump [#12]
17. Number of Ix2 in map [4] [can be 0/skipped if no h2 map]
18. String size of Ix2 [8]
19. Ix2 strings dump [#15]
20. HPairs[0] dump [num ref by 4]
21. HPairs[1] dump [num ref by 4]

## Input File Formats

### Reference Sequences (for BUILD mode)

- Format: FASTA
- Requirements:
  - Linearized format: exactly one line per header and one subsequent line for the entire sequence
  - No word wrap or blocks of sequence allowed
  - Use the `lingenome` tool to create compliant reference files from a directory of FASTA files

Example:
```
>Reference1
ATCGATCGATCGATCGATCG...
>Reference2
GCTAGCTAGCTAGCTAGCTA...
```

### Taxonomy Mapping File (for BUILD mode)

- Format: Tab-delimited text file
- Columns:
  1. Exact name of the reference sequence (without ">" character or file extensions)
  2. Semicolon-delimited list of strings denoting taxonomic (or other hierarchical) ranks, from most general to most specific
  3. (Optional) Another taxonomic (or other hierarchical) ranking in the same format as the second column

Example:
```
Ref1    Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli    Function1;Function2;Function3
Ref2    Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis    Function1;Function4;Function5
```

### Query Sequences (for ALIGN mode)

- Formats supported: FASTA or FASTQ
- Can be gzipped or uncompressed
- For FASTA format, sequences must be linearized (no line breaks within sequences)
- Maximum sequence length: 100MB

## Build Mode Details

### Compression Levels

- 0: No compression (default for small marker-gene databases)
- 1: Light compression (keep k-mers if the characters before the prefix are "A" then not "A")
- 2: Recommended for large whole-genome databases (keep k-mers if the characters before the prefix are "A" then "C")

### K-mer Size

The k-mer size (`-k` option) affects sensitivity and specificity:
- Longer k-mers: More specific but less sensitive
- Shorter k-mers: More sensitive but less specific
- Default: 29 (balances specificity and sensitivity)

### Log Output Option (`--log-out`)

When using the `--log-out` option in BUILD mode, XTree generates a tab-delimited log file with the following columns for each reference sequence:

1. Reference: Name of the reference sequence
2. UsableLen: Length of the sequence usable for k-mer generation (total length minus k-1 and ambiguous bases)
3. TotalKmers: Total number of k-mers in the sequence
4. UniqKmers: Number of k-mers unique to this sequence within the entire database
5. UniqKmers_SC: Number of k-mers unique to this sequence, including those appearing multiple times within it
6. AdamantineKmers (WeakKmers): Number of initially unique k-mers that become non-unique when considering single-nucleotide mutations

Note: The truly "adamantine" k-mers, those remaining unique even after mutation, would be the difference between UniqKmers and WeakKmers.

## Align Mode Details

### Confidence Threshold

The `--confidence` option sets the confidence level for LCA assignment:
- Scale: 0-1 (lowest to highest confidence)
- Interpretation (at compression level 2):
  - 0.25: ~80% statistical confidence
  - 0.5: 95%+ statistical confidence
  - 0.125: ~50% statistical confidence
- Recommended: 0.25 for whole-genome databases

### Redistribution Algorithm

The `--redistribute` option enables a maximum-likelihood algorithm to assign ambiguously mapping query sequences to the most likely reference or taxonomy.
- Uses a multi-pass algorithm to refine assignments
- Can be limited to a single pass with `--fast-redistribute`

### Lowest Common Ancestor (LCA) Usage

The LCA algorithm is applied in the following scenarios:

1. When generating the taxonomy output (`--tax-out`):
   - If `--redistribute` is not used, LCA is applied when a query matches multiple taxa at the chosen confidence level.
   - Broader-ranked taxa appear for queries that can't be confidently assigned to a specific taxon.

2. During per-query processing:
   - Used when a query doesn't have a clear, confident match to a single reference or taxon.
   - Applied when there are multiple top-scoring references or when the best match's k-mer count ratio is below the confidence threshold.

3. LCA calculation can be optimized using the `--shallow-lca` option.

### Adamantium Mode

- Enabled by default, can be disabled with `--no-adamantium`
- Identifies "adamantine" k-mers that remain unique even after single-nucleotide mutations
- Provides additional columns in the coverage output for more robust analysis

## Output Formats

### Per-query Output (`--perq-out`)

Tab-delimited file with the following columns:
1. Query sequence name
2. Reference sequence name
3. K-mers matched at assigned taxonomic level [most matched; other matched]
4. LCA taxonomy at the specified confidence level
5. Total k-mers matched

Example:
```
Query1  Ref1    [100;50]    Bacteria;Proteobacteria;Gammaproteobacteria    150
Query2  Ref2    [80;40]     Bacteria;Firmicutes;Bacilli    120
```

### Reference Output (`--ref-out`)

Tab-delimited file with two columns:
1. Reference sequence name
2. Total number of query sequences that matched that reference

Example:
```
Ref1    500
Ref2    300
```

### Taxonomy Output (`--tax-out`)

Tab-delimited file with two columns:
1. Taxonomy
2. Count

- If `--redistribute` is used, only full-rank taxa will appear
- If `--redistribute` is not used, LCA algorithm is applied, and broader-ranked taxa may appear

Example:
```
Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli    300
Bacteria;Firmicutes;Bacilli    200
```

### Coverage Output (`--cov-out`)

Tab-delimited file with the following columns:

1. Reference: Reference sequence name
2. Ref_len: Length of the reference sequence
3. Total_kmers: Total number of k-mers in the reference sequence
4. Total_unique: Number of unique k-mers in the reference sequence
5. Kmers_found: Total k-mers from query sequences that mapped to this reference
6. Unique_kmers_found: Unique k-mers from query sequences that mapped to this reference
7. Kmers_covered: Number of k-mers in the reference sequence that experienced at least one matching k-mer from a query sequence
8. Unique_kmers_covered: Number of unique k-mers in the reference sequence that experienced at least one matching k-mer from a query sequence
9. Proportion_covered: Decimal fraction of the reference sequence that experienced at least one matching k-mer from a matching read
10. Unique_proportion_covered: Decimal fraction of unique k-mers in the reference sequence that experienced at least one matching k-mer from a matching read
11. Reads_covered: Number of query sequences (not k-mers) that piled up onto the reference sequence
12. Exp: Expected coverage based on a Poisson distribution
13. Suspicion: A measure of potential contamination or strain-level variation
14. Coverage_est: Estimated coverage of the reference sequence
15. Num_reads_est: Estimated number of reads that mapped to the reference sequence

When "adamantium" mode is enabled (default), additional columns are included:

16. Adamantine_tot: Total number of adamantine k-mers in the reference sequence
17. Adamantium_found: Number of adamantine k-mers found in the query sequences
18. Adamantium_covered: Number of adamantine k-mers covered by query sequences
19. Adamantium_prop: Proportion of adamantine k-mers covered
20. XCov_adamantium: Cross-coverage estimate using adamantine k-mers
21. Reads_covered_adamantium: Estimated number of reads covered using adamantine k-mers

### Orthogonal Output (`--orthog-out`)

Tab-delimited file with three columns:
1. First hierarchy (e.g., taxonomy)
2. Second hierarchy (e.g., function)
3. Count

Example:
```
Bacteria;Proteobacteria;Gammaproteobacteria    Function1;Function2    150
Bacteria;Firmicutes;Bacilli    Function1;Function4    100
```

## Performance Optimization

- Use `--threads` to set the number of threads for parallel processing
- For large databases, use compression level 2 (`--comp 2`)
- The `--copymem` option can improve performance by copying the database to local memory
- Use `--shallow-lca` for faster LCA calculation in certain cases

## Limitations and Considerations

- Maximum read length is set to 100MB
- Large whole-genome databases may require significant memory
- The confidence threshold interpretation may change in future versions to match common expectations (e.g., 0.80 meaning 80% confidence)

## Error Handling

- The program includes basic error checking for input file formats and command-line arguments
- It will print error messages and exit if it encounters issues like malformed input files or invalid options

## Conclusion

XTree (CrossTree) is a powerful and flexible tool for metagenomic analysis, offering a balance between speed and accuracy. Its various output options and optimization features make it suitable for a wide range of metagenomic studies and large-scale sequence analysis projects.

Key features include:
- Fast k-mer-based sequence lookup and tallying
- Genome coverage analysis and maximum-likelihood read redistribution
- Cross-tabulation of annotation hierarchies (e.g., taxonomy and function)
- Per-query annotation and statistics
- Support for multiple reporting granularities

Users should pay close attention to input file requirements, particularly for reference sequences and taxonomy mapping in BUILD mode. The various output formats provide flexibility in analyzing results, from per-query details to broad taxonomic and functional summaries.

By understanding the details provided in this documentation, researchers can effectively leverage XTree's capabilities to process and analyze large-scale metagenomic datasets, gaining valuable insights into community composition and functional potential.
