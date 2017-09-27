### UWho?
UTree is an ultra-fast (up to 2x faster than Kraken) metagenomic profiling tool. Its database and RAM usage is an order of magnitude smaller than Kraken's, and it has greater specifity/precision with comparable recall after Bracken-like redistribution of lower-rank hits. It is distributed as a single static binary on Linux (and compiles for Mac and Windows) with no dependencies. It has a flexible compression parameter that allows the user to tune the size of the database at the expense of precision and recall. On a RefSeq representative genomes database containing over 5,000 prokaryotic genomes and over 7,000 viruses, the highest-compressed (L4) database consumes 500MB of memory and storage, making it feasible to run on a mobile device or phone, including running on long-read sequence data. The standard-compressed database (L2) of the same takes under 8GB of storage and RAM during alignment, and assigns dynamically interpolated taxonomy to reads at the rate of 16 million reads per minute on a 32-core Ivy Bridge server, or ~2 million reads per minute on an 8GB RAM Macbook Air. 

## To search a database:
[v1.5a] usage: xtree-searchGG compTree.ctr fastaToSearch.fa output.txt [threads] [RC]

Example: utree-searchGG rep82.gg.ctr seqs.fna classifications.txt 96 RC

Info:
The program will only use about 8.5GB of RAM (although it may appear to reserve more with many threads). 
The inputs are the CTR (compressed tree) and linearized FASTA file (alternating lines of header and sequence
with no line breaks in the sequences). The output is a simple text file listing the queries and their
interpolated classification, followed by (tabular) the number of k-mers matched in that read, the number among 
these that map to unique taxa records, then either a star (if only one unique record matched) or a per-level
support listing from kingdom to strain (k, p, c, o, f, g, s, t) in the format SUPPORTING_K-MERS;BAYESIAN_RANGE
and can be interpreted by simple division (2;3 is 66.7% confidence at a given level). 
Queries up to 16 megabases in length are supported, so you can assign taxonomy to whole bacterial genomes and long reads too.

e.g:
```
<QUERY>	<TAXONOMY>	<TOTAL_K-MERS_IN_READ>	<SUPPORT_KINGDOM>	<SUPPORT_PHYLUM>	<SUPPORT_CLASS>	<SUPPORT_ORDER>	<SUPPORT_FAMILY>	<SUPPORT_GENUS>	<SUPPORT_SPECIES>	<SUPPORT_STRAIN>
Even.40M.1.ninja.100000.4_318   k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__;t__    21      3       21;21   21;21   21;21   21;21   21;21   21;21   2;3     0;0
```

## To build a database:
[v2.0 SigNature Edition] usage: utree-buildGG input_fasta.fa labels.map output.ubt threads(0=auto) [complevel]
[v1.5a] usage: xtree-compress preTree.ubt compTree.ctr


Example: 
```
utree-buildGG myDatabase.fasta gg_format_taxonomy.txt temp.ubt 0 2
utree-compress temp.ubt myDatabase.ctr
rm temp.ubt
```


## Info:
This 2-step procedure constructs a database out of arbitrary reference sequences and their associated taxonomy. 
By default, UTree will not bother creating nodes below phylum level as they are too vague. The complevel value (2)
controls the sparsity of sampling the reference genomes. Larger values up to 4 (default is 1) decrease database 
size at the cost of sensitivity. Can be set to 0 for ultra-dense complete sampling. A value of 2 means the 
database size will be approximately equal to 1/3 the size of the original FASTA file. 1 = same size, 4 = 1/40.
Because the final compressed database is a complete, in-RAM searchable database, RAM requirements during search 
will almost exactly mirror the storage requirements for the final CTR database. 

File formats:
The fasta and taxonomy map formatting must be as follows, with sequence records up to 256 megabases (NO LINE BREAKS):
myDatabase.fasta (typical linearized fasta without line breaks):
>my first reference sequence
GACGATGCTAGCTGATCGATCGTGACTGCATGCTCAGTCGA
>my second reference sequence 
AGCGACGTAGCTGAGCA


gg_format_taxonomy.txt (seq name [spaces included!], tab, 8-level taxonomy with no spaces after semicolons):
```
my first reference sequence	k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__;t__
my second reference sequence	k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Sulfitobacter;s__Sulfitobacter_mediterraneus;t__Sulfitobacter_mediterraneus_KCTC_32188
```

## To cite: See release: https://github.com/knights-lab/UTree/releases/tag/v2.0c

## Misc
ERRATA/MISC/TODO:
- Should fail more gracefully if 0 leaves added
- COMPRESS makes trees way way bigger if original data is super small (e.g. toy tree goes from ~10kb to >60 MB)
- SEARCH can be sped up with more efficient I/O

TO COMPILE RANK-SPECIFIC (GCC and/or ICC and/or OpenMP+Clang+GNU extensions required to build):
```
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D BUILD -o utree-build
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D COMPRESS -o xtree-compress
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D SEARCH -o xtree-search
```


AND/OR RANK-FLEXIBLE (the COMPRESS stays the same):
```
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D BUILD_GG -o utree-buildGG
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D SEARCH_GG -o xtree-searchGG
```

OPTIONAL COMPILER FLAGS (USE ON ALL BUILDS):
-D IXTYPE=uint32_t (if you have more than 64,000 unique labels or expect to extrapolate that many)
-D PACKSIZE=64 (if you want to use 64-mers instead of the default 32-mers. 4,8,and 16 are also valid k here)
-D PFBITS=26 (for larger desktops; only affects build. Basically this lets the program take more RAM to build a DB faster. Even numbers up to 30 are also possible for SUPER SERVERS with > 128GB RAM)


QUERY BEHAVIOR COMPILER FLAGS (RANK-SPECIFIC)

-D SLACK=X

-D SPARSITY=Y

per query controls.

So for a given query, recall that the program slides along the query one base at a time (full k-1 overlap) and checks that k-mer against the utree database to see if it uniquely matches something. 


SLACK is the number of times more the majority assignment must appear than the next most-count assignment.
So if a query matches 5 times to one species and 2 times to another, the ratio 5/2 is 2.5, which is less than a slack of 3, so its assignment is considered chimeric and it is not assigned. This only applies to rank-specific voting; optimal aufbau voting is implemented in search-gg (rank-flexible)


SPARSITY is the number of bases that have to elapse in the query before a series of identical consecutive calls can be counted twice.

So with SPARSITY=4, if a query matches species X, then slides a base and matches X, and another, and another, every time matching X in all 4 slides, it's only counted as one match to X. A series of matches can be interrupted and reset by matching to Y (a different species). A rule of thumb would be to set this between 1/4 and 1/8 of the database overlap size. Again, this only applies to rank-specific voting.
