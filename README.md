### Download and run
Download the newest zip file with the binaries for your system from our release page [here](https://github.com/knights-lab/UTree/releases/latest). Unix (Mac and Linux) users: You may also have to set the files as executable with chmod +x utree-*

### Compilation (optional, only if you want to adjust internal parameters)

This set of compilation instructions are for Windows (msys/BASH environment) and Linux users only. Mac users will not be able to run these commands unless they have a real version of gcc (not the default clang). 

```
# Clone UTree
git clone https://github.com/knights-lab/UTree.git

# Make UTree
cd UTree
make
```

Instructions:
- Needs linearized, newline terminated (including trailing newline at end of file) FASTA format and a simple tab delimited map of [header] to [taxon].
- Mapping file specifics: [Taxon] can be any string without tab or null or newline. If doing rank-flexible (MAKE_GG or SEARCH_GG), it additionally needs to exactly respect the greengenes formatting standard including single-letter taxon labels, underscores, and spaces after the delimiting semicolons. A single label can map to multiple fasta sequences as long as the fasta sequences share a name (up to the terminating delimiter discussed below). It does not need to be sorted or deduplicated.
- FASTA file specifics: parsing stops at first encountered of the following: null, newline, underscore. Don't include tabs either as they will make the corresponding map unreadable (it is tab delimited). The string up to the end will be the label the program searches for in the mpping file to identify taxonomy information for the sequences. All sequences need to be accounted for in the map.

Took a trip down memory lane, brought up the code, and lo and behold... all sorts of small goodies to tweak. So here's a new version. With all info in the same place this time hopefully.
CHANGES:
- Supports up to 268 megabases per line in a FASTA file. That's bigger than human chromosome 1. Let me know if you need the limit removed
- Better defaults for compiler flags; a couple removed
- Slightly better usage, error messages (still needs a rewrite to be a "real"/"pure" program)

ERRATA/MISC/TODO:
- Should fail more gracefully if 0 leaves added
- COMPRESS makes trees way way bigger if they're small (e.g. toy goes from ~10kb to >60 MB)
- SEARCH is single threaded; multi-thread this guy!

TO COMPILE RANK-SPECIFIC:
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D BUILD -o utree-build
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D COMPRESS -o xtree-compress
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D SEARCH -o xtree-search

AND/OR RANK-FLEXIBLE (the COMPRESS stays the same):
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D BUILD_GG -o utree-buildGG
gcc -std=gnu11 -m64 -O3 itree.c -fopenmp -D SEARCH_GG -o xtree-searchGG

OPTIONAL COMPILER FLAGS (USE ON ALL BUILDS):
-D IXTYPE=uint32_t (if you have more than 64,000 unique labels or expect to extrapolate that many)
-D PACKSIZE=64 (if you want to use 64-mers instead of the default 32-mers. 4,8,and 16 are also valid k here)
-D PFBITS=26 (for larger desktops; only affects build. Basically this lets the program take more RAM to build a DB faster. Even numbers up to 32 are also possible for SUPER SERVERS; I'd suggest 28 for tminx)

QUERY BEHAVIOR COMPILER FLAGS
-D SLACK=X
-D SPARSITY=Y
per query controls.
So for a given query, recall that the program slides along the query one base at a time (full k-1 overlap) and checks that k-mer against the utree database to see if it uniquely matches something. 

SLACK is the number of times more the majority assignment must appear than the next most-count assignment.
So if a query matches 5 times to one species and 2 times to another, the ratio 5/2 is 2.5, which is less than a slack of 3, so its assignment is considered chimeric and it is not assigned. This only applies to rank-specific voting; optimal aufbau voting is implemented in search-gg (rank-flexible)

SPARSITY is the number of bases that have to elapse in the query before a series of identical consecutive calls can be counted twice.
So with SPARSITY=4, if a query matches species X, then slides a base and matches X, and another, and another, every time matching X in all 4 slides, it's only counted as one match to X. A series of matches can be interrupted and reset by matching to Y (a different species). A rule of thumb would be to set this between 1/4 and 1/8 of the database overlap size. Again, this only applies to rank-specific voting.
