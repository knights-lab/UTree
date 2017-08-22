##To search a database:
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
Queries up to 16 megabases in length are supported, so you can assign taxonomy to whole genomes and long reads too.


e.g:
<QUERY>	<TAXONOMY>	<TOTAL_K-MERS_IN_READ>	<SUPPORT_KINGDOM>	<SUPPORT_PHYLUM>	<SUPPORT_CLASS>	<SUPPORT_ORDER>	<SUPPORT_FAMILY>	<SUPPORT_GENUS>	<SUPPORT_SPECIES>	<SUPPORT_STRAIN>
Even.40M.1.ninja.100000.4_318   k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__;t__    21      3       21;21   21;21   21;21   21;21   21;21   21;21   2;3     0;0

##To build a database:
[v2.0 SigNature Edition] usage: utree-buildGG input_fasta.fa labels.map output.ubt threads(0=auto) [complevel]
[v1.5a] usage: xtree-compress preTree.ubt compTree.ctr

Example: 
utree-buildGG myDatabase.fasta gg_format_taxonomy.txt temp.ubt 0 2
utree-compress temp.ubt myDatabase.ctr
rm temp.ubt

Info:
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
my first reference sequence	k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__;t__
my second reference sequence	k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Sulfitobacter;s__Sulfitobacter_mediterraneus;t__Sulfitobacter_mediterraneus_KCTC_32188

