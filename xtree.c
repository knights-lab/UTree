// Gabriel Al-Ghalith. Efficient characterization of orthogonal metagenomic
//  taxonomy and pathway coverage with CrossTree. 2018.
#define VER "CrossTree v0.92i by Gabe"
#define VNO 1
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <stdint.h>
#ifndef _OPENMP
	#define omp_get_max_threads() 1
	#define omp_get_thread_num() 0
	#include <time.h>
	#define omp_get_wtime() ((double)clock()/CLOCKS_PER_SEC)
	#define omp_set_num_threads(x) x
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#ifndef kmer_t
	#define kmer_t uint32_t
#endif
#ifndef rix_t
	#define rix_t uint32_t
#endif
#define AMBIG 4
#include <math.h>
#include <zlib.h>

void * huge_malloc(size_t n) {
	void *ptr = 0;
	posix_memalign(&ptr, 1 << 21, n);
	madvise(ptr, n, MADV_HUGEPAGE);
	return ptr;
}
void * huge_calloc(size_t n) {
	void *ptr = huge_malloc(n);
	memset(ptr,0,n);
	return ptr;
}
const uint8_t CONV[32]  = {4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4};
const uint8_t RCONV[32] = {4,3,4,2,4,4,4,1,4,4,4,4,4,4,4,4,4,4,4,4,0,0,4,4,4,4,4,4,4,4,4,4};

#pragma pack(1)
typedef struct {
	kmer_t sfx;
	rix_t rix;
} KPod;

//typedef struct {uint32_t h1, h2;} hpair_t;

int strcmp_ptr(const void* a, const void* b)
	{ return strcmp(*(char**)a,*(char**)b); }

typedef union {uint32_t a[4]; __uint128_t n; uint8_t b[16];} MasterBin_t;
int binCmp(const void *a, const void *b) {
	MasterBin_t *A = *(MasterBin_t **)a, *B = *(MasterBin_t **)b;
	uint64_t valA = *(uint64_t *)(A->b+4), valB = *(uint64_t *)(B->b+4);
	return valA < valB ? -1 : valB < valA;
}

static inline uint32_t prefix_to_num(char *s, int nl, int *err, uint32_t kshift) {
	uint32_t nib = 0, k = kshift;
	for (uint32_t i = 0; i < nl; ++i, k-=2) {
		uint32_t c = CONV[31 & s[i]];
		if (c == AMBIG) {*err = 1; return i;}
		nib |= c << k;
	}
	return nib;
}
static inline uint32_t prefix_to_num_RC(char *s, int nl) {
	uint32_t nib = 0, k = 0, c;
	for (uint32_t i = 0; i < nl; ++i, k+=2) 
		c = RCONV[31 & s[i]],
		nib |= c << k;
	return nib;
}

static inline kmer_t dna_to_num(char *s, int nl, int *err, kmer_t kshift) {
	kmer_t nib = 0, k = kshift;
	for (uint32_t i = 0; i < nl; ++i, k-=2) {
		kmer_t c = CONV[31 & s[i]];
		if (c == AMBIG) {*err = 1; return i;}
		nib |= c << k;
	}
	return nib;
}
static inline kmer_t dna_to_num_RC(char *s, int nl) {
	kmer_t nib = 0, k = 0, c;
	for (uint32_t i = 0; i < nl; ++i, k+=2) 
		c = RCONV[31 & s[i]],
		nib |= c << k;
	return nib;
}

static inline uint64_t binsearch_str(char **Strings, char *key, uint64_t N) {
	uint64_t lo = 0, hi = N;
	while (lo < hi) {
		uint64_t mid = lo + ((hi-lo) >> 1);
		int cmp = strcmp(key,Strings[mid]);
		if (cmp > 0) lo = mid+1;
		else if (cmp < 0) hi = mid;
		else return mid;
	}
	return -1;
}
// Key string delimited by tab, null, newline...
static inline uint64_t binsearch_str_d(char **Strings, char *key, uint64_t N) {
	uint64_t lo = 0, hi = N;
	while (lo < hi) {
		uint64_t mid = lo + ((hi-lo) >> 1);
		char *a = key, *b = Strings[mid];
		while (*b && *a==*b) ++a, ++b;
		if (!*b && (!*a || *a == '\n' || *a=='\t')) return mid;
		if (*a < *b) hi = mid;
		else lo = mid+1;
	}
	return -1;
}

// returns first possible match as if null existed at position key_len, matched as far as possible
// A helper function will translate position via Map_h* minus ranks missing. 
static inline uint64_t binsearch_str_L(char **Strings, char *key, uint64_t N, uint32_t key_len) {
	uint64_t lo = 0, hi = N;
	while (lo < hi) {
		uint64_t mid = lo + ((hi - lo) >> 1);
		int cmp = 0;
		//strcmp(key,Strings[mid]);
		uint32_t i = 0; char *ref = Strings[mid];
		for (; i < key_len && ref[i]; ++i) 
			if (key[i] != ref[i]) break;
		if (i == key_len /* && key[i] */ && ref[i]) cmp = -1;
		else cmp = key[i] - ref[i];
		//printf("Compared key:\n%s\nagainst ref:\n%s\nuntil pos %u: %d [lo=%lu, mid=%lu, hi=%lu]\n",
		//	key,ref,key_len,cmp,lo,mid,hi);
		//if (!cmp) return mid;
		if (cmp > 0) lo = mid + 1;
		else if (cmp < 0) hi = mid;
		else return mid;
	}
	return lo; 
}
static inline int cmpfunc(const void *a, const void *b) {
	return *(uint64_t*)a < *(uint64_t*)b ? -1 :
		*(uint64_t*)b < *(uint64_t*)a;
}
static inline int u16cmp(const void *a, const void *b) {
	return *(uint16_t*)a < *(uint16_t*)b ? -1 : 
		*(uint16_t*)b < *(uint16_t*)a;
}
static inline int kpackcmp(const void *a, const void *b) {
	KPod *k1 = (KPod *)a, *k2 = (KPod *)b;
	if (k1->sfx < k2->sfx) return -1;
	if (k2->sfx < k1->sfx) return 1;
	if (k1->rix < k2->rix) return -1;
	if (k2->rix < k1->rix) return 1;
	return 0;
}

/* static inline KPod * WBS_k(KPod *KP, uint64_t Lx, kmer_t k) {
	KPod *p=KP;
	while (Lx) {
		size_t w = (Lx >> 1) + 1;
		if (p[w].sfx < k) p+=w, Lx-=w; 
		else if (p[w].sfx == k) return p+w;
		else Lx = w-1;
	}
	return p->sfx==k ? p : 0;
} */

static inline uint64_t LBS_k(KPod *KP, uint64_t N, kmer_t k) {
	uint64_t L = 0, R = N;
	while (L < R) {
		uint64_t m = (L + R) >> 1;
		if (KP[m].sfx < k) L = m+1;
		else R = m;
	}
	return KP[L].sfx == k ? L : -1;
}

static inline uint64_t get_queries(gzFile in, uint8_t **QBucket, char **HBucket, uint8_t *head, uint8_t *line, int qChunk, uint64_t szmax) {
	uint64_t nq = 0; uint8_t *QB_ptr = *QBucket; char *H_ptr = *HBucket;
	uint8_t *eol; int len;
	while (nq < qChunk && QB_ptr - *QBucket <= szmax) {
		if (!gzgets(in,head,(1 << 20)-1)) break;
		eol = strchr(head,'\n');
		*eol = 0; 
		len = eol-head+1;
		memcpy(H_ptr,head+1,len-1);
		HBucket[nq] = H_ptr;
		H_ptr += len-1;
		
		if (!gzgets(in,line,(1 << 28)-1)) break;
		eol = strchr(line,'\n');
		*eol = 0; 
		len = eol-line+1;
		memcpy(QB_ptr,line,len);
		QBucket[nq] = QB_ptr;
		QB_ptr += len;
		if (*head == '@' && (!gzgets(in,line,1 << 28) || 
			!gzgets(in,line,1 << 28))) break;
		++nq;
	}
	return nq;
}

#define USAGE1 "USAGE: xtree {BUILD,ALIGN} [options]\n  "
#define USAGE2 USAGE1 "Options for both BUILD and ALIGN, with args: {seqs,log-out,threads,db}\n"
#define USAGE3 USAGE2 "BUILD Options\n  With args: {map,comp,k,db-out} <arg>\n"
#define USAGE4 USAGE3 "ALIGN Options\n  With args: {confidence,perq-out,ref-out,tax-out,cov-out,orthog-out}\n"
#define USAGE USAGE4 "  Without args: {redistribute,shallow-lca,copymem}"
int main(int argc, char *argv[]) {
	puts(VER);
	char *dbPath = 0, *seqPath = 0, *ixPath = 0, *covPath = 0;
	char *perqPath = 0, *taxPath = 0, *orthogPath = 0, *refPath = 0, *logPath = 0;
	int doBuild = 0, threads = omp_get_max_threads(), comp = 0, kchoice = 0;
	int doFullLCA = 1, doRedist = 0, doFastRedist = 0, doCopyMem = 0;
	double conf = 0.33; // Reasonable default for compression lv 2
	uint32_t nUniqMatches = 0;
	
	// Robust parser
	for (int a = 1; a < argc; ++a) {
		if (!strcmp(argv[a],"BUILD")) doBuild = 1; 
		else if (!strcmp(argv[a],"--map")) ixPath = argv[++a];
		else if (!strcmp(argv[a],"--comp")) comp = atoi(argv[++a]); //I
		else if (!strcmp(argv[a],"--k")) kchoice = atoi(argv[++a]); //I
		
		else if (!strcmp(argv[a],"ALIGN")) doBuild = 0; 
		else if (!strcmp(argv[a],"--confidence")) {
			double ctemp = atof(argv[++a]);
			if (ctemp <= 1) conf = ctemp, printf("Setting confprop %f\n",conf);
			else nUniqMatches = ctemp, printf("Setting min uniq ref matches to %u\n",nUniqMatches);
		}
		else if (!strcmp(argv[a],"--perq-out")) perqPath = argv[++a];
		else if (!strcmp(argv[a],"--ref-out")) refPath = argv[++a];
		else if (!strcmp(argv[a],"--tax-out")) taxPath = argv[++a];
		else if (!strcmp(argv[a],"--cov-out")) covPath = argv[++a];
		else if (!strcmp(argv[a],"--orthog-out")) orthogPath = argv[++a];
		else if (!strcmp(argv[a],"--redistribute")) doRedist = 1; //NA
		else if (!strcmp(argv[a],"--fast-redistribute")) doRedist = 1, doFastRedist = 1; //NA
		else if (!strcmp(argv[a],"--shallow-lca")) doFullLCA = 0; //NA
		else if (!strcmp(argv[a],"--copymem")) doCopyMem = 1; //NA
		
		
		// Options for both BUILD and ALIGN
		else if (!strcmp(argv[a],"--seqs")) seqPath = argv[++a]; 
		else if (!strcmp(argv[a],"--log-out")) logPath = argv[++a];
		else if (!strcmp(argv[a],"--threads")) threads = atoi(argv[++a]);
		else if (!strcmp(argv[a],"--db") || !strcmp(argv[a],"--db-out")) 
			dbPath = argv[++a];
		
		else {printf("Unrecognized option: %s\n",argv[a]); exit(1);}
	}
	threads = threads > 256 ? 256 : threads;
	omp_set_num_threads(threads); 
	printf("Using %d thread(s)\n",threads);
	if (argc < 4) {puts(USAGE); exit(1);}
	
	if (doBuild) {
		uint32_t PL = 13, SL = sizeof(kmer_t)*4;
		uint64_t K = PL+SL;
		if (comp) printf("Setting compression level to %d\n",comp);
		
		if (kchoice) K = kchoice;
		SL = K - PL;
		if (K < PL || !SL || SL > sizeof(kmer_t)*4) {printf("Bad K! [%lu]\n",K); exit(1);}
		printf("Building DB with K=%lu [PL %d, SL %d]\n",K,PL,SL);
		
		uint32_t kpre_shf = PL*2-2;
		kmer_t kpst_shf = SL*2-2;
		uint32_t pre_bshf = 32-(PL*2), pre_bshf_2 = pre_bshf+2;
		kmer_t pst_bshf = sizeof(kmer_t)*8-(SL*2), pst_bshf_2 = pst_bshf+2;

		FILE *in = fopen(seqPath,"rb");
		if (!in) {printf("ERROR: bad FASTA input: %s\n",seqPath); exit(2);}
		struct stat sb; int fno = fileno(in); fstat(fno,&sb);
		uint64_t fsz = sb.st_size;
		char *Raw = mmap(0, fsz+16, PROT_READ, MAP_SHARED, fno, 0);
		madvise(Raw,fsz,MADV_SEQUENTIAL);
		
		// We need: num records, num valid k-mers (each!), header locs
		if (Raw[0] != '>') {puts("Uh oh. Input FASTA looks fishy."); exit(2);}
		
		// Set up the tank for input -- up to 1 billion refs allowed
		uint64_t nbins = 1<<(2*PL),
			*Offsets = huge_malloc(((uint64_t)1<<30)*sizeof(*Offsets)),
			*Nibs = huge_calloc((nbins+1)*sizeof(*Nibs));
		if (!Offsets || !Nibs) {puts("ERROR:OOM Offsets"); exit(3);}
		uint32_t ns = 0;
		double wtime = omp_get_wtime();
		#pragma omp parallel for
		for (uint64_t z = 0; z < fsz; ++z) {
			if (Raw[z] > 64 && Raw[z-1] == '\n') {
				uint32_t ix;
				#pragma omp atomic capture
				ix = ns++;
				uint64_t x = z, y = z;
				while (Raw[y] && Raw[y] != '\n') ++y;
				
				uint32_t num, slideSafe = 0;
				Offsets[ix] = x; // Offsets contain start of sequence
				while (x + K <= y) {
					int err=0, a = 0; 
					if (slideSafe) {
						uint32_t c = CONV[31 & Raw[x+PL-1]];
						if (c==AMBIG) {x+=PL; slideSafe=0; continue;}
						num = (num << pre_bshf_2 >> pre_bshf) | c;
					} else {
						num = prefix_to_num(Raw+x,PL,&err,kpre_shf);
						if (err) {x+= num+1; slideSafe = 0; continue;}
						slideSafe = 1;
					}
					while (a < comp && !CONV[31 & Raw[x+a-comp]]) ++a; // compression
					if (a == comp)
						#pragma omp atomic
						++Nibs[num];
					++x;
				}
			}
		}
		printf("There were %u records here (%f s)\n",ns,omp_get_wtime()-wtime);
		if (ns > 65535 && sizeof(rix_t) == 2) 
			{puts("ERROR: too many refs (>65K)"); exit(2);}
		//Offsets = realloc(Offsets,ns*sizeof(*Offsets));
		qsort(Offsets,ns,sizeof(*Offsets),cmpfunc);
		
		// Create the data structures
		uint64_t totWords = 0;
		for (uint64_t i = 0; i < nbins; ++i)
			totWords += Nibs[i];
		printf("In total, we need a structure that is %lu large.\n",totWords);
		printf("Now allocating %f GB of RAM...\n",(double)totWords*sizeof(KPod)/1073741824);
		wtime = omp_get_wtime();
		
		KPod *KGrid = huge_malloc(totWords*sizeof(*KGrid));
		uint64_t *KIx = huge_calloc((nbins+1)*sizeof(*KIx));
		if (!KIx) {puts("OOM:KIx"); exit(3);}
		for (uint64_t i = 1; i < nbins; ++i) 
			KIx[i] = KIx[i-1] + Nibs[i-1];
		for (uint64_t i = 0; i < nbins; ++i) Nibs[i] = KIx[i];
		
		#ifdef WSL
		for (uint64_t i = 0; i < totWords; i+=4096/sizeof(*KGrid)) 
			KGrid[i].sfx=1;
		#endif
		
		uint32_t mask32 = ((uint64_t)1 << (2*PL)) - 1;
		kmer_t maskK = ((kmer_t)1 << (2*SL)) - 1;
		if (SL==sizeof(kmer_t)*4) maskK = -1;
		#pragma omp parallel for schedule(dynamic)
		for (uint32_t i = 0; i < ns; ++i) {
			uint64_t x = Offsets[i];
			uint64_t y = x; while (Raw[y] && Raw[y] != '\n') ++y;
			int slideSafe = 0; uint32_t num; kmer_t kmer;
			while (x + K <= y) {
				int err=0, a = 0;
				if (slideSafe) {
					uint32_t c = CONV[31 & Raw[x+PL-1]];
					if (c==AMBIG) {x+=PL; slideSafe=0; continue;}
					//num = (num << pre_bshf_2 >> pre_bshf) | c;
					num = ((num << 2) | c) & mask32;
					kmer_t k = CONV[31 & Raw[x+K-1]];
					if (k==AMBIG) {x+=K; slideSafe=0; continue;}
					//kmer = (kmer << pst_bshf_2 >> pst_bshf) | k;
					kmer = ((kmer << 2) | k) & maskK;
				} else {
					num = prefix_to_num(Raw+x,PL,&err,kpre_shf);
					if (err) {x+= num+1; slideSafe=0; continue;}
					kmer = dna_to_num(Raw+x+PL,SL,&err,kpst_shf);
					if (err) {x+= kmer+1; slideSafe=0; continue;}
					slideSafe = 1;
				}
				
				while (a < comp && !CONV[31 & Raw[x+a-comp]]) ++a; // compression
				if (a==comp) {
					uint64_t pod_ix;
					#pragma omp atomic capture
					pod_ix = KIx[num]++;
					
					KGrid[pod_ix] = (KPod){kmer,i};
				}
				++x;
			}
		}
		printf("Time: %f\n",omp_get_wtime()-wtime);
		uint64_t numK = 0;
		#pragma omp parallel for schedule(dynamic,1024) reduction(+:numK)
		for (uint64_t i = 0; i < nbins; ++i) {
			if (KIx[i] <= Nibs[i]) continue;
			KPod *start = KGrid + Nibs[i];
			size_t num = KIx[i]-Nibs[i];
			qsort(start,num,sizeof(*start),kpackcmp);
			numK += num;
		}
		printf("There were %lu k-mers.\n",numK);
		
		printf("Some stats on distributions!\n");
		uint64_t n_dupe = 0, n_multi = 0;
		#pragma omp parallel for schedule(dynamic,1024) reduction(+:n_dupe,n_multi)
		for (uint64_t i = 0; i < nbins; ++i) {
			if (KIx[i] <= Nibs[i]) continue;
			for (uint64_t j = Nibs[i]+1; j < KIx[i]; ++j) 
				if (KGrid[j].sfx==KGrid[j-1].sfx) {
					++n_multi;
					if (KGrid[j].rix==KGrid[j-1].rix) ++n_dupe;
				}
		}
		printf("Exact duplicates: %lu; across refs: %lu\n",n_dupe,n_multi);
		
		// Now write the results
		/* File structure:
		   1. Version byte and rix_t size [1]
		   2. Size of prefix [1]
		   3. Size of suffix [1]
		   4. Size of kmer_t [1]
		   5. Num refs [4]
		   6. Num kmers [8]
		   7. All prefix indices [1 << (2 * #2) x 8]
		   8. All kmer data [(#2 + #4) * #6]
		   9. Size of string data [8]
		   10. All strings [#9]
		   
		   11. Number of Ix1 in map [4] [0 means skip rest of file]
		   12. String size for Ix1 [8]
		   13. Ix1 strings dump [#12]
		   14. Number of Ix2 in map [4] [can be 0/skipped if no h2 map]
		   15. String size of Ix2 [8]
		   16. Ix2 strings dump [#15] 
		   //17. hpair_t dump [num ref by 8]
		   17. HPairs[0] dump [num ref by 4]
		   18. HPairs[1] dump [num ref by 4]
		*/
		
		uint64_t fileSz = 0, stringSz = 0;
		#pragma omp parallel for reduction(+:stringSz)
		for (uint32_t i = 0; i < ns; ++i) {
			uint64_t x = Offsets[i];
			uint64_t y = x; while (Raw[y] != '>') --y;
			stringSz += x-y - 1; // -1 for the '>' we're on, -1 '\n', but +1 '\0'
		}
		fileSz = 24 + sizeof(*Nibs)*(nbins+1) + sizeof(*KGrid)*numK + stringSz;
		printf("Initial file size = %lu\n",fileSz+4);
		
		FILE *db = fopen(dbPath,"wb");
		if (!db) {puts("I/O error: invalid db output file"); exit(1);}
		setvbuf(db, 0, _IOFBF, 1<<22);
		wtime = omp_get_wtime();
		fputc((VNO << 4) | sizeof(rix_t),db); 
		fputc(PL,db); fputc(SL,db); fputc(sizeof(kmer_t),db);
		size_t wrote = 4;
		wrote += fwrite(&ns,sizeof(ns),1,db); 
		wrote += fwrite(&numK,sizeof(numK),1,db);
		uint64_t tally = 0;
		for (uint64_t i = 0; i < nbins+1; ++i) {
			wrote += fwrite(&tally,sizeof(*Nibs),1,db);
			tally += KIx[i]-Nibs[i];
		}
		printf("First write: %f s.\n",omp_get_wtime()-wtime);
		wtime = omp_get_wtime();
		for (uint64_t i = 0; i < nbins; ++i) {
			if (KIx[i] <= Nibs[i]) continue;
			uint64_t num = KIx[i]-Nibs[i];
			fwrite(KGrid + Nibs[i], sizeof(*KGrid), num,db);
		}
		printf("Second write: %f s.\n",omp_get_wtime()-wtime);
		
		wtime = omp_get_wtime();
		fwrite(&stringSz,sizeof(stringSz),1,db);
		for (uint32_t i = 0; i < ns; ++i) {
			uint64_t x = Offsets[i];
			uint64_t y = x; while (Raw[y] != '>') --y;
			fwrite(Raw+y+1,1,x-y-2,db);fputc(0,db);
		}
		printf("Third write: %f s.\n",omp_get_wtime()-wtime);
		
		// TODO: also write statistics; for each genome, total K etc
		FILE *out = fopen(logPath,"wb");
		if (!out) printf("No log file specified; won't produce tally\n");
		else {
			uint32_t *TotK_m = huge_calloc((uint64_t)ns*sizeof(*TotK_m));
			uint32_t *TotUniq_m = huge_calloc((uint64_t)ns*sizeof(*TotUniq_m));
			#pragma omp parallel
			{
				int tid = omp_get_thread_num();
				uint32_t *TotK = TotK_m, *TotUniq = TotUniq_m;
				
				if (tid) 
					TotK = huge_calloc((uint64_t)ns*sizeof(*TotK)),
					TotUniq = huge_calloc((uint64_t)ns*sizeof(*TotUniq));
				
				#pragma omp for schedule(dynamic,1024)
				for (uint64_t i = 0; i < nbins; ++i) {
					if (KIx[i] <= Nibs[i]) continue; // empty
					uint32_t ambig = 0;
					uint64_t end = Nibs[i]+(KIx[i]-Nibs[i]), nd;
					kmer_t thisK = KGrid[Nibs[i]].sfx+1;
					
					for (uint64_t j = Nibs[i]; j < end; j += nd) {
						rix_t rix = KGrid[j].rix;
						
						// If new k-mer, check if ambig, store max value
						if (KGrid[j].sfx != thisK) {
							thisK = KGrid[j].sfx; 
							ambig = 0; 
							for (uint64_t k = j+1; k < end && KGrid[k].sfx == thisK; ++k)
								ambig |= KGrid[k].rix ^ rix;
						}
						
						// Find number of in-ref copies
						nd = 1;
						for (uint64_t k = j+1; k < end && 
						KGrid[k].sfx == thisK && KGrid[k].rix == rix; ++k) ++nd;
						
						// Increment the appropriate variables
						if (!ambig) TotUniq[rix] += nd;
						TotK[rix] += nd;
					}
				}
				#pragma omp critical
				if (tid) for (uint32_t i = 0; i < ns; ++i) 
						TotK_m[i] += TotK[i], TotUniq_m[i] += TotUniq[i];
			}
			fprintf(out,"Reference\tTotalKmers\tUniqKmers\n");
			for (uint64_t i = 0; i < ns; ++i) {
				uint64_t x = Offsets[i];
				uint64_t y = x; while (Raw[y] != '>') --y;
				fwrite(Raw+y+1,1,x-y-2,out); // write the ref name
				fprintf(out,"\t%u\t%u\n",TotK_m[i],TotUniq_m[i]);
			}
		}
		
		// Handle the H1/H2 mappings
		if (!ixPath) { // If no ixMap was provided, finish up.
			uint32_t zeroRef = 0;
			fwrite(&zeroRef,sizeof(zeroRef),1,db);
			exit(0); 
		}
		
		wtime = omp_get_wtime();
		
		// Read the mapping file in. Up to 3 columns are used (first is index)
		// proof of concept -- read whole file, sort all 3 (perhaps parallel) by ptrs
		// Then perform classical deduplication and data structure creation 
		free(Nibs); free(KGrid); free(KIx); 
		
		// New observation: can determine a priori number and placement of bins 
		// for all possible interpolations, given sorted map. Also search/store!
		
		//gzFile map = gzopen(ixPath,"rb");
		//int mapFsz = gzread(map,wholeMap,(uint64_t)1 << 38);
		FILE *map = fopen(ixPath,"rb");
		if (!map) {fprintf(stderr,"Can't open map: %s\n",ixPath); exit(2);}
		
		uint64_t sz = 0; 
		fseeko(map,0,SEEK_END); sz = ftello(map); rewind(map);
		if (sz < 2) {fprintf(stderr,"ERR: map malformatted\n"); exit(2);}
		char *wholeMap = huge_calloc(sz+16); 
		size_t mapFsz = fread(wholeMap,1,sz,map);
		if (mapFsz != sz) {puts("BAD MAP!"); exit(10101);}
		
		uint64_t nL = 0;
		if (wholeMap[sz-1]!='\n') ++nL;
		
		#pragma omp parallel for reduction(+:nL)
		for (uint64_t i = 0; i < sz; ++i) 
			nL += wholeMap[i] == '\n';
		
		printf("Map contained %lu lines.\n",nL);
		//typedef struct {char *str, *h1, *h2;} map_str_t;
		//map_str_t *MapStrs = calloc(nL,sizeof(*MapStrs));
		char **RefStr = calloc(nL,sizeof(*RefStr)),
			**H1Str = calloc(nL,sizeof(*H1Str)),
			**H2Str = calloc(nL,sizeof(*H2Str));
		char *map_ptr = wholeMap;
		int ncol = 2; uint64_t numTimesH2showed = 0;
		for (uint64_t i = 0; i < nL; ++i) {
			RefStr[i] = map_ptr;
			while (*map_ptr != '\t' && *map_ptr != '\n') ++map_ptr;
			if (*map_ptr == '\n') {puts("Bad map! Need >1 columns!"); exit(2);}
			*map_ptr++ = 0;
			H1Str[i] = map_ptr;
			while (*map_ptr != '\t' && *map_ptr != '\n') ++map_ptr;
			if (*map_ptr == '\n') {*map_ptr++ = 0; ncol = 1; /* printf("1col line %lu\n",i); */ continue;}
			*map_ptr++ = 0;
			H2Str[i] = map_ptr;
			++numTimesH2showed;
			while (*map_ptr != '\n') ++map_ptr;
			*map_ptr++ = 0;
		}
		wholeMap[sz-1] = 0; // in case it's a newline
		printf("Detected %d columns (a second showed up %lu times). Parsed. [%f]\n",ncol,numTimesH2showed,omp_get_wtime()-wtime);
		
		wtime = omp_get_wtime();
		uint64_t nuniq_ref = 0, nuniq_h1 = 0, nuniq_h2 = 0;
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				double wtimeL = omp_get_wtime();
				qsort(RefStr,nL,sizeof(*RefStr),strcmp_ptr);
				for (uint64_t i = 1; i < nL; ++i) // dedupe!
					if (strcmp(RefStr[i],RefStr[i-1])) 
						RefStr[nuniq_ref++] = RefStr[i-1];
				RefStr[nuniq_ref++] = RefStr[nL-1];
				printf("Name sort complete [%f]\n",omp_get_wtime()-wtimeL);
			}
			#pragma omp section
			{
				double wtimeL = omp_get_wtime();
				qsort(H1Str,nL,sizeof(*H1Str),strcmp_ptr);
				for (uint64_t i = 1; i < nL; ++i) // dedupe!
					if (strcmp(H1Str[i],H1Str[i-1])) 
						H1Str[nuniq_h1++] = H1Str[i-1];
				H1Str[nuniq_h1++] = H1Str[nL-1];
				printf("H1 sort complete [%f]\n",omp_get_wtime()-wtimeL);
			}
			#pragma omp section
			{
				double wtimeL = omp_get_wtime();
				if (ncol > 1) {
					qsort(H2Str,nL,sizeof(*H2Str),strcmp_ptr);
					for (uint64_t i = 1; i < nL; ++i) // dedupe!
						if (strcmp(H2Str[i],H2Str[i-1])) 
							H2Str[nuniq_h2++] = H2Str[i-1];
					H2Str[nuniq_h2++] = H2Str[nL-1];
				}
				printf("H2 sort complete [%f]\n",omp_get_wtime()-wtimeL);
			}
		}
		
		printf("All sorting complete. [%f]\n",omp_get_wtime()-wtime);
		printf("Unique: %lu refs, %lu H1's, %lu H2's.\n",
			nuniq_ref, nuniq_h1, nuniq_h2);
		wtime = omp_get_wtime();
		
		// Go thru each ref and match it to the list. 
		// Then lookup the h1 and h2 strings by proxy (past null).
		// Get the ids of the h1 and h2's from their respective lists.
		// Then add the pair [h1,h2] to the internal-ref-length struct.
		//hpair_t *HPairs = malloc(sizeof(*HPairs)*ns);
		uint32_t *HPairs[2] = { huge_calloc((uint64_t)ns*sizeof(*HPairs[0])),
			huge_calloc((uint64_t)ns*sizeof(*HPairs[1])) };
		#pragma omp parallel for schedule(dynamic,16)
		for (uint32_t i = 0; i < ns; ++i) {
			char *ref = Raw+Offsets[i];
			while (*ref != '>') --ref;
			uint64_t refmatch = binsearch_str_d(RefStr, ++ref, nuniq_ref);
			if (refmatch == (uint64_t)-1) {
				fprintf(stderr,"ERR: Map missing '%.*s'\n",
					(int)(Raw+Offsets[i] - ref - 1), ref); exit(2);
			}
			char *h = RefStr[refmatch];
			uint64_t h1match = 0, h2match = 0;
			while (*h) ++h;
			h1match = binsearch_str(H1Str,++h,nuniq_h1);
			if (h1match == (uint64_t)-1) {puts("INTERNAL ERROR H1"); exit(9);}
			if (ncol > 1) {
				while (*h) ++h;
				h2match = binsearch_str(H2Str,++h,nuniq_h2);
				if (h2match == (uint64_t)-1) {puts("INTERNAL ERROR H2"); exit(9);}
			}
			//HPairs[i] = (hpair_t){h1match,h2match};
			HPairs[0][i] = h1match, HPairs[1][i] = h2match;
		}
		
		// Should we make the taxonomy sublevel-map now or later? 
		// (In theory it can be made during read-in but...)
		// Update: now made during read-in!
		
		// Need to dump:
		// 1. String data (for both headers only; tracking length)
		// 2. The hpair_t dump -- update: now the two separate H arrays
		
		/*
		   11. Number of Ix1 in map [4] [0 means skip rest of file]
		   12. String size for Ix1 [8]
		   13. Ix1 strings dump [#12]
		   14. Number of Ix2 in map [4] [can be 0/skipped if no h2 map]
		   15. String size of Ix2 [8]
		   16. Ix2 strings dump [#15] 
		   //17. hpair_t dump [num ref by 8]
		   17. HPairs[0] dump [num ref by 4]
		   18. HPairs[1] dump [num ref by 4]
		*/
		
		
		stringSz = 0;
		#pragma omp parallel for reduction(+:stringSz)
		for (uint64_t i = 0; i < nuniq_h1; ++i) 
			stringSz += strlen(H1Str[i]) + 1;
		
		fileSz += 4 + 8 + stringSz;
		fwrite(&nuniq_h1,4,1,db); // # 11
		fwrite(&stringSz,sizeof(stringSz),1,db); // # 12
		for (uint64_t i = 0; i < nuniq_h1; ++i) 
			fwrite(H1Str[i],1,strlen(H1Str[i])+1,db); // #13
		
		fileSz += 4; // H2
		fwrite(&nuniq_h2,4,1,db); // #14
		stringSz = 0;
		if (ncol > 1) {
			#pragma omp parallel for reduction(+:stringSz)
			for (uint64_t i = 0; i < nuniq_h2; ++i) 
				stringSz += strlen(H2Str[i]) + 1;
			
			fwrite(&stringSz,sizeof(stringSz),1,db); // #15
			for (uint64_t i = 0; i < nuniq_h2; ++i) 
				fwrite(H2Str[i],1,strlen(H2Str[i])+1,db); // #16
		} else fwrite(&nuniq_h2,8,1,db); // 15 & 16 w/no col2
		fileSz += 8 + stringSz;
		
		fileSz += (uint64_t)ns*sizeof(*HPairs[0]);
		fwrite(HPairs[0],sizeof(*HPairs[0]),ns,db); // #17
		if (nuniq_h2) fileSz += (uint64_t)ns*sizeof(*HPairs[1]),
			fwrite(HPairs[1],sizeof(*HPairs[1]),ns,db); // #18
		
		printf("Final filesize: %lu\n",fileSz);
		
		exit(0);
	}
	
	/// Now handle the parsing and searching. 
	// TODO: make parser size-aware, not "num sequences" fixed
	
	if (ixPath) puts("WARNING: map file only applicable during DB BUILD");
	//  Read the database in.
	FILE *db = fopen(dbPath,"rb");
	if (!db) {puts("ERROR: bad input"); exit(2);}
	struct stat sb; int fno = fileno(db); fstat(fno,&sb);
	uint64_t fsz = sb.st_size, place = 0;
	char *Raw = mmap(0, fsz+16, PROT_READ, MAP_SHARED | 
		MAP_POPULATE, fno, 0);
	madvise(Raw,fsz,MADV_WILLNEED);
	double wtime = omp_get_wtime();
	if (doCopyMem) {
		puts("Copying database into local memory...");
		char *copied = huge_malloc(fsz+16);
		if (!copied) {puts("Can't do copymem on this system. Exiting..."); exit(3);}
		memcpy(copied,Raw,fsz+1);
		munmap(Raw,fsz+16);
		Raw = copied;
		printf("Copied into local memory [%f]\n",omp_get_wtime()-wtime);
		wtime = omp_get_wtime();
	}
	uint32_t ver = Raw[0] >> 4, rixSz = (uint8_t)Raw[0] & 15,
		PL = Raw[1], SL = Raw[2], ktSz = Raw[3];
	uint32_t numRef = *(uint32_t *)(Raw+4);
	uint64_t numK = *(uint64_t *)(Raw+8);
	printf("DBv: %d, rixSz = %d, PL = %d, SL = %d, ktSz = %d\n",
		ver, rixSz, PL, SL, ktSz);
	printf("Number of refs = %u, kmers = %lu\n",numRef, numK);
	if (sizeof(kmer_t)!=ktSz || sizeof(rix_t) != rixSz) 
		{puts("ERROR: wrong K or R size(s) for this DB"); exit(2);}
	place = 16; //#1-6
	
	// Initialize stats
	uint64_t K = PL+SL;
	uint32_t kpre_shf = PL*2-2;
	kmer_t kpst_shf = SL*2-2;
	uint32_t pre_bshf = 32-(PL*2), pre_bshf_2 = pre_bshf+2;
	kmer_t pst_bshf = sizeof(kmer_t)*8-(SL*2), pst_bshf_2 = pst_bshf+2;

	// Read the prefix array. 
	uint64_t *Nibs = (uint64_t *)(Raw+place);
	uint64_t nbins = (uint64_t)1 << (2*PL);
	place += (nbins+1) * sizeof(*Nibs); //#7
	//printf("DEBUG: nbins: %lu, Nibs[nbins] = %lu\n",nbins,Nibs[nbins]);
	KPod *KGrid = (KPod *)(Raw + place);
	printf("Size of kpod = %ld\n",sizeof(*KGrid));
	place += numK*sizeof(*KGrid); //#8
	
	// Read the ref names
	uint64_t stringSz = *(uint64_t *)(Raw + place);
	place += sizeof(stringSz); //#9
	char *RefRaw = Raw + place;
	place += stringSz; //#10
	char **RefNames = malloc(sizeof(*RefNames)*numRef); // defer
	//printf("String size = %lu\n",stringSz);
	//printf("String 1 = %s\n",RefRaw);
	
	// Read the h1 and h2 lists
	uint32_t nuniq_h1 = 0, nuniq_h2 = 0;
	char *H1Raw = 0, *H2Raw = 0, **HStr[2] = {0,0}; // **H1Str = 0, **H2Str = 0;
	//hpair_t *HPairs = 0;
	uint32_t *HPairs[2] = { 0,0 };
	nuniq_h1 = *(uint32_t *)(Raw+place);
	//printf("nuniq_h1 = %u\n",nuniq_h1);
	place += sizeof(nuniq_h1); //#11
	if (nuniq_h1) {
		HStr[0] = malloc((uint64_t)nuniq_h1*sizeof(*HStr[0]));
		stringSz = *(uint64_t *)(Raw + place);
		place += sizeof(stringSz); //#12
		H1Raw = Raw + place;
		place += stringSz; //#13
		
		nuniq_h2 = *(uint32_t *)(Raw+place);
		place += sizeof(nuniq_h2); //#14
		if (nuniq_h2) HStr[1] = malloc((uint64_t)nuniq_h2*sizeof(*HStr[1]));
		stringSz = *(uint64_t *)(Raw + place);
		place += sizeof(stringSz); //#15
		
		H2Raw = Raw + place;
		place += stringSz; //#16
		HPairs[0] = (uint32_t *)(Raw + place);
		place += (uint64_t)numRef * sizeof(*HPairs[0]); // #17
		if (nuniq_h2) HPairs[1] = (uint32_t *)(Raw + place);
		place += !nuniq_h2? 0 : (uint64_t)numRef * sizeof(*HPairs[1]); // #18
	}
	uint32_t NUniqH[2] = {nuniq_h1, nuniq_h2};
	printf("Read file of size = %lu (h1: %u, h2: %u)\n",place,nuniq_h1,nuniq_h2);
	
	// Make data structures to contain bin mappings
	// (first think how you want to count and create the interpolated
	// bins to be accessed -- ideally you'd specify an interpolation 
	// and get the unique (interpolated) index back.
	uint32_t **LBins[2] = {calloc(4096,sizeof(*LBins[0])), calloc(4096,sizeof(*LBins[1]))};
	//uint32_t **Tbins_h1 = calloc(4096,sizeof(*Tbins_h1));
	//uint32_t **Tbins_h2 = calloc(4096,sizeof(*Tbins_h2));
	
	// Parse the strings in parallel
	if (!nuniq_h1) 
		printf("WARNING: No taxonomy was included during DB formation\n");
	//FILE *debug = fopen("debug.txt","wb"); // DEBUG ONLY
	#pragma omp parallel sections
	{
		#pragma omp section
		{
		for (uint32_t i = 0; i < numRef; ++i) // Ref names
			RefNames[i] = RefRaw, RefRaw = strchr(RefRaw,0)+1;
		

		/* FILE *log = fopen("Dbg.txt","wb");
		for (uint32_t i = 0; i < numRef; ++i)
			fprintf(log,"%u\t%s\n",i,RefNames[i]);
		fflush(log); fclose(log); */
		//exit(10101);
		}

		#pragma omp section
		{
			//char **H1Str = HStr[0]; //uint32_t **Tbins_h1 = LBins[0];
			for (uint32_t i = 0; i < nuniq_h1; ++i) // H1 names
				HStr[0][i] = H1Raw, H1Raw = strchr(H1Raw,0)+1;
			for (uint32_t i = 0; i < nuniq_h1; ++i) {
				char *ref = HStr[0][i], *ptr = ref-1;
				int lv = 0;
				while (ptr = strchr(ptr+1,';')) {
					//fprintf(debug,"Searching: %s until %lu...\n --> %.*s\n",ref,ptr-ref,(int)(ptr-ref),ref);
					int64_t find = binsearch_str_L(HStr[0],ref,nuniq_h1,ptr-ref);
					//fprintf(debug,"%.*s\t%ld\n",(int)(ptr-ref),ref,find);
					if (!LBins[0][lv]) {
						LBins[0][lv] = calloc(nuniq_h1,sizeof(*LBins[0][lv]));
						for (uint32_t j = 0; j < nuniq_h1; ++j) LBins[0][lv][j]=-1;
					}
					LBins[0][lv++][i] = find;
				}
			}
		}
		#pragma omp section
		{
			for (uint32_t i = 0; i < nuniq_h2; ++i) // H2 names
				HStr[1][i] = H2Raw, H2Raw = strchr(H2Raw,0)+1;
			for (uint32_t i = 0; i < nuniq_h2; ++i) {
				char *ref = HStr[1][i], *ptr = ref-1;
				int lv = 0;
				while (ptr = strchr(ptr+1,';')) {
					int64_t find = binsearch_str_L(HStr[1],ref,nuniq_h2,ptr-ref);
					if (!LBins[1][lv]) {
						LBins[1][lv] = calloc(nuniq_h2,sizeof(*LBins[1][lv]));
						for (uint32_t j = 0; j < nuniq_h2; ++j) LBins[1][lv][j]=-1;
					}
					LBins[1][lv++][i] = find;
				}
			}
		}
	}
	
	// debug: print everything out.
	/* FILE *debug = fopen("debug.txt","wb");
	for (uint32_t i = 0; i < numRef; ++i) 
		fprintf(debug,"%s\t%s\t%s\n",RefNames[i],H1Str[HPairs[i].h1],H2Str[HPairs[i].h2]);
	fprintf(debug,"AND NOW THE H1 ACTION BOYS:\n");
	for (uint32_t i = 0; i < nuniq_h1; ++i) 
		fprintf(debug,"%s\t%u\t%u\n",H1Str[i],Tbins_h1[i+1],Tbins_h1[i+1]-Tbins_h1[i]);
	fprintf(debug,"AND NOW THE H2 ACTION BOYS:\n");
	if (nuniq_h2) for (uint32_t i = 0; i < nuniq_h2; ++i) 
		fprintf(debug,"%s\t%u\t%u\n",H2Str[i],Tbins_h2[i+1],Tbins_h2[i+1]-Tbins_h2[i]);
	
	// Find all the combos in each grid!!
	fprintf(debug,"AND NOW THE SEARCH IS ON!\n");
	for (uint32_t i = 0; i < nuniq_h1; ++i) {
		char *ref = H1Str[i], *ptr = ref;
		while (ptr = strchr(ptr,';')) {
			//printf("Searching: %s until %lu...\n",ref,ptr-ref-1);
			int64_t find = binsearch_str_L(H1Str,ref,nuniq_h1,ptr-ref-1);
			int nsemi = 0; for (char *p = ptr; *p; ++p) nsemi += *p==';';
			fprintf(debug,"%.*s\t%ld\t%ld\n",(int)(ptr-ref),ref,find,Tbins_h1[find+1]-nsemi);
			//printf("%.*s\t%ld\n",(int)(ptr-ref-1),ref,find);
			++ptr;
		}
		//int64_t find = binsearch_str(H1Str,ref,nuniq_h1);
		int64_t find = binsearch_str_L(H1Str,ref,nuniq_h1,9999);
		int64_t find2 = binsearch_str(H1Str,ref,nuniq_h1);
		if (find2 != find) {puts("ERROR FIND"); exit(6);}
		fprintf(debug,"%s\t%ld\t%ld\n",ref,find,Tbins_h1[find+1]);
	}
	
	exit(1); */
	
	#ifdef WSL
	uint16_t wsum = 0;
	for (uint64_t i = 0; i < numK; i+=4096/sizeof(*KGrid)) 
		wsum += KGrid[i];
	printf("4k Checksum: %d\n",wsum);
	#endif
	
	uint32_t *QueryAligns = calloc(numK,sizeof(*QueryAligns)),
		*FullQueryAligns = calloc((uint64_t)numRef,sizeof(*FullQueryAligns));
	
	
	// Open those queries up and align 'em
	KPod *EndPod = KGrid + Nibs[nbins] + 1;
	uint64_t n_raw = 0, n_filt = 0, n_matchedF = 0, n_matchedR = 0;
	uint32_t mask32 = ((uint64_t)1 << (2*PL)) - 1;
	kmer_t maskK = ((kmer_t)1 << (2*SL)) - 1;
	if (SL==sizeof(kmer_t)*4) maskK = -1;
	
	uint32_t preposit = PL*2-2, pstposit = SL*2-2;
	wtime = omp_get_wtime();
	printf("\n* Beginning alignment...\n");
	
	gzFile in;
	if (!strcmp(seqPath,"-")) in = gzdopen(fileno(stdin),"rb");
	else in = gzopen(seqPath,"rb");
	if (!in) {puts("ERROR: bad input fast[a/q][.gz]"); exit(2);}
	uint8_t *lineO = calloc(16 + (1 << 28),1), *line = lineO+16, 
		*head = malloc(1 << 20);
	uint32_t qChunk = 1<<16;
	uint8_t **QBucket = malloc(sizeof(*QBucket)*qChunk);
	char **HBucket = malloc(sizeof(*HBucket)*qChunk);
	*QBucket = calloc((uint64_t)131072*qChunk,sizeof(**QBucket));
	*HBucket = calloc((uint64_t)131072*qChunk,sizeof(**HBucket));
	
	// Make a bucket for temporary LCA & ref votes, per query
	uint64_t maxQsz = 1 << 27; //27; // arbitrary? Set a limit?
	uint64_t masterBnSz = (uint64_t)1 << 32;
	uint64_t masterLstSz = (uint64_t)1 << 31;

	typedef struct {uint32_t p; uint64_t s;} SBin_t;
	
	SBin_t **SBins = malloc(sizeof(*SBins)*threads);   // For ixing matches
	int32_t **RBins = malloc(sizeof(*RBins)*threads),  // For tally stores
		**TBins = malloc(sizeof(*TBins)*threads);
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		//printf("I'm thread %d / %d\n",tid, threads);
		int counter = 0;
		counter += !!(SBins[tid] = huge_malloc(maxQsz*sizeof(*SBins)));
		counter += !!(TBins[tid] = huge_malloc(maxQsz*sizeof(*TBins)));
		counter += !!(RBins[tid] = huge_calloc((uint64_t)numRef*sizeof(*RBins)));
		if (counter != 3) printf("Error allocating bins on thread %d (counter = %d)\n",tid,counter);
	}
	// For the capitalist bins (variable-length bins depending on hits). Thread-local.
	uint64_t **C_ixs = huge_malloc(sizeof(*C_ixs)*threads); //[3];   // for rix, h1, h2 cap arrays
	uint64_t **C_szs = huge_malloc(sizeof(*C_szs)*threads); //[3];   // the current size of the bins
	uint64_t ***C_bins = huge_malloc(sizeof(*C_bins)*threads); //[3]; // the bin arrays themselves
	//uint64_t C_MAX = (uint64_t)1 << 21;
	uint64_t C_INIT = (uint64_t)1 << 12;
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		C_ixs[tid] = calloc(3,sizeof(**C_ixs));
		C_szs[tid] = calloc(3,sizeof(**C_szs)); // enables thread resizing
		C_bins[tid] = calloc(3,sizeof(**C_bins));
		for (int i = 0; i < 3; ++i) 
			C_bins[tid][i] = malloc(C_INIT*sizeof(***C_bins)),
			C_szs[tid][i] = C_INIT-1;
	}
	// For "master" bins -- number of bins equal to number of queries. Single global storage.
	//typedef union {uint32_t a[4]; __uint128_t n; uint8_t b[16];} MasterBin_t;
	MasterBin_t *MasterBin = calloc(masterBnSz,sizeof(*MasterBin)); // bin IX // 32
	if (!MasterBin) {printf("ERROR: failed to allocate master bin!\n"); exit(3);}
	// a[0] = rix, a[1] = h1, a[2] = h2, a[3] = ...?
	uint8_t  **MasterIx = calloc(3,sizeof(*MasterIx)); // which thread's bin
	uint64_t **MasterList = calloc(3,sizeof(*MasterList)); 
	for (int i = 0; i < 3; ++i) {

		int counter = 0;
		counter += !!(MasterList[i] = calloc(masterLstSz,sizeof(*MasterList[i]))); //31
		counter += !!(MasterIx[i] = calloc(masterLstSz,sizeof(*MasterIx[i]))); //31
		if (counter != 2) printf("Error allocating master list (counter = %d)\n",counter);
	}

	FILE *outq = 0;
	if (perqPath) {
		outq = fopen(perqPath,"wb");
		if (!outq) {puts("ERROR: can't open per-q output file!"); exit(2);}
	}
	
	uint64_t nq, NQ = 0, nAligns = 0;

	while (nq = get_queries(in,QBucket,HBucket,head,line,qChunk,INT32_MAX)) {
		printf("Processed %lu queries\r",NQ); fflush(stdout);
		if (nq + NQ >= INT32_MAX) {puts("Exceeded 2B queries; stopping"); nq = 0; break;}
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			SBin_t *SBin = SBins[tid];
			int32_t *RBin = RBins[tid]; uint32_t *TBin = TBins[tid];
			uint64_t *C_ix = C_ixs[tid], *C_sz = C_szs[tid],
				**C_bin = C_bins[tid];
			#pragma omp for schedule(dynamic,1) reduction(+:n_raw,n_filt,n_matchedF,n_matchedR)
			for (uint64_t q = 0; q < nq; ++q) {
				
				uint8_t *Inf = QBucket[q];
				char *Qhed = HBucket[q];
				uint64_t x = 0, y = strlen(Inf);
				if (y >= maxQsz) 
					y = maxQsz-1,
					printf("\nQuery %lu exceed maximum length of %lu; truncating\n",q+NQ,maxQsz-1);
				int slideSafe = 0; 
				uint32_t num, numR, c, tix = 0;
				kmer_t kmer, kmerR, k;
				while (x + K <= y) {
					++n_raw;
					uint64_t slen, seed;
					if (slideSafe) {
						c = CONV[31 & Inf[x+PL-1]];
						if (c==AMBIG) {x+=PL; slideSafe=0; continue;}
						num = ((num << 2) | c) & mask32;
						k = CONV[31 & Inf[x+K-1]];
						if (k==AMBIG) {x+=K; slideSafe=0; continue;}
						kmer = ((kmer << 2) | k) & maskK;
						
					} else {
						uint32_t err = 0;
						num = prefix_to_num(Inf+x,PL,&err,kpre_shf);
						if (err) {x+= num+1; slideSafe=0; continue;}
						kmer = dna_to_num(Inf+x+PL,SL,&err,kpst_shf);
						if (err) {x+= kmer+1; slideSafe=0; continue;}
					}		
					
					slen = Nibs[num+1]-Nibs[num];
					if (!slen) goto RC_S;
					
					seed = LBS_k(KGrid+Nibs[num],slen,kmer);
					if (seed==(uint64_t)-1) goto RC_S;
					//if (seed>=slen || seed+Nibs[num] >= numK) goto RC_S;
					
					//#pragma omp atomic
					//QueryAligns[seed-KGrid]++; // TODO: write to temp first, then vote
					//if (tix >= maxQsz-1) {puts("\nWARNING: memory load high\n"); exit(0);}
					SBin[tix++] = (SBin_t){num,seed+Nibs[num]};
					
					++n_matchedF;
					RC_S:
					++n_filt;
					
					if (slideSafe) 
						c = RCONV[31 & Inf[x+K-1]],
						numR = (numR >> 2) | (c << preposit),
						k = RCONV[31 & Inf[x+SL-1]],
						kmerR = (kmerR >> 2) | (k << pstposit);
					else numR = prefix_to_num_RC(Inf+x+SL,PL),
						kmerR = dna_to_num_RC(Inf+x,SL);
					slen = Nibs[numR+1]-Nibs[numR];
					slideSafe = 1;
					
					if (!slen) {++x; continue;}
					seed = LBS_k(KGrid+Nibs[numR],slen,kmerR);
					if (seed==(uint64_t)-1) 
					//if (seed>=slen || seed+Nibs[numR] >= numK)
						{++x; continue;}
					
					//#pragma omp atomic
					//QueryAligns[seed-KGrid]++; // TODO: write to temp first, then vote
					SBin[tix++] = (SBin_t){numR,seed+Nibs[numR]};
					
					++n_matchedR; // bear in mind -- this is a raw k-mer tally.
					++x;
				}
				
				// Resize the query tally bins if needed (rix, h1, h2)
				if (doRedist) for (int j = 0; j < 3; ++j) if (tix+3 >= C_sz[j]-C_ix[j]) 
					//{printf("OOM: [%lu] C_ix[%d]\n",NQ+q,j); exit(3);}
					C_bin[j] = C_bins[tid][j] = realloc(C_bin[j],(C_sz[j]=2*(C_sz[j]+tix+3))*sizeof(*C_bin[j]));
						//madvise(C_bin[j], C_sz[j], MADV_HUGEPAGE);
				
				if (doRedist) for (int j = 0; j < 3; ++j) 
					MasterList[j][NQ+q] = C_ix[j],
					MasterIx[j][NQ+q] = tid;
				MasterBin[q+NQ].n = (__uint128_t)-1; // initialize MasterBin (ref, h1, h2, ?)
				if (!tix) {
					if (perqPath) fprintf(outq,"%s\tNo matches found\n",Qhed);
					if (doRedist) for (int j = 0; j < 3; ++j) 
						C_bin[j][C_ix[j]++] = -1;
					continue; // no alignments, no worries
				}
				
				/// Per-query processing (reporting, taxonomy tally, etc)
				// Vote for which reference to store
				int tempix = 0;
				// Go through and tally refs into storage bin RBin (init'ed to 0)
				for (uint32_t i = 0; i < tix; ++i) {
					uint32_t pfx = SBin[i].p;
					uint64_t ix = SBin[i].s;
					
					uint64_t hardstop = Nibs[pfx+1];
					rix_t prev_rix = -1; // set to non-this for dupe detection
					kmer_t prev_sfx = KGrid[ix].sfx;
					for (uint64_t j = ix; j < hardstop; ++j) {
						if (KGrid[j].sfx != prev_sfx) break;
						rix_t rix = KGrid[j].rix;
						//printf("Q %s: [%u/%u, piece %lu], %u / %s\n",Qhed,i,tix,j-ix,rix,RefNames[rix]);
						if (rix == prev_rix) continue; // don't double-count refs
						if (!RBin[rix]) TBin[tempix++] = rix;
						++RBin[rix];
						prev_rix = rix;
					}
				} // All unique references for this query are now in TBin by name, and RBin by count.

				// Resize the query tally bins if needed (rix, h1, h2)
				if (doRedist) for (int j = 0; j < 3; ++j) if (tempix+3 >= C_sz[j]-C_ix[j]) 
					//{printf("OOM: [%lu] C_ix[%d]\n",NQ+q,j); exit(3);}
					C_bin[j] = C_bins[tid][j] = realloc(C_bin[j],(C_sz[j]=2*(C_sz[j]+tempix+3))*sizeof(*C_bin[j]));
						//madvise(C_bin[j], C_sz[j], MADV_HUGEPAGE);
				
				// Consider all references that matched, uniquely
				uint32_t max = 0, max2 = 0, sum = 0;
				rix_t maxRix = -1, maxRix2 = -1;
				for (int i = 0; i < tempix; ++i) { // tally to determine max and 2nd max ref
					uint32_t rix = TBin[i];
					if (RBin[rix] > max || RBin[rix]==max && rix < maxRix) 
						max2 = max, maxRix2 = maxRix,
						max = RBin[rix], maxRix = rix;
					else if (RBin[rix] > max2)
						max2 = RBin[rix], maxRix2 = rix;
				}

				// Early terminate if no ref, or one ref got N or more assignments (controled by cmdline):
				if (!tempix || maxRix == -1 || max < nUniqMatches) {
					if (perqPath) fprintf(outq,"%s\tNo matches found\n",Qhed);
					if (doRedist) for (int j = 0; j < 3; ++j) 
						C_bin[j][C_ix[j]++] = -1;
					for (int i = 0; i < tempix; ++i) RBin[TBin[i]] = 0;
					continue; // no alignments, no worries
				}
				
				if (covPath) for (uint32_t i = 0; i < tix; ++i) { // coverage
					uint64_t ix = SBin[i].s;
					uint64_t hardstop = Nibs[SBin[i].p+1];
					kmer_t prev_sfx = KGrid[ix].sfx;
					for (uint64_t j = ix; j < hardstop; ++j) {
						if (KGrid[j].sfx != prev_sfx) break;
						rix_t rix = KGrid[j].rix;
						if (RBin[rix] == max) {
							#pragma omp atomic
							++QueryAligns[ix]; // TODO: change "ix" to "j" ?
						}
					}
				}
				for (int i = 0; i < tempix; ++i) { // ref-level set and clear
					uint32_t rix = TBin[i];
					if (RBin[rix] == max) {
						if (covPath) {
							#pragma omp atomic
							++FullQueryAligns[rix]; // # all tied refs
						}
						if (doRedist) C_bin[0][C_ix[0]++] = rix;
					}
					RBin[rix] = 0; // clear the bin
				}
				if (doRedist) C_bin[0][C_ix[0]++] = -1;
				
				uint32_t finalRix = maxRix;
				MasterBin[q+NQ].a[0] = finalRix;
				char *finalT1 = "", *finalT2 = "";
				//uint32_t f1_l = UINT16_MAX, f2_l = UINT16_MAX; // full length
				char *finalT[2] = {finalT1,finalT2};
				int32_t finalL[2] = {UINT16_MAX,UINT16_MAX};
				//if (!max2 || (!doFullLCA && max > max2 && (double)max/tix >= conf)) { // TODO: revisit this earlyterm (specifically conf)
				if (!max2 || (max > max2 && (double)max/tix >= conf)) { // TODO: revisit this earlyterm (specifically conf)
					// early terminate -- call this taxon and don't do LCA
					if (HStr[0]) finalT[0] = HStr[0][HPairs[0][maxRix]];
					if (taxPath) MasterBin[q+NQ].a[1] = HPairs[0][maxRix];
					if (HStr[1]) {
						finalT[1] = HStr[1][HPairs[1][maxRix]];
						if (taxPath) MasterBin[q+NQ].a[2] = HPairs[1][maxRix];
					}
					if (doRedist) {
						C_bin[1][C_ix[1]++] = HStr[0] ? HPairs[0][maxRix] : -1;
						C_bin[1][C_ix[1]++] = -1;
						if (HStr[1]) C_bin[2][C_ix[2]++] = HPairs[1][maxRix],
							C_bin[2][C_ix[2]++] = -1;
					}
				} else for (int H = 0; H < 2 && HStr[H]; ++H) { // Interpolate the taxonomies. 
					// Consider max level first!
					tempix = 0;
					// Go through and tally refs into storage bin RBin (init'ed to 0)
					for (uint32_t i = 0; i < tix; ++i) {
						uint32_t pfx = SBin[i].p;
						uint64_t ix = SBin[i].s;
						
						uint64_t hardstop = Nibs[pfx+1];
						kmer_t prev_sfx = KGrid[ix].sfx;
						for (uint64_t j = ix; j < hardstop; ++j) {
							if (KGrid[j].sfx != prev_sfx) break;
							rix_t rix = KGrid[j].rix;
							rix_t h = HPairs[H][rix];
							if (!RBin[h]) TBin[tempix++] = h;
							if (RBin[h] >= 0) RBin[h] = -RBin[h] - 1;
						}
						for (int j = 0; j < tempix; ++j)
							RBin[TBin[j]] = RBin[TBin[j]] < 0 ? -RBin[TBin[j]] : RBin[TBin[j]];
					}
					int32_t h_max1 = 0, h_maxIx1 = -1,
						h_max2 = 0, h_maxIx2 = -1;
					// Consider all references that matched, uniquely
					for (int i = 0; i < tempix; ++i) { // tally to determine max and 2nd max ref
						uint32_t rix = TBin[i];
						if (RBin[rix] > h_max1 || RBin[rix]==h_max1 && rix < h_maxIx1) 
							h_max2 = h_max1, h_maxIx2 = h_maxIx1,
							h_max1 = RBin[rix], h_maxIx1 = rix;
						else if (RBin[rix] > h_max2)
							h_max2 = RBin[rix], h_maxIx2 = rix;
						//RBin[rix] = 0; // clear the bin
					}
					for (int i = 0; i < tempix; ++i) { // add best ref to redistribute pile; clear Ref count hash (RBin)
						uint32_t rix = TBin[i];
						if (doRedist && RBin[rix] == h_max1) // TODO: revisit adding a specificity control in here (i.e. only if !h_max2 or h_max1/h_max2 > 1.5...)
							C_bin[H+1][C_ix[H+1]++] = rix;
						RBin[rix] = 0; // clear the bin
					}
					if (doRedist) C_bin[H+1][C_ix[H+1]++] = -1;
					// Report if good....
					//if (perqPath) fprintf(outq,"-->%s\t[%u]\t%s\t%u\t%s\t%u\n",Qhed,tix,H1Str[h1_maxIx1],h1_max1,H1Str[h1_maxIx2],h1_max2);
						
					if (!h_max2 || (!doFullLCA && h_max1 > h_max2 && (double)h_max1/tix >= conf))  // TODO: revisit this early terminate (particularly conf)
					//if (!h_max2 || (h_max1 > h_max2 && (double)h_max1/tix >= conf))  // TODO: revisit this early terminate (particularly conf)
						finalT[H] = HStr[H][h_maxIx1];
					else { // We are in for pain. Do a full aufbau
						uint32_t agreed = tix; // Start out all agreeing at 0 semicolons!
						uint32_t ag_thres = conf*tix; 
						uint32_t winner = -1, winLv = -1;
						// TODO: implement blank-aware aufbau (remove blank "" taxa from consideration and total)
						//printf("%u matches in query %s\n",tix,Qhed);
						for (int semi = 1; agreed >= ag_thres; ++semi) {
							if (!LBins[H][semi-1]) break;
							agreed = 0; 
							int tempix = 0;
							for (uint32_t i = 0; i < tix; ++i) {
								uint64_t hardstop = Nibs[SBin[i].p+1];
								kmer_t prev_sfx = KGrid[SBin[i].s].sfx;
								for (uint64_t j = SBin[i].s; j < hardstop; ++j) {
									if (KGrid[j].sfx != prev_sfx) break;
									uint32_t h = HPairs[H][KGrid[j].rix];
									char *ref = HStr[H][h];
									if (LBins[H][semi-1][h]==-1) continue;
									uint32_t find = LBins[H][semi-1][h];
									//printf("[%u:%lu] Searching [lv %d]: %.*s\n --> Found: %ld: %s\n",
									//	i,j,semi,(int)(ptr-ref),ref,find,H1Str[find]);
									if (!RBin[find]) TBin[tempix++] = find;
									if (RBin[find] >= 0) RBin[find] = -RBin[find] - 1;
								}
								for (int j = 0; j < tempix; ++j)
									RBin[TBin[j]] = RBin[TBin[j]] < 0 ? -RBin[TBin[j]] : RBin[TBin[j]];
							}
							// The winner will be the count that is >= ag_thres.
							uint32_t local_max = 0, local_max2 = 0;
							uint32_t local_winner = 0;
							for (int i = 0; i < tempix; ++i) { // go thru unique options 
								uint32_t thisTax = TBin[i];
								if (RBin[thisTax] >= ag_thres) {
									if (RBin[thisTax] > local_max)
										local_max2 = local_max,
										local_max = RBin[thisTax],
										local_winner = thisTax;
									else if (RBin[thisTax] > local_max2) 
										local_max2 = RBin[thisTax];
								}
								RBin[thisTax] = 0; // also reset RBin
							}
							if (local_max > local_max2 && local_max >= ag_thres) 
								agreed = local_max, winner=local_winner, winLv = semi;
							//printf("Purity at lv %d: %u/%u [2nd: %u] (%s)\n --> Winner %d @ lv %d [%s]\n",
							//	semi,local_max,tix,local_max2,local_winner? H1Str[local_winner] : "[NONE]",winner,winLv,H1Str[winner]);
						}
						
						// All levels have been explored. Now let's expand the winner's string. 
						if (winner != -1) {
							
							if (perqPath) {
								int L = 0; char *p = HStr[H][winner]-1;
								for (int s = 0; s < winLv; ++s) 
									p = strchr(p+1,';');
								L = p - HStr[H][winner];
								finalL[H] = L;
								finalT[H] = HStr[H][winner];
							}
							if (taxPath) MasterBin[q+NQ].a[H+1] = winner + winLv * NUniqH[H];
						}
					}
				}
				
				DONE_TAXA:NULL;
				#pragma omp atomic
				++nAligns;
				if (perqPath) fprintf(outq,"%s\t%s\t[%u,%u]\t%.*s\t%.*s\t%u\n",Qhed,
					finalRix!=-1? RefNames[finalRix] : "",max,max2,finalL[0],finalT[0],
					finalL[1],finalT[1],tix);
			}
		}
		NQ += nq;
	}
	gzclose(in);
	printf(" - Total k-mers: %lu, non-error: %lu, matchedF: %lu, matchedR: %lu\n",
		n_raw, n_filt, n_matchedF, n_matchedR);
	printf("--> Successfully aligned %lu/%lu queries [%f].\n",nAligns,NQ,omp_get_wtime()-wtime);
	
	
	if (doRedist) {
		wtime = omp_get_wtime();
		printf("\n* Performing capitalist redistribution...\n");
		// The goal is to go multi-pass through each cap bin
		uint64_t Sizes[3] = {(refPath? numRef: 0),NUniqH[0],NUniqH[1]};
		int Pass[3] = {0};
		for (int i = 0; i < 3; ++i) {
			if (!Sizes[i]) continue;
			uint64_t *List = MasterList[i];
			uint8_t *Ix = MasterIx[i];
			uint64_t *Tally = huge_calloc(Sizes[i]*sizeof(*Tally)),
				*NextTally = huge_calloc(Sizes[i]*sizeof(*NextTally));
			if (!Tally || !NextTally) {puts("ERROR: OOM TALLY"); exit(3);}
			uint64_t debugSum_q = 0;
			for (uint64_t q = 0; q < NQ; ++q) {
				//if (MasterBin[q].a[0]==-1) continue;
				uint64_t *Bin = C_bins[Ix[q]][i]+List[q];
				debugSum_q += *Bin != -1;
				while (*Bin != -1) ++Tally[*Bin++];
			}
			//for (uint64_t j = 0; j < Sizes[i]; ++j) NextTally[j] = Tally[j];

			// Skip over blanks (disable with "--redistribute-blanks")
			// First, figure out which (if any) taxa are blanks. If multiple, stop at first.
			//printf("Beginning assignment, header %d\n",i); fflush(stdout);
			uint64_t firstIx = -1, next = 0;
			
			if (i > 0) {
				#pragma omp parallel for reduction (min:firstIx)
				for (uint64_t j = 0; j < Sizes[i]; ++j) 
					if (!HStr[i-1][j][0]) firstIx = firstIx < j ? firstIx : j;
				next = firstIx + 1; if (next >= Sizes[i]) next = 0;
				//printf("--> Found blank taxonomy string at H %d entry ix %lu (next: %s)\n",
				//	i,firstIx,HStr[i-1][next]); fflush(stdout);
			}
			
			uint64_t debugSum = 0;
			for (uint64_t j = 0; j < Sizes[i]; ++j) debugSum += Tally[j];
			//printf("Sum %lu [qsum %lu]; Doing passes now...\n",debugSum,debugSum_q); fflush(stdin); 
			uint64_t changes = -1, conv = NQ/100000, maxPass = 100;
			if (doFastRedist) maxPass = 1;
			for (int pass = 0; pass < maxPass && changes > conv; ++pass) {
				changes = 0;
				for (uint64_t q = 0; q < NQ; ++q) {
					//if (MasterBin[q].a[0]==-1) continue;
					uint64_t *Bin = C_bins[Ix[q]][i]+List[q];
					uint64_t maxTally = 0, whichMax = -1;
					while (*Bin != -1) {
						//if (Tally[*Bin] > maxTally) 
						if (Tally[*Bin] > maxTally && (*Bin != firstIx || whichMax == -1)) 
							maxTally = Tally[*Bin], whichMax = *Bin;
						++Bin;
					}
					
					if (whichMax != -1) ++NextTally[whichMax];
				}
				for (uint64_t j = 0; j < Sizes[i]; ++j) 
					changes += abs((int64_t)Tally[j] - (int64_t)NextTally[j]);
				//printf(" - pass: %d, changes: %lu\n",p, changes);
				
				uint64_t *t = Tally; Tally = NextTally; NextTally = t;
				memset(NextTally,0,sizeof(*NextTally)*Sizes[i]);
				++Pass[i];
			}
			
			//printf("Assigning distros...\n"); fflush(stdout);
			
			// Now assign final distributions
			for (uint64_t q = 0; q < NQ; ++q) {
				uint64_t *Bin = C_bins[Ix[q]][i]+List[q];
				uint64_t maxTally = 0, whichMax = -1;
				while (*Bin != -1) {
					if (Tally[*Bin] > maxTally && (*Bin != firstIx || whichMax == -1)) 
						maxTally = Tally[*Bin], whichMax = *Bin;
					++Bin;
				}
				MasterBin[q].a[i] = whichMax;
			}
			free(Tally); free(NextTally);
		}
		printf("--> Capitalist redistribution complete (%d,%d,%d) [%f]\n",
			Pass[0],Pass[1],Pass[2],omp_get_wtime()-wtime);
	}
	
	if (refPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing reference tally table...\n");
		FILE *taxout = fopen(refPath,"wb");
		if (!taxout) {printf("ERROR: couldn't write to %s!\n",refPath); exit(2);}
		
		uint64_t *TaxTally = huge_calloc((uint64_t)numRef*sizeof(*TaxTally));
		for (uint64_t i = 0; i < NQ; ++i) // multithread?
			if (MasterBin[i].a[0] != -1) ++TaxTally[MasterBin[i].a[0]];
		for (uint32_t i = 0; i < numRef; ++i) if (TaxTally[i]) 
			fprintf(taxout,"%s\t%lu\n",RefNames[i],TaxTally[i]);
		free(TaxTally);
		fclose(taxout);
		printf("--> Reference table written [%f]\n",omp_get_wtime()-wtime);
	}
	
	if (taxPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing taxonomy tally table...\n");
		FILE *taxout = fopen(taxPath,"wb");
		if (!taxout) {printf("ERROR: couldn't write to %s!\n",taxPath); exit(2);}
		// One option is to identify the max value for H1 and H2's values in MasterBin,
		// and allocate accordingly
		uint32_t maxH[2] = {0};
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[1] != -1 && MasterBin[q].a[1] > maxH[0])
				maxH[0] = MasterBin[q].a[1];
			if (MasterBin[q].a[2] != -1 && MasterBin[q].a[2] > maxH[1])
				maxH[1] = MasterBin[q].a[2];
		}
		//printf("Maximum tid for H1: %u, H2: %u\n",maxH[0],maxH[1]);
		
		// Indices >= NUniqH signify interpolation at lv = i/NUniqH
		
		// Create a bucket to tally the taxa
		for (int H = 0; H < 2; ++H) {
			if (!maxH[H]) continue;
			//printf("Number of unique headers to consider: %u, maxH = %u\n",NUniqH[H],maxH[H]);
			uint64_t *TaxTally = huge_calloc((uint64_t)(maxH[H]+1)*sizeof(*TaxTally));
			for (uint64_t i = 0; i < NQ; ++i) // multithread?
				if (MasterBin[i].a[H+1] != -1) ++TaxTally[MasterBin[i].a[H+1]];
			uint32_t lv = 0, nextLv = NUniqH[H];
			for (uint32_t i = 0; i <= maxH[H]; ++i) if (TaxTally[i]) { // 
				while (i >= nextLv) // Control depth of taxonomy reporting
					nextLv += NUniqH[H], ++lv;
				if (!lv) fprintf(taxout,"%s\t%lu\n",HStr[H][i],TaxTally[i]);
				else { // the interpolation case
					char *s = HStr[H][i-(nextLv-NUniqH[H])], *orig = s;
					for (uint32_t semi = 0; semi < lv; ++s)
						semi += *s==';';
					fprintf(taxout,"%.*s\t%lu\n",(uint32_t)(s-orig-1),orig,TaxTally[i]);
				}
			}
			free(TaxTally);
		}
		fclose(taxout);
		printf("--> Taxonomy written [%f]\n",omp_get_wtime()-wtime);
	} // end taxonomy tally output
	
	if (orthogPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing orthogonal tally table...\n");
		if (!NUniqH[0] || !NUniqH[1]) 
			{printf("--> ERROR: Orthogonalizing requires 2 taxonomies.\n"); exit(2);}
		FILE *taxout = fopen(orthogPath,"wb");
		if (!taxout) {printf("ERROR: couldn't write to %s!\n",orthogPath); exit(2);}
		
		uint32_t maxH[2] = {0};
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[1] != -1 && MasterBin[q].a[1] > maxH[0])
				maxH[0] = MasterBin[q].a[1];
			if (MasterBin[q].a[2] != -1 && MasterBin[q].a[2] > maxH[1])
				maxH[1] = MasterBin[q].a[2];
		}
		
		//#define PRIME 1299827
		#define PRIME 4969
		uint64_t HashTable[PRIME+2] = {0};
		
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[0]==-1) continue;
			uint64_t val = *(uint64_t *)(MasterBin[q].b+4);
			++HashTable[(val % PRIME)+2];
		}
		for (uint64_t i = 1; i <= PRIME+1; ++i) HashTable[i] += HashTable[i-1];
		printf(" - Number added: %lu\n",HashTable[PRIME+1]);
		MasterBin_t **SortedIX = huge_calloc(HashTable[PRIME+1]*sizeof(*SortedIX));
		for (uint64_t q = 0; q < NQ; ++q) {
			if (MasterBin[q].a[0]==-1) continue;
			uint64_t val = *(uint64_t *)(MasterBin[q].b+4);
			SortedIX[HashTable[(val % PRIME)+1]++] = MasterBin + q;
		}
			
		for (uint64_t h = 0; h < PRIME; ++h) {
			uint64_t range = HashTable[h+1]-HashTable[h];
			MasterBin_t **ThisBinP = SortedIX+HashTable[h];
			if (range) qsort(ThisBinP,range,sizeof(*ThisBinP),binCmp);
		}
		// print result
		for (uint64_t h = 0; h < PRIME; ++h) {
			uint64_t range = HashTable[h+1]-HashTable[h];
			MasterBin_t **ThisBinP = SortedIX+HashTable[h];
			if (range) {
				uint64_t last = *(uint64_t *)(ThisBinP[0]->b+4);
				uint64_t tally = 0;
				for (uint64_t i = 0; i < range; ++i) {
					uint64_t val = *(uint64_t *)(ThisBinP[i]->b+4);
					if (val != last || i == range-1) { // commit
						uint32_t h1 = ThisBinP[i-1]->a[1], 
							h2 = ThisBinP[i-1]->a[2],
							lv1=h1/NUniqH[0], lv2=h2/NUniqH[1],
							L1 = 0, L2 = 0;
						for (int s = 0; s < lv1; ++L1) 
							s += HStr[0][h1][L1]==';';
						for (int s = 0; s < lv2; ++L2) 
							s += HStr[1][h2][L2]==';';
						L1 = L1 ?: UINT16_MAX;
						L2 = L2 ?: UINT16_MAX;
						fprintf(taxout,"%.*s\t%.*s\t%lu\n",L1-1,HStr[0][h1],L2-1,HStr[1][h2],tally);
						tally = 0;
					}
					++tally; 
					last = val; 
				}
			}
		}
		
		fclose(taxout);
		printf("--> Taxonomy written [%f]\n",omp_get_wtime()-wtime);
		
	}
	
	if (covPath) {
		wtime = omp_get_wtime();
		printf("\n* Printing coverage table...\n");
		uint32_t *TotK_m = huge_calloc((uint64_t)numRef*sizeof(*TotK_m));
		uint32_t *TotUniq_m = huge_calloc((uint64_t)numRef*sizeof(*TotUniq_m));
		uint64_t *FoundK_m = huge_calloc((uint64_t)numRef*sizeof(*FoundK_m));
		uint64_t *FoundUniq_m = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq_m));
		uint32_t *PropK_m = huge_calloc((uint64_t)numRef*sizeof(*PropK_m));
		uint32_t *PropUniq_m = huge_calloc((uint64_t)numRef*sizeof(*PropK_m));
		
		#pragma omp parallel
		{
			int tid = omp_get_thread_num();
			uint32_t *TotK = TotK_m, *TotUniq = TotUniq_m, 
				*PropK=PropK_m, *PropUniq=PropUniq_m;
			uint64_t *FoundK = FoundK_m, *FoundUniq = FoundUniq_m;
			
			if (tid) 
				TotK = huge_calloc((uint64_t)numRef*sizeof(*TotK)),
				TotUniq = huge_calloc((uint64_t)numRef*sizeof(*TotUniq)),
				FoundK = huge_calloc((uint64_t)numRef*sizeof(*FoundK)),
				FoundUniq = huge_calloc((uint64_t)numRef*sizeof(*FoundUniq)),
				PropK = huge_calloc((uint64_t)numRef*sizeof(*PropK)),
				PropUniq = huge_calloc((uint64_t)numRef*sizeof(*PropK));
			
			#pragma omp for schedule(dynamic,1024)
			for (uint64_t i = 0; i < nbins; ++i) {
				if (Nibs[i+1] <= Nibs[i]) continue; // empty
				uint32_t ambig = 0;
				uint64_t end = Nibs[i+1], nd;
				uint32_t mv = 0;
				kmer_t thisK = KGrid[Nibs[i]].sfx+1;
				
				for (uint64_t j = Nibs[i]; j < end; j += nd) {
					rix_t rix = KGrid[j].rix;
					if (rix >= numRef) {
						printf("ERROR: rix at bin %lu, nib %lu = %lu, >= numRef [%lu]. nbins = %lu\n",
							i, j, rix, numRef,nbins);
						continue;
					}
					
					// If new k-mer, check if ambig, store max value
					if (KGrid[j].sfx != thisK) {
						thisK = KGrid[j].sfx; 
						ambig = 0; mv = QueryAligns[j];
						for (uint64_t k = j+1; k < end && KGrid[k].sfx == thisK; ++k)
							mv = mv > QueryAligns[k] ? mv : QueryAligns[k],
							ambig |= KGrid[k].rix ^ rix;
					}
					
					// Find number of in-ref copies
					nd = 1;
					for (uint64_t k = j+1; k < end && 
					KGrid[k].sfx == thisK && KGrid[k].rix == rix; ++k) ++nd;
					
					// Increment the appropriate variables
					if (!ambig) {
						TotUniq[rix] += nd;
						FoundUniq[rix] += mv;
						PropUniq[rix] += mv < nd ? mv : nd;
					}
					TotK[rix] += nd;
					FoundK[rix] += mv;
					PropK[rix] += mv < nd ? mv : nd;
				}
			}
			#pragma omp critical
			if (tid) {
				for (uint32_t i = 0; i < numRef; ++i) {
					TotK_m[i] += TotK[i], TotUniq_m[i] += TotUniq[i],
					FoundK_m[i] += FoundK[i], FoundUniq_m[i] += FoundUniq[i],
					PropK_m[i] += PropK[i], PropUniq_m[i] += PropUniq[i];
				}
			}
		}
		
		printf(" - Done with coverage loop. [%f]\n",omp_get_wtime() - wtime);
		printf(" - Beginning accumulation...\n");
		
		uint64_t totFoundKUniq = TotUniq_m[numRef-1];
		for (uint64_t i = 1; i < numRef; ++i) totFoundKUniq+=TotUniq_m[i-1];
		printf(" - Total k-mers found: %lu (%lu unique)\n",numK,totFoundKUniq);
		
		//wtime = omp_get_wtime();
		
		FILE *out = fopen(covPath,"wb");
		if (!out) {puts("ERROR: can't open output file!"); exit(2);}
		setvbuf(out, 0, _IOFBF, 1<<22);
		fprintf(out, "Reference\tKmers_found\t");
		fprintf(out, "Unique_kmers_found\tKmers_covered\tUnique_kmers_covered\t");
		fprintf(out, "Proportion_covered\tUnique_proportion_covered\tReads_covered\n");
		//FILE *log = fopen("DbgPost.txt","wb");
		for (uint32_t i = 0; i < numRef; ++i) {
			//fprintf(log,"%u\t%s\n",i,RefNames[i]);

			if (!FoundK_m[i]) continue;
			fprintf(out,"%s\t%lu\t%lu\t%u\t%u",RefNames[i],//TotK_m[i],TotUniq_m[i],
				FoundK_m[i],FoundUniq_m[i],PropK_m[i],PropUniq_m[i]);
			fprintf(out,"\t%f\t%f\t%u\n",(double)PropK_m[i]/TotK_m[i],
				(double)PropUniq_m[i]/TotUniq_m[i],FullQueryAligns[i]); 
		}
		printf("--> Coverage table written [%f]\n",omp_get_wtime()-wtime);
		fclose(out);
	} // end coverage output
	
	printf("\nAll outputs written.\n");
	return 0;
}
