/* UTree -- ultrafast unique k-mer mapper by Gabe. 
Copyright 2015-2017 Knights Lab, Regents of the University of Minnesota.
This software is released under the GNU Affero General Public License (AGPL) v3.0.
*/
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h> 
#include <fcntl.h>
#include <math.h>
#include <omp.h>

off_t lseek(int fd, off_t offset, int whence); 
/* D options:
 -D BUILD, BUILD_GG, SEARCH, SEARCH_GG, COMPRESS
 NO_COUNT
 IXTYPE=uint16_t, CNTTYPE=uint32_t, PACKSIZE=32[4,8,16,32,64]
 PFBITS=16[24 max for laptops, 26-28 max tb servers]
 SLACK=3, SPARSITY=4


*/
/*  RadixalTriecrobium. A tree utility for k-mer manipulation.
 Gutted the API because nobody (would) use(d) it
*/

#ifdef _OPENMP
#include <omp.h>
#endif
#define NO_COUNT
#ifndef IXTYPE
	#define IXTYPE uint16_t
#endif
#ifndef CNTTYPE
	#define CNTTYPE uint32_t
#endif
#ifndef PACKSIZE
	#define PACKSIZE 32
#endif
// Hash of trees
#ifndef PFBITS
	#define PFBITS 24
#endif
#if PFBITS>16
	#define PFTYPE uint32_t
#else
	#define PFTYPE uint16_t
#endif
#define KHASH_SIZE (WTYPE)((WTYPE)1 << PFBITS)
#define BALBUFSZ 65536
#if PACKSIZE==64
	#define WTYPE __uint128_t
	#define STYPE WTYPE
#elif PACKSIZE==32
	#define WTYPE uint64_t
	#define STYPE WTYPE
#elif PACKSIZE==16
	#define WTYPE uint32_t
	#define STYPE uint16_t
#elif PACKSIZE==8
	#define WTYPE uint16_t
	#define STYPE uint8_t
#elif PACKSIZE==4
	#define WTYPE uint8_t
	#define STYPE WTYPE
#endif
const char* TYPEARR[17] = {"NA","uint8_t", "uint16_t","NA","uint32_t",
"NA","NA","NA","uint64_t","NA","NA","NA","NA","NA","NA","NA","__uint128_t"};

unsigned critical_cutoff = 2;

typedef struct KMerY KMerY;
struct 
/* __attribute__ ((__packed__))  */
KMerY {
	STYPE word;
	IXTYPE ix;
	KMerY *left, *right;
};

typedef struct SampX SampX;
struct /* __attribute__ ((__packed__))  */SampX {
	char *name;
	IXTYPE ix;
	SampX *left, *right;
};

// Transformation of DNA to numeral
WTYPE *C2Xb; _Bool INIT_K = 0;
char *X2C = "ACGTNNNNNNNNNNNNNNNN";
char *X2C_RC = "TGCANNNNNNNNNNNNNNNN";

size_t UBAL_THRES=7, FIRETHRES = 1000000; 

// Comment out the below to disable memanage
KMerY ***UBANK = 0; 
size_t UBANK_MAXK = 10000, UBANK_INITBINS = 100;
size_t *UBANK_BIN = 0, *UBANK_BINCNT = 0, *UBANK_IX = 0;

char WORDTEMP[PACKSIZE+1] = {0};
IXTYPE BAD_IX = (IXTYPE)-1;
IXTYPE EMPTY_IX = (IXTYPE)-2;
KMerY BAD_NODE; 

// initialize DNA converters
void initConverter() {
	if (!INIT_K) {
	C2Xb = malloc(256*sizeof(WTYPE));
	for (int i = 0; i < 256; ++i) C2Xb[i] = 255;
	C2Xb['a'] = 0; C2Xb['A'] = 0; 
	C2Xb['c'] = 1; C2Xb['C'] = 1; 
	C2Xb['g'] = 2; C2Xb['G'] = 2;
	C2Xb['t'] = 3; C2Xb['T'] = 3;
	BAD_NODE = (KMerY){0,BAD_IX,0,0};
	INIT_K = 1;
	}
}
typedef struct /* __attribute__ ((__packed__))  */{
	WTYPE word;
	IXTYPE ix;
} WordIxPair;

typedef struct {
	KMerY **Roots; // prefix cache
	size_t *BalanceThreshes, *NumsInserted, *TotalCounts;  // parallel to prefix cache
	size_t numInserted, totalCount; // global
	int numThreads; // only for searching
	SampX *Samps; // sample tree
	IXTYPE sampIX; // current IX (num of samples currently in)
	char **SampStrings; // array of sample names
	uint64_t *SampCnts; // array of sample counts
	uint8_t *semicolons;	
	size_t sampStringSz; // size of above array
	size_t queuedClumps;
	WordIxPair *Pairs; // Word and IX in here
	uint64_t *BinIx;
	char *Dump;
} UTree;

// Comment in the below and out the one after to disable memanage
/*inline KMerY * u_xalloc(STYPE word, IXTYPE ix) {
	KMerY *new = malloc(sizeof(*new));
	*new = (KMerY){word,ix,0,0};
	return new;
} */

inline KMerY * u_xalloc(STYPE word, IXTYPE ix) { // KBANK,KBANK_INITBINS,KBANK_MAXK,KBANK_BIN,KBANK_IX
	#ifdef MINRAM
		#define BUMP_MEM 1.1
	#else
		#define BUMP_MEM 1.25
	#endif 
	KMerY *Kptr = UBANK[0][UBANK_BIN[0]] + UBANK_IX[0];
	#ifdef NO_COUNT
		*Kptr = (KMerY){word,ix,0,0};
	#else
		*Kptr = (KMerY){word,1,ix,0,0};
	#endif
	if (++UBANK_IX[0] == UBANK_MAXK) { // reset the ix, increment bin
		UBANK_IX[0] = 0;
		if (++UBANK_BIN[0] == UBANK_BINCNT[0]) { // resize bin array
			UBANK[0] = realloc(UBANK[0],
				sizeof(*UBANK[0])*(UBANK_BINCNT[0]*=BUMP_MEM));
			if (!UBANK[0]) { puts("ERROR: xalloc 1"); exit(3); }
			for (size_t x=UBANK_BINCNT[0]/BUMP_MEM; x<UBANK_BINCNT[0]; ++x) {
				UBANK[0][x] = malloc(UBANK_MAXK*sizeof(*UBANK[0][x]));
				if (!UBANK[0][x]) { puts("ERROR: xalloc 2"); exit(3); }
			}
		}
	}
	return Kptr;
}


#define EXPAND_SAMPARRAYS() \
if (++ktree->sampIX >= ktree->sampStringSz) { \
	ktree->SampStrings = realloc(ktree->SampStrings, \
		sizeof(*ktree->SampStrings)*(ktree->sampStringSz *=2)); \
	DO_SAMPCNTS() \
}
#define EXPAND_SAMPCNTS() \
ktree->SampCnts = realloc(ktree->SampCnts, \
	sizeof(*ktree->SampCnts)*(ktree->sampStringSz)); \
for (long i=ktree->sampStringSz/2; i < ktree->sampStringSz; ++i) \
	ktree->SampCnts[i] = 0;

#define ADDSAMP() \
SampX *tree = ktree->Samps; \
int cmp = strcmp(str,tree->name); \
do { \
	if (cmp > 0) {\
		if (!tree->right) { \
			tree->right = malloc(sizeof(*tree->right)); \
			tree->right->name = malloc(strlen(str)+1); \
			strcpy(tree->right->name,str); \
			tree->right->right = 0, tree->right->left = 0; \
			EXPAND_SAMPARRAYS() \
			ktree->SampStrings[ktree->sampIX] = tree->right->name; \
			return tree->right->ix = ktree->sampIX; \
		} \
		tree = tree->right; \
	} \
	else if (cmp < 0) { \
		if (!tree->left) { \
			tree->left = malloc(sizeof(*tree->left)); \
			tree->left->name = malloc(strlen(str)+1); \
			strcpy(tree->left->name,str); \
			tree->left->right = 0, tree->left->left = 0; \
			EXPAND_SAMPARRAYS() \
			ktree->SampStrings[ktree->sampIX] = tree->left->name; \
			return tree->left->ix = ktree->sampIX; \
		} \
		tree = tree->left; \
	} \
} while (cmp=strcmp(str,tree->name)); \
return tree->ix;

#define DO_SAMPCNTS() {}
inline IXTYPE addSampleUd(UTree *ktree, char *str) { ADDSAMP() }
inline IXTYPE addSampleU(UTree *ktree, char *str) {
	if (!ktree->Samps) {
		ktree->Samps = malloc(sizeof(*ktree->Samps));
		char *strT = malloc(strlen(str)+1);
		strcpy(strT,str);
		*ktree->Samps = (SampX){strT,0,0,0};
		ktree->SampStrings = malloc(sizeof(*ktree->SampStrings)*ktree->sampStringSz);
		//ktree->SampCnts = calloc(ktree->sampStringSz,sizeof(*ktree->SampCnts));
		*ktree->SampStrings = strT;
		return 0;
	}
	ADDSAMP()
}
#undef DO_SAMPCNTS
#define DO_SAMPCNTS EXPAND_SAMPCNTS
static IXTYPE addSampleUdX(UTree *ktree, char *str) { ADDSAMP() }

// returns node status: -1 = just went bad, 0 = no change, 1 = just created
inline int xeTreeU(KMerY *tree, WTYPE word, IXTYPE ix) { 
	do {
		if (word > tree->word) { // go right
			if (!tree->right) {
				tree->right = u_xalloc(word,ix);
				return 1;
			}
			tree = tree->right; 
		}
		else if (word < tree->word) { // go left
			if (!tree->left) {
				tree->left = u_xalloc(word,ix);
				return 1;
			}
			tree = tree->left;
		}
	} while (word != tree->word);
	
	if (ix != tree->ix) { 
		if (tree->ix == BAD_IX) return 0; // already bad
		tree->ix = BAD_IX; 
		return -1; 
	}
	return 0;
}
// returns node status as above; 2+ = converted to new ix
static long xeTreeU_RF(KMerY *tree, WTYPE word, IXTYPE ix, UTree *ut) { 
	do {
		if (word > tree->word) { // go right
			if (!tree->right) {
				tree->right = u_xalloc(word,ix);
				return 1;
			}
			tree = tree->right; 
		}
		else if (word < tree->word) { // go left
			if (!tree->left) {
				tree->left = u_xalloc(word,ix);
				return 1;
			}
			tree = tree->left;
		}
	} while (word != tree->word);
	if (ix != tree->ix) {
		if (tree->ix >= EMPTY_IX) return 0; // already bad
		char *old = ut->SampStrings[tree->ix], *oldOrig = old;
		char *new = ut->SampStrings[ix];
		
		unsigned numP = 0, ixP = 0;
		while (*old == *new) {
			if (*old == ';') ++numP, ixP = old - oldOrig;
			++old, ++new;
		}
		if (numP < critical_cutoff) { tree->ix = BAD_IX; return -1; }
		char buffer[UINT16_MAX+1];
		memcpy(buffer,oldOrig,ixP);
		buffer[ixP] = 0;
		IXTYPE new_ix = addSampleUd(ut,buffer);
		
		//printf("--> New IX = %u [old=%u]\n",new_ix,tree->ix);
		tree->ix = new_ix;
		//exit(3);
		return 2+(long)new_ix;
	}
	return 0;
}

// Modify tree and set sparsity pointer based on whether
// a k-mer conflicts with one already in the tree.

inline int parapluieU(UTree *ktree, WTYPE wordO, IXTYPE ix) {
	int z = 1, denom = 1;
	PFTYPE prefix = wordO >> ((PACKSIZE << 1) - PFBITS);
	KMerY *tree = ktree->Roots[prefix];
	STYPE word = (wordO << PFBITS) >> PFBITS;
	if (tree->ix == EMPTY_IX) {
		*tree = (KMerY){word,ix,0,0};
		ktree->NumsInserted[prefix] = 1;
		ktree->TotalCounts[prefix] = 1;
		return z/denom;
	}
	do { // PARAPLUIE
		if (word > tree->word) { 
			if (!tree->right) { 
				tree->right = u_xalloc(word,ix); 
				++ktree->TotalCounts[prefix]; 
				++ktree->NumsInserted[prefix]; 
				return z/denom; 
			} 
			tree = tree->right; 
		} 
		else if (word < tree->word) { 
			if (!tree->left) { 
				tree->left = u_xalloc(word,ix); 
				++ktree->TotalCounts[prefix]; 
				++ktree->NumsInserted[prefix]; 
				return z/denom; 
			} 
			tree = tree->left; 
		} 
	} while (word != tree->word); 
	if (ix == tree->ix) return 0;
	tree->ix = BAD_IX; 
	++ktree->TotalCounts[prefix]; 
	return z/denom; /*2*z/(denom+1); */
}

void traceBalanceU(KMerY *tree, KMerY **array, size_t *ix) {
	if (tree->left) traceBalanceU(tree->left, array, ix);
	array[(*ix)++] = tree; // if on top, DFS. If mid, IOS, if bot: LFS
	if (tree->right) traceBalanceU(tree->right, array, ix);
}
void traceBalancePurgeU(KMerY *tree, KMerY **array, size_t *ix, size_t *cnt) {
	if (tree->left) traceBalancePurgeU(tree->left, array, ix, cnt);
	if (tree->ix < EMPTY_IX) array[(*ix)++] = tree; 
	if (tree->right) traceBalancePurgeU(tree->right, array, ix,cnt);
}


void buildBalanceLU(KMerY *tree, KMerY **array, size_t sz);
void buildBalanceRU(KMerY *tree, KMerY **array, size_t sz);
#define BUILDBALANCEU() \
	if (!sz) { \
		CHILD = *array; \
		CHILD->left = 0; \
		CHILD->right = 0; \
		return; \
	} \
	size_t ix = sz >> 1; \
	CHILD = array[ix]; \
	if (ix) buildBalanceLU(CHILD,array,ix-1); \
	else CHILD->left = 0; \
	buildBalanceRU(CHILD,array+(ix+1), sz-(ix+1));

// set a branch of the given tree, and recurse with that branch as root
void buildBalanceLU(KMerY *tree, KMerY **array, size_t sz) {
	#define CHILD tree->left
		BUILDBALANCEU()
	#undef CHILD
}
void buildBalanceRU(KMerY *tree, KMerY **array, size_t sz) {
	#define CHILD tree->right
		BUILDBALANCEU()
	#undef CHILD
}

char * print_bignum(WTYPE n, char *str) {
	//char str[40] = {0}; // log10(1 << 128) + '\0'
	if (n == 0) return str;
	char *s = str + 39;
	while (n != 0) {
		*--s = "0123456789"[n % 10];
		n /= 10;
	}
	return s;
}

void traceTreeBU(PFTYPE i, KMerY *tree, FILE *of) {
	if (tree->left) traceTreeBU(i, tree->left, of);
	if (tree->ix < EMPTY_IX) {
		WTYPE word = ((WTYPE)i << (2*PACKSIZE - PFBITS)) + tree->word;
		fwrite(&word,sizeof(word),1,of);
		fwrite(&tree->ix,sizeof(IXTYPE),1,of);
	}	
	if (tree->right) traceTreeBU(i, tree->right, of);
}
void traceTreeBUf(PFTYPE i, KMerY *tree, uint64_t *cnts, FILE *of) {
	if (tree->left) traceTreeBUf(i, tree->left, cnts, of);
	if (tree->ix < EMPTY_IX) {
		WTYPE word = ((WTYPE)i << (2*PACKSIZE - PFBITS)) + tree->word;
		++cnts[tree->ix];
		fwrite(&word,sizeof(word),1,of);
		fwrite(&tree->ix,sizeof(IXTYPE),1,of);
	}	
	if (tree->right) traceTreeBUf(i, tree->right, cnts, of);
}

inline KMerY *simpleBalanceU(KMerY *tree, size_t sz) {
	//if (sz < 3) return tree;
	size_t ix = 0;
	//KMerY **array = malloc(sizeof(*array) * sz);
	KMerY *array[sz];
	traceBalanceU(tree, array, &ix);
	ix = --sz >> 1;
	tree = array[ix];
	buildBalanceLU(tree, array, ix - 1);
	buildBalanceRU(tree, array + (ix+1), sz - (ix+1));
	//free(array);
	return tree;
}

/************************************************
          USER-ACCESSIBLE FUNCTIONS
************************************************/

void UT_addWordIx(UTree *utree, WTYPE word, IXTYPE ix) {
	PFTYPE prefix = word >> ((PACKSIZE << 1) - PFBITS); 
	STYPE suffix = (word << PFBITS) >> PFBITS; 
	++utree->TotalCounts[prefix];
	if (prefix >= KHASH_SIZE) {puts("ERROR"); exit(5);}
	//printf("prefix=%llu (%llu)\n",(uint64_t)prefix, (uint64_t)KHASH_SIZE);
	if (utree->Roots[prefix]->ix == EMPTY_IX) {  
		*utree->Roots[prefix] = 
		(KMerY){suffix,ix,0,0}, utree->NumsInserted[prefix]=1; 
		return; 
	} 
	if (xeTreeU(utree->Roots[prefix],suffix,ix)==1) ++utree->NumsInserted[prefix]; 
	if (utree->NumsInserted[prefix] == utree->BalanceThreshes[prefix]) { 
		utree->Roots[prefix] = simpleBalanceU(utree->Roots[prefix], 
			utree->NumsInserted[prefix]); 
		utree->BalanceThreshes[prefix] = utree->BalanceThreshes[prefix] > 1000 ? 0 : 
			(utree->BalanceThreshes[prefix]+1)*2-1; 
	} 
}

void UT_addWordIxRF(UTree *utree, WTYPE word, IXTYPE ix) {
	PFTYPE prefix = word >> ((PACKSIZE << 1) - PFBITS); 
	STYPE suffix = (word << PFBITS) >> PFBITS; 
	++utree->TotalCounts[prefix];
	if (utree->Roots[prefix]->ix == EMPTY_IX) {  
		*utree->Roots[prefix] = 
		(KMerY){suffix,ix,0,0}, utree->NumsInserted[prefix]=1; 
		return; 
	} 
	if (xeTreeU_RF(utree->Roots[prefix],suffix,ix,utree)==1) ++utree->NumsInserted[prefix]; 
	if (utree->NumsInserted[prefix] == utree->BalanceThreshes[prefix]) { 
		utree->Roots[prefix] = simpleBalanceU(utree->Roots[prefix], 
			utree->NumsInserted[prefix]); 
		utree->BalanceThreshes[prefix] = utree->BalanceThreshes[prefix] > 1000 ? 0 : 
			(utree->BalanceThreshes[prefix]+1)*2-1; 
	} 
}

inline size_t crBST(char *key, size_t sz, char **String) {
	char **p = String;
	while (sz) {
		size_t w = sz >> 1; 
		char *ref_s = *(p+w+1), *key_s = key;
		
		while (*ref_s == *key_s++) if (!*ref_s++) return p+w+1-String; 
		if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	char *ref_s = *p, *key_s = key;
	while (*ref_s == *key_s++) if (!*ref_s++) return p - String;
	return -1;
	//return p - String; // replace last 3 lines for unsave ver
}

inline int xcmp(str1, str2) register const char *str1, *str2; {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}

int xcmpP(const void *str1, const void *str2) {
	return xcmp(**(char ***)str1, **(char ***)str2);
}


size_t UT_parseSampFastaExternOSFA(UTree *utree, char* filename, char* labels, 
int ixCol, int lblCol, uint32_t lv, int doGG) {
	FILE *fp = fopen(filename, "rb"), *fp2 = fopen(labels, "rb");
	if (fp == 0 || fp2 == 0) { puts("Invalid input file(s)"); exit(1); }
	// read in the second file
	fseek(fp2,0,SEEK_END); size_t filesz = ftell(fp2); rewind(fp2);
	char *dump = malloc(filesz+1); if (!dump) {puts("Map err."); exit(2);}
	fread(dump,filesz,1,fp2); fclose(fp2);
	dump[filesz] = 0;
	printf("Parsed map. %llu bytes",filesz);
	if (!filesz) {puts("\nInput map empty."); exit(1);}
	size_t lines = 0;
	char *ptr = dump; do if (*ptr == '\n') ++lines; while (*++ptr);
	if (*(ptr-1)!='\n') ++lines;
	printf(", %llu lines.\n",lines);
	//printf("Taxonomy labels (raw) = %llu\n",lines); 
	char **lblList = malloc(sizeof(*lblList)*lines);
	char **ixList = malloc(sizeof(*ixList)*lines);
	if (!lblList || !ixList) {puts("Map list err."); exit(2);}
	//FILE *outfile = fopen("outfile_temp_db.txt","wb");
	ptr = dump; 
	void (*addWord)(UTree *, WTYPE, IXTYPE) = doGG ? &UT_addWordIxRF : &UT_addWordIx; 
	if (ixCol < lblCol) for (size_t i = 0; i < lines; ++i) { //default
		int j = 0;
		for (; j < ixCol; ++j) {
			//puts("Called non-default");
			while (*++ptr != '\t') if (!*ptr) {printf("Err tab0: %llu\n",i); exit(2);}; 
			*ptr++ = 0; // skip to correct tab & nullify
		}
		if (*ptr == '\n' || *ptr == '\r') { 
			printf("ERROR: map line %llu\nBlank indices are NOT ALLOWED.\n",i);
			exit(2);
		}
		ixList[i] = ptr; // store first string as key
		if (!ptr) {printf("Stored null in ixList: %llu\n",i); exit(2);}
		for (; j < lblCol; ++j) {
			if (*ptr == '\t') {printf("map: extra tab, line %llu\n",i); exit(2);}
			while (*++ptr != '\t') if (!*ptr) {printf("Err tab1: %llu\n",i); exit(2);}; 
			*ptr++ = 0; // nullify after tab reached
		}
		if (*ptr == '\n' || *ptr == '\r') { 
			printf("\nERROR: map line %llu\nBlank labels are NOT ALLOWED.\n",i+1);
			exit(2);
		}
		lblList[i] = ptr; // store second string as taxon
		if (!ptr) {printf("Stored null in lblList: %llu\n",i); exit(2);}
		while (*ptr != '\n') { //} && *ptr != '\t') {
			if (!*ptr) {printf("Err line counter: %llu\n",i); exit(2);}
			if (*ptr == '\r' || *ptr == '\t') *ptr = 0;
			ptr++;
		}
		if (!*ptr) {printf("Malformatted map before line %llu",i+1); exit(2);}
		*ptr++=0;
	}
	else for (size_t i = 0; i < lines; ++i) { // TODO: update this?
		int j=0; for (; j < lblCol; ++j) {
			while (*++ptr != '\t'); *ptr++ = 0; // skip to correct tab & nullify
		}
		lblList[i] = ptr; // store first string as key
		for (; j < ixCol; ++j) {
			while (*++ptr != '\t'); *ptr++ = 0; // skip to correct tab & nullify
		}
		ixList[i] = ptr; // store second string as taxon
		while (*++ptr != '\n') if (*ptr=='\t') *ptr = 0; // get to end of line
		*ptr++ = 0; // replace newline with null
	}
	// sorting by pointers is done here
	char ***Pointers = malloc(sizeof(*Pointers)*lines);
	if (!Pointers) {puts("Map ptr err."); exit(2);}
	for (size_t i = 0; i < lines; ++i) Pointers[i] = &ixList[i];
	qsort(Pointers,lines,sizeof(*Pointers),xcmpP);
	char **ixSorted = malloc(sizeof(*ixSorted)*lines),
		 **lblSorted = malloc(sizeof(*lblSorted)*lines);
	if (!ixSorted || !lblSorted) {puts("Map srtd err."); exit(2);}
	for (size_t i = 0; i < lines; ++i) {
		ixSorted[i] = *Pointers[i];
		lblSorted[i] = lblList[Pointers[i] - ixList];
	}
	//free(Pointers), free(ixList), free(lblList);
	--lines; // lowerbound on search function is max index, not total no
	size_t ns=0, LINELEN = 268435456; // 256MB lines
	char *line = malloc(LINELEN + 1), *origLine = line;
	if (!line) {puts("FASTA parse: memory error."); exit(2);}
	const uint32_t k1 = PACKSIZE - 1; uint32_t kv = k1 + lv;
	while (++ns, line = fgets(line,LINELEN,fp)) { 
		char *src = line + 1; // sample name parsed to generate ix
		src = strchr(src,'\n'); 
		if (src) memset(src,'\0',1); 
		size_t pre_ix = crBST(line+1,lines, ixSorted); 
		//printf("Line %llu, preIX=%llu \n",ns,pre_ix);
		//fprintf(outfile,"Tag is %s, which maps to %llu\n",line+1,pre_ix);
		if (pre_ix == -1) {printf("Error: taxon map incomplete (line %u)\n",ns); exit(4);}
		IXTYPE ix = addSampleU(utree,lblSorted[pre_ix]);
		//fprintf(outfile,"Tag is %s, which maps to %s, which is ix %u\n",line+1,lblSorted[pre_ix],ix);
		if (!(line = fgets(line,LINELEN,fp))) // encode sequence
			{ printf("Error parsing FASTA (1pass): %llu",ns); exit(2); }
		src = line;
		register uint32_t length = strlen(src);
		if (src[length-1] == '\n') --length; // lop off newline(s)
		if (src[length-1] == '\r') --length; // supports every platform!
		WTYPE w = 0;
		// add clumps and indices to tree
		for (uint32_t i = kv; i < length; ++i) {
			//WTYPE m = *(C2Xb+src[i]);
			if (lv >= 1) {
				if (C2Xb[src[i-kv]] != 0) continue;
				if (lv >= 2) {
					if (C2Xb[src[i-kv+1]] != 2) continue;
					if (lv >= 3) {
						if (C2Xb[src[i-kv+2]] != 1) continue;
						if (lv >= 4) {
							if (C2Xb[src[i-kv+3]] != 3) continue;
						}
					}
				}
			}
			w = 0;
			for (uint32_t j = i - k1, p = j; j <= i; ++j) {
				if (C2Xb[src[j]] == 255) {i += j - p + lv; goto ENDER;}
				w <<= 2u, w |= C2Xb[src[j]];
			}
			addWord(utree, w, ix);
			ENDER:NULL;
		}
	}
	size_t totNodes = 0;
	#pragma omp parallel for reduction(+:totNodes)
	for (size_t i = 0; i < KHASH_SIZE; ++i) 
		totNodes += utree->NumsInserted[i];
	printf("Done with sequence parse: %llu k-mers made\n",totNodes);
	if (!totNodes) {puts("Error: no k-mers. Bad input/params!"); exit(2); }
	fclose(fp); free(origLine); // purge buffers
	free(dump);
	return ns-1; // number of sequences parsed
}

typedef struct WordCountPair WordCountPair;
struct WordCountPair {
	char *word;
	uint16_t count;
	WordCountPair *left, *right;
};
WordCountPair * WCP_insert(WordCountPair *tree, WordCountPair **wcp, char *string, uint16_t count) {
	if (!tree->count) {
		*tree = (WordCountPair){string,count,0,0};
		return tree;
	}
	int cmp = strcmp(string,tree->word);
	do {
		if (cmp > 0) {
			if (!tree->right) {
				tree->right = ++*wcp;
				*tree->right = (WordCountPair){string,count,0,0};
				return tree->right;
			}
			tree = tree->right;
		}
		else if (cmp < 0) {
			if (!tree->left) {
				tree->left = ++*wcp;
				*tree->left = (WordCountPair){string,count,0,0};
				return tree->left;
			}
			tree = tree->left;
		}
	} while (cmp = strcmp(string,tree->word));
	tree->count += count;
	return tree;
}
// Code-name "Humpty_dumpty". Compile with -Wall


// Define search decompression macros 
#define ADDR_AT(raw,i) ((raw) + SZ*(i))
#define WORD_AT(raw,i,pre) ((*(WTYPE*)ADDR_AT(raw,i) & MASK) | pre)
#define SUFFIX_AT(raw,i) (*(WTYPE*)ADDR_AT(raw,i) & MASK)
#define IX_AT(raw,i) *(IXTYPE*)(ADDR_AT(raw,i)+CMPWDSZ)

#define FIRST_WORD(raw,pre) ((*(WTYPE*)(raw) & MASK) | pre)
#define FIRST_SUFFIX(raw) (*(WTYPE*)(raw) & MASK)
#define FIRST_IX(raw) *(IXTYPE*)((raw)+CMPWDSZ)
#define IXDIST(raw,new) ((new)-(raw))/SZ
#define PREFIX(word) (word & PMASK)
#define PREFIX_L(word) (word >> SXBITS)
#define SUFFIX(word) (word & MASK)
#define PREFIX_AT32(ix,BinIx) ((WTYPE)uWBS32(BinIx,ix,0,BINRANGE) << SXBITS)

// Define search compression constants
#define CMP 3
#define PXBITS (CMP * 8)
const int SZ= sizeof(WTYPE) + sizeof(IXTYPE) - CMP, 
	CMPWDSZ = sizeof(WTYPE) - CMP;
const unsigned int NUMBINS = (1 << PXBITS) + 1;
const int SXBITS = PACKSIZE*2 - PXBITS;
const size_t BINRANGE = 1 << PXBITS;

WTYPE MASK = 0, PMASK = 0;

inline char * xtSuffixBS( char array[], size_t size, WTYPE sx) {
	char *p=array;
	while (size) {
		size_t w = size >> 1;
		if (SUFFIX_AT(p,w+1) <= sx) p+=(w+1)*SZ, size-=w+1; 
		else size =w;
	}
	return FIRST_SUFFIX(p)==sx ? p : 0;
}

inline uint32_t uWBS32(uint64_t *ixList, uint64_t key, uint64_t low, uint64_t high) {
	uint64_t middle; //, low = 0, high = range;
	while (low < high) {
		middle = low + ((high - low) >> 1);
		if (key > ixList[middle]) low = middle + 1;
		else high = middle;
	}
	if (ixList[low] > key) --low;
	return low; 
}

inline IXTYPE XT_getIX32(UTree *utree, WTYPE word) {
	uint64_t *BinIx = utree->BinIx;
	WTYPE qprefix = PREFIX_L(word), qsuff = SUFFIX(word);
	//printf("Prefix=%u, suffix = %u\n",(unsigned)qprefix,(unsigned)qsuff);
	size_t start_i = BinIx[qprefix], end_i = BinIx[qprefix+1];
	//printf("--> Start = %llu, end = %llu\n",start_i, end_i);
	if (start_i >= end_i) return BAD_IX;
	
	char *found = xtSuffixBS(ADDR_AT(utree->Dump,start_i), end_i-start_i-1, qsuff);
	return found ? FIRST_IX(found) : BAD_IX;
}

int readSamplesFPdelim(UTree *ktree, FILE *fp, char delim);
UTree *XT_read32(char *db, char delim) {
	FILE *dp = fopen(db,"rb");
	if (!dp) { puts("Invalid DB file"); exit(0); }
	uint64_t metadata[4] = {0};
	size_t numRead = fread(metadata, sizeof(*metadata), 4, dp);
	if (numRead < 4 || !metadata[3]) {puts("Tree malformatted."); exit(0); }
	// Check file for compatibility with this compilation config
	#ifdef NO_COUNT
		#define CNT_SIZE 0
	#else
		#define CNT_SIZE sizeof(CNTTYPE)
	#endif
	uint64_t numNodes = metadata[3]; 
	if (metadata[0] != sizeof(WTYPE) || metadata[1] != CNT_SIZE ||
		metadata[2] != sizeof(IXTYPE)) {
			printf("ERROR. Input tree requires PACKSIZE=%u, CNTTYPE=%s, IXTYPE=%s\n", 
				metadata[0] << 2,TYPEARR[metadata[1]], TYPEARR[metadata[2]]);
		exit(0);
	}
	// Prepare the bin delimiters 
	// See constants declared above
	if (numNodes < UINT32_MAX) puts("Using 32-bit counters");
	else puts("Holey smokes, a tree of over 4 billion k-mers. Here goes...");
	uint64_t *BinIx = malloc(NUMBINS*sizeof(*BinIx));
	int ixSize = numNodes < UINT32_MAX ? sizeof(uint32_t) : sizeof(uint64_t);
	numRead = 0;
	for (size_t i = 0; i < NUMBINS; ++i) numRead += fread(&BinIx[i],ixSize,1,dp);
	//numRead = fread(BinIx, sizeof *BinIx,NUMBINS,dp);
	printf("%llu elements read.\n",numRead);
	
	// Read in the tree
	printf("Nodes in input tree: %llu (PACKSIZE=%u, CNTTYPE=%s, IXTYPE=%s, SZ=%d)\n",
		numNodes, metadata[0] << 2, TYPEARR[metadata[1]], TYPEARR[metadata[2]],SZ);
	char *Dump = calloc(numNodes * SZ + 32,1); // add extra WTYPE buffer on end
	numRead = fread(Dump, SZ, numNodes, dp);
	if (numRead != numNodes) {puts("Error in reading tree."); exit(3);}
	printf("Read %llu nodes.\n",numRead);
	
	// Prepare annotation
	UTree *utree = calloc(1,sizeof *utree);
	utree->sampStringSz = 10;
	
	int annotated = readSamplesFPdelim(utree, dp, delim);
	if (!annotated) puts("No annotation found in tree file."); 
	fclose(dp);
	
	// ATCGGAAAGAATCCCT
	char *cmask = malloc(sizeof(WTYPE)); //, *cpmask = malloc(sizeof(WTYPE));
	for (int i = 0; i < CMPWDSZ; ++i) cmask[i] = 0xff;
	for (int i = CMPWDSZ; i < sizeof(WTYPE); ++i) cmask[i] = 0x00;
	MASK = *(WTYPE*)cmask, PMASK = ~MASK; // register ? const ?
	
	// Save in tree object
	initConverter();
	utree->Dump = Dump;
	utree->BinIx = BinIx;
	utree->totalCount = numNodes;
	
	// Check consistency (node total matches summed sub-trees, max index matches max taxonomy string
	if (BinIx[NUMBINS-1] != numNodes) 
		printf("Warning: detected nodes %u != %u\n",BinIx[NUMBINS-1],numNodes);
	#ifdef DEBUG
	printf("Initiating full debug consistency checks...\n");
	uint32_t maxIX = 0, outOfBounds = 0, real = utree->sampIX+1; 
	#pragma omp parallel for reduction(max:maxIX) reduction(+:outOfBounds)
	for (WTYPE i = 0; i < numNodes; ++i) {
		//XT_getIX32(utree, i)
		IXTYPE ix = IX_AT(Dump,i);
		if (ix > maxIX) maxIX = ix;
		if (annotated && ix >= real) ++outOfBounds;
	}
	if (outOfBounds) printf("Warning: Maximum index = [%u / %u[+1]] (%u out of bounds)\n",maxIX,real,outOfBounds);
	maxIX = 0, outOfBounds = 0;
	uint64_t outOfRange = 0, nuts = 0;
	#pragma omp parallel for reduction(max:maxIX) reduction(+:outOfBounds,outOfRange,nuts)
	for (uint64_t i = 0; i < NUMBINS; ++i) {
		if (BinIx[i+1] <= BinIx[i]) continue;
		nuts += BinIx[i+1] - BinIx[i];
		for (uint64_t j = BinIx[i]; j < BinIx[i+1]; ++j) {
			
			IXTYPE ix = IX_AT(Dump,j);
			if (j >= numNodes) puts("Dang.");
			if (ix > maxIX) maxIX = ix;
			//if (PREFIX_AT32(ix,BinIx) != i) ++outOfRange;
			if (annotated && ix >= real) ++outOfBounds;
			
		}
	}
	printf("nuts = %llu\n",nuts);
	if (outOfBounds) printf("Warning R2: Maximum index = [%u / %u[+1]] (%u out of bounds)\n",maxIX,real,outOfBounds);
	if (outOfRange) printf("Warning: There were %u nodes with reported ix out of search range\n",outOfRange);
	#endif
	
	puts("Tree read.");
	return utree;
}

typedef struct {char *s; uint32_t n;} String_Count_t;
static int byStr(const void *A, const void *B) {
	return strcmp(((String_Count_t *)A)->s, ((String_Count_t *)B)->s);}
static inline size_t XT_doSearch32(UTree *utree, char* filename, char* outfile, int doCollapse, int lv, int doRC) {
	FILE *fp = fopen(filename, "rb"), *fpo = fopen(outfile, "wb");
	if (fp == NULL) { puts("Invalid input files"); exit(1); }
	static const size_t LINELEN = 16777216; // 16MB lines
	uint64_t li = 0, goodFinds = 0;
	char RC[256];
	memset(RC,'N',256);
	RC['A'] = RC['a'] = 'T'; RC['C'] = RC['c'] = 'G';
	RC['G'] = RC['g'] = 'C'; RC['T'] = RC['t'] = 'A';
	char *PADTX[] = {"k__;p__;c__;o__;f__;g__;s__;t__",
					";p__;c__;o__;f__;g__;s__;t__",
					";c__;o__;f__;g__;s__;t__",
					";o__;f__;g__;s__;t__",
					";f__;g__;s__;t__",
					";g__;s__;t__",
					";s__;t__",
					";t__",
					""};
	
	// Cache important variables
	char *Dump = utree->Dump;
	uint64_t *BinIx = utree->BinIx, totalCount = utree->totalCount;
	IXTYPE maxIX = utree->sampIX + 1;
	char **SampStrings = utree->SampStrings;
	uint8_t *semicolons = utree->semicolons;
	const int k1 = PACKSIZE - 1; int kv = k1 /*+ lv*/;
	
	#define XT_INITIATE_WS() \
	char *line = malloc(LINELEN*2 + 2), *line2 = malloc(LINELEN*2 + 2), stop=0; \
	if (!line || !line2) {fputs("OOM:lineIn\n",stderr); exit(3);} \
	line[LINELEN-1] = 0, line[LINELEN*2+1] = 0; \
	line2[LINELEN-1] = 0, line2[LINELEN*2+1] = 0; \
	uint64_t ns = 0; \
	for (;;) { \
		_Pragma("omp critical") \
		{ \
			line = fgets(line,LINELEN,fp); \
			if (line) { \
				line2 = fgets(line2,LINELEN,fp); \
				if (!line2) { fprintf(stderr,"ERROR: can't read sequence L %u\n",ns); exit(2); } \
			} else stop = 1; \
		} \
		if (stop) break; \
		_Pragma("omp atomic capture") \
		ns = ++li; \
		if (!(ns & 1048575)) printf("Searched %llu queries...\n",ns); \
		char *src = line; \
		if (*src != '>') {fprintf(stderr,"ERROR: no header '>' [L %u, '%s']\n",ns,src); exit(2);} \
		while (*++src && *src != ' ' && *src != '\n'); \
		line[src-line] = 0; \
		/* fprintf(fpo,"%s\t",line+1); print sample name */ \
		/* if (ns > 729500) printf("SEQUENCE %u: %s\n",ns,line+1); */ \
		src = line2; \
		if (*src == '>') {fprintf(stderr,"ERROR: sequence begins '>' [L %u, '%s']\n",ns,src); exit(2);} \
		register int length = strlen(src); \
		if (!length) {fprintf(stderr,"ERROR: empty query line %u\n",ns); exit(2);} \
		if (src[length-1] == '\n') --length; /* lop off newline(s) */ \
		if (src[length-1] == '\r') --length; /* supports every platform! */ \
		if (doRC) { \
			src[length] = 'N'; \
			src[(length << 1)+1] = 0; \
			for (int x = length+1; x <= length << 1; ++x) \
				src[x] = RC[src[length+length-x]]; \
			length = (length << 1) + 1; \
			/* printf("Sequence with recvomp is: %s\n",src); */ \
		} \
		WTYPE w = 0; \
		size_t foundUniq = 0; \
		size_t maxSemis = 0; //, nextValid = PACKSIZE;
	#define XT_FINALIZE_WS() }
	#define XT_WORD_SEARCH() \
	/* search for words in tree */ \
	/* printf("Conidering: \n"); */ \
	for (int i = kv, z = -4; i < length; ++i) { \
		/*if (lv >= 1) { \
			if (C2Xb[src[i-kv]] != 0) continue; \
			if (lv >= 2) { \
				if (C2Xb[src[i-kv+1]] != 2) continue; \
				if (lv >= 3) { \
					if (C2Xb[src[i-kv+2]] != 1) continue; \
					if (lv >= 4) { \
						if (C2Xb[src[i-kv+3]] != 3) continue; \
					} \
				} \
			} \
		} \*/ \
		int j; \
		if (i < z + kv) w <<= (i-z-1) << 1, j = z + 1; \
		else w = 0, j = i - k1; \
		for (int p = j; j <= i; ++j) { \
			if (C2Xb[src[j]] == 255) {i += j - p /*+ lv*/, z = 0; break;} \
			w <<= 2u, w |= C2Xb[src[j]]; \
		} \
		if (j <= i) continue; \
		z = i; \
		IXTYPE ix = XT_getIX32(utree, w); \
		if (ix < maxIX) { \
			++foundUniq; \
			XT_PREP_VOTE() \
		} \
	}
	// These are prototypes of the requisite XT_PREP_VOTE
	#define XT_FULLVOTE() AllTheKingsHorses[foundUniq-1] = ix;
	#define XT_EARLYTERMINATE() \
		/* For early termination, no voting */ \
		++goodFinds; \
		fprintf(fpo,"%s\t%s\n",line+1,SampStrings[ix]); \
		break;
	#define XT_HEURISTICVOTE() \
		/* heuristic vote: break ties at finest taxon */ \
		if (semicolons[ix] > maxSemis) \
			maxSemis = semicolons[ix], \
			kingsMen = 0, *AllTheKingsHorses = ix; \
		else if (semicolons[ix] == maxSemis) \
			AllTheKingsHorses[++kingsMen] = ix;
	#define XT_SHALLOWVOTE() \
		/* simple vote: max taxid wins */ \
		i += PACKSIZE/SPARSITY - 1;  \
		AllTheKingsHorses[kingsMen++] = ix;
	#ifndef TOLERANCE_THRESHOLD
		#define TOLERANCE_THRESHOLD 2
	#endif
	#ifndef SLACK
		#define SLACK 2
	#endif
	#ifndef SPARSITY
		#define SPARSITY 4
	#endif
	//XT_INITIATE_WS()
	if (doCollapse == -1) {
		XT_INITIATE_WS()
		#define XT_PREP_VOTE() XT_EARLYTERMINATE()
		XT_WORD_SEARCH()
		#undef XT_PREP_VOTE
		XT_FINALIZE_WS()
	}
	else if (!doCollapse) {
		IXTYPE *AllTheKingsHorses = malloc(LINELEN*2*sizeof(IXTYPE));
		uint32_t *Hashes = calloc(maxIX,sizeof(*Hashes));
		int kingsMen = 0, humpty, togetherAgain;
		if (!AllTheKingsHorses || !Hashes) {
			fputs("ERROR: Init out of memory.\n",stderr); exit(3);}
		XT_INITIATE_WS()
		#define XT_PREP_VOTE() XT_SHALLOWVOTE()
		kingsMen = 0;
		XT_WORD_SEARCH()
		#undef XT_PREP_VOTE
		if (foundUniq) { // voting system, ultrafast
			++goodFinds;
			if (!kingsMen++) fprintf(fpo,"%s\t%s\n",line+1,SampStrings[*AllTheKingsHorses]);
			else {
				for (int i = 0; i < kingsMen; ++i)
					++Hashes[AllTheKingsHorses[i]];
				int most = 0, secondMost = 0;
				IXTYPE mostIX, secondMostIX = BAD_IX;
				for (int i = 0; i < kingsMen; ++i) {
					if (Hashes[AllTheKingsHorses[i]] > most)
						secondMost = most, secondMostIX = mostIX,
						mostIX = AllTheKingsHorses[i],
						most = Hashes[AllTheKingsHorses[i]];
					else if (Hashes[AllTheKingsHorses[i]] > secondMost)
						secondMostIX = AllTheKingsHorses[i],
						secondMost = Hashes[secondMostIX];
					Hashes[AllTheKingsHorses[i]] = 0;
				}
				/* if ( most < TOLERANCE_THRESHOLD || (secondMost && secondMost > most - SLACK 
					&& secondMostIX != mostIX) ) --goodFinds; */
				if ( most < TOLERANCE_THRESHOLD || most < SLACK*secondMost) --goodFinds;
				else 
					fprintf(fpo,"%s\t%s\t%f\t%d\n",line+1,SampStrings[mostIX],(double)1-(double)secondMost/most,most);
			}
		}
		//fprintf(fpo,"\n");
		XT_FINALIZE_WS()
	}
	else {
		#pragma omp parallel reduction(+:goodFinds)
		{
			// Full aufbau
			uint32_t *Hashes = calloc(maxIX,sizeof(*Hashes));
			String_Count_t *Tax_Cnt = malloc(maxIX*sizeof(*Tax_Cnt));
			char Taxon[UINT16_MAX+1] = {0};
			IXTYPE *AllTheKingsHorses = malloc(LINELEN*2*sizeof(IXTYPE));
			int kingsMen = 0, humpty, togetherAgain;
			if (!(Tax_Cnt && Hashes && Taxon && AllTheKingsHorses)) 
				{fputs("OOM:Tax_Cnt\n",stderr); exit(3);}
			
			XT_INITIATE_WS()
			#define XT_PREP_VOTE() XT_FULLVOTE()
			XT_WORD_SEARCH()
			
			#undef XT_PREP_VOTE
			#ifndef TAXACUT
			#define TAXACUT 4
			#endif
			if (foundUniq) {
				++goodFinds;
				//printf("[%u] K's found = %u\n",ns,foundUniq);
				if (foundUniq == 1) {
					fprintf(fpo,"%s\t%s\t1\t1\t*\n",line+1,SampStrings[*AllTheKingsHorses]); continue; }
				for (uint32_t i = 0; i < foundUniq; ++i) 
					++Hashes[AllTheKingsHorses[i]];
				uint32_t uix = 0;
				for (uint32_t i = foundUniq, t; i; --i) if (Hashes[t=AllTheKingsHorses[i-1]]) 
					Tax_Cnt[uix++] = (String_Count_t){SampStrings[t],Hashes[t]},
					Hashes[t] = 0;
				if (uix == 1) {
					fprintf(fpo,"%s\t%s\t%u\t1\t*\n",line+1,SampStrings[*AllTheKingsHorses],foundUniq); continue; }
				qsort(Tax_Cnt, uix, sizeof(*Tax_Cnt), byStr);
				
				//for (uint32_t i = 0; i < uix; ++i) printf("[%u] %s\n",Tax_Cnt[i].n,Tax_Cnt[i].s);
				uint32_t cutoff = foundUniq - foundUniq/TAXACUT, lv = 0, 
					st = 0, ed = uix, dv = -1, orun = foundUniq, sl, ol; // sl[INT16_MAX+1], ol[8]={0};
				cutoff += foundUniq >> 1 >= cutoff;
				for (;;) {
					uint32_t run = Tax_Cnt[st].n, td = dv, adj = 0; 
					//printf("-->lv %u, st = %u, ed = %u, dv = %u, cutoff = %u, run = %u...",lv,st,ed,dv,cutoff,run);
					for (uint32_t z = st + 1; z < ed; ++z) {
						char *s1 = Tax_Cnt[z-1].s, *s2 = Tax_Cnt[z].s;
						if (!s1[dv+(dv==-1)]) { // non-aufbau -- just reset run ('='), st=z, continue;
							run = Tax_Cnt[z].n, st = z;
							orun -= Tax_Cnt[z-1].n;
							cutoff = orun - orun/TAXACUT;
							cutoff += orun >> 1 >= cutoff;
							//--lv; // offset increment -- stay on this level with narrowed range
							continue;
						} 
						for (td = dv + 1; s1[td] && s1[td] == s2[td]; ++td) 
							if (s1[td]==';') break;
						if (s1[td] == s2[td]) run += Tax_Cnt[z].n;
						else if ((!s1[td] && s2[td]==';') || ((s1[td]==';' || !s1[td]) && s1[td-1]=='_')) // aufbau
							run = Tax_Cnt[z].n, st = z, 
							orun -= Tax_Cnt[z-1].n, 
							cutoff = orun - orun/TAXACUT,
							cutoff += orun >>1 >= cutoff;
						else if (run >= cutoff) {ed = z; break;}
						else run = Tax_Cnt[z].n, st = z;
					}
					sl = run, ol=orun; //*100/orun;
					if (run < cutoff) break;
					if (st + 1 >= ed) {
						#ifdef DEBUG
						if (!ed || ed > uix) printf("ERROR 038: [%u] ed is %u [0 %u]\n",ns,ed,uix); 
						else 
						#endif
							if (Tax_Cnt[ed-1].n >= cutoff) dv = -2, lv = INT16_MAX;
						break;
					}
					if (!Tax_Cnt[ed-1].s[td] || Tax_Cnt[ed-1].s[td]==';') ++lv, sl = run, ol = orun; //*100/orun;
					orun = run;
					dv = td;
					cutoff = run - run/TAXACUT;
					cutoff += run >> 1 >= cutoff;
				}
				char *toPrint = dv == -1 ? "" : dv == -2 ? Tax_Cnt[ed-1].s : Taxon;
				if (dv < -2) memcpy(Taxon,Tax_Cnt[ed-1].s,dv), Taxon[dv] = 0; 
				#ifdef DEBUG
				if (lv > INT16_MAX && dv < -2) {
					printf("Warning: LV %u interpolated on seq %u [%s----%s] %u, %u\n",lv,ns,toPrint,Tax_Cnt[ed-1].s,st,ed); 
					//for (uint32_t i = 0; i < uix; ++i) printf("<%u>--> [%u] %s\n",ns,Tax_Cnt[i].n,Tax_Cnt[i].s);
				}
				#endif
				//printf("--> %s---%s\n\n",toPrint,PADTX[lv]);
				fprintf(fpo,"%s\t%s\t%u\t%u\t%u;%u\n",line+1,toPrint,foundUniq,uix,sl,ol); 
					//,sl[1],ol[1],sl[2],ol[2],sl[3],ol[3],sl[4],ol[4],sl[5],ol[5],sl[6],ol[6],sl[7],ol[7]);
			}
			//fprintf(fpo,"\n");
			XT_FINALIZE_WS()
		}
	}
	//XT_FINALIZE_WS()
	fclose(fp); fclose(fpo); //free(line); // purge buffers
	//free(Hashes); //free(KingsMen);
	printf("Good finds: %llu\n",goodFinds);
	return li; // number of sequences parsed
}

UTree * UT_createTree(int numThreads) {
	initConverter();
	// Initialize memanage (comment out the below block to disable memanage)
	UBANK = malloc(sizeof(*UBANK));
	UBANK_BIN = calloc(1,sizeof(*UBANK_BIN));
	UBANK_BINCNT = malloc(sizeof(*UBANK_BINCNT));
	UBANK_IX = calloc(1,sizeof(*UBANK_IX));
	UBANK[0] = malloc(UBANK_INITBINS*sizeof(*UBANK[0]));
	UBANK_BINCNT[0] = UBANK_INITBINS;
	for (size_t j=0; j < UBANK_INITBINS; ++j) {
		UBANK[0][j] = malloc(UBANK_MAXK*sizeof(*UBANK[0][j])); // init this bin's kmers
		if (!UBANK[0][j]) {puts("error: xallocY 0"); exit(3); }
	} 
	
	// Initialize tree specific variables
	UTree *myTree = malloc(sizeof(*myTree));
	KMerY **Roots = malloc(KHASH_SIZE*sizeof(*Roots)); //prefix table
	size_t *NumsInserted = calloc(KHASH_SIZE,sizeof(*NumsInserted)),
		*TotalCounts = calloc(KHASH_SIZE,sizeof(*TotalCounts)),
		*BalanceThreshes = calloc(KHASH_SIZE,sizeof(*BalanceThreshes));
	
	// Populate initialized tree heads
	for (size_t i = 0; i < KHASH_SIZE; ++i) {
		Roots[i] = malloc(sizeof(*Roots[i]));
		*Roots[i] = (KMerY){0,EMPTY_IX,0,0};
		BalanceThreshes[i] = UBAL_THRES;
	}
	
	// Add variables to myTree
	myTree->BalanceThreshes = BalanceThreshes;
	myTree->NumsInserted = NumsInserted;
	myTree->TotalCounts = TotalCounts;
	myTree->Roots = Roots;
	myTree->numThreads = numThreads;
	myTree->numInserted = 0, myTree->totalCount = 0;
	myTree->Samps = 0, myTree->sampIX = 0;
	myTree->sampStringSz = 10, myTree->SampStrings = 0;
	myTree->SampCnts = 0;
	myTree->semicolons = 0;
	myTree->queuedClumps = 0;
	myTree->Pairs = malloc(FIRETHRES*sizeof(*myTree->Pairs));
	return myTree;
}

#define READ_ADD_SAMPLES() { \
	size_t ns=0, LINELEN = 268435456; \
	char *line = malloc(LINELEN + 1); line[LINELEN] = 0; \
	while (++ns, line = fgets(line,LINELEN,fp)) { \
		/* printf("%s",line); continue; */ \
		char *src = line; \
		/* *src != '_' && *src != ' ' && */ \
		while (*src != '\t') ++src; \
		memset(src,'\0',1); \
		/* printf("sampLine = %s [%llu]\n",line,(uint64_t)atol(src+1)); continue; */ \
		if (!ktree->Samps) { \
			ktree->Samps = malloc(sizeof(*ktree->Samps)); \
			/* int tab = -1; while (++tab, *(line+tab) != '\t'); \
			line[tab] = 0; */ \
			/*int len = strlen(line)+1; */ \
			char *str = malloc(src-line+1); \
			strcpy(str,line); \
			*ktree->Samps = (SampX){str,0,0,0}; \
			ktree->SampStrings = malloc(sizeof(*ktree->SampStrings)*ktree->sampStringSz); \
			ktree->SampCnts = calloc(ktree->sampStringSz,sizeof(*ktree->SampCnts)); \
			*ktree->SampStrings = str; \
			EXTRA_ADD() \
		} \
		else { \
			\
			ADD_FUNC(ktree,line); \
			EXTRA_ADD() \
			/* ktree->SampCnts[ktree->sampIX] = (CNTTYPE)atol(src+1); */ \
		} \
	} \
	free(line); \
}

#define EXTRA_ADD_NORMAL() {}
#define EXTRA_ADD_DELIM() { \
	if (ktree->sampStringSz > oldSize) \
		ktree->semicolons = realloc(ktree->semicolons, \
			sizeof(*ktree->semicolons)*ktree->sampStringSz), \
		oldSize = ktree->sampStringSz; \
	uint8_t semis = 0; \
	char *taxP = line - 1; \
	while (*++taxP) if (*taxP == delim) \
		if(*(taxP-1) != '_') ++semis; else break; \
	if (semis == 6 && *(taxP-1) != '_') ++semis; \
	ktree->semicolons[ktree->sampIX] = semis; \
}


int readSamplesFPdelim(UTree *ktree, FILE *fp, char delim) {
	if (!fp) { puts("Invalid input file"); exit(1); } 
	if (!delim) {
		#define EXTRA_ADD() EXTRA_ADD_NORMAL() \
			ktree->SampCnts[ktree->sampIX] = (uint64_t)atol(src+1);
		#define ADD_FUNC addSampleUdX
			READ_ADD_SAMPLES()
		#undef ADD_FUNC
		#undef EXTRA_ADD
		return 1;
	}
	if (!ktree->semicolons) ktree->semicolons = 
		malloc(sizeof(*ktree->semicolons)*ktree->sampStringSz);
	size_t oldSize = ktree->sampStringSz;
	#define EXTRA_ADD() EXTRA_ADD_DELIM() \
		ktree->SampCnts[ktree->sampIX] = (uint64_t)atol(src+1);
	#define ADD_FUNC addSampleUdX
		READ_ADD_SAMPLES()
	#undef EXTRA_ADD
	#undef ADD_FUNC
	return 1;
}

int UT_writeSamples(UTree *utree, char *filename) {
	FILE *of = fopen(filename,"wb");
	if (!of) { puts("Invalid output file"); return 0; }
	for (size_t i = 0; i < utree->sampIX+1; ++i) 
		fprintf(of,"%s\t%llu\n",utree->SampStrings[i],(uint64_t)utree->SampCnts[i]);
	fclose(of);
	return 1;
}

void XT_cmp32(char *filename, char *outfile) {
	FILE *dp = fopen(filename,"rb");
	if (!dp) { puts("Invalid input filename"); exit(0); }
	uint64_t metadata[4] = {0};
	size_t numRead = fread(metadata, sizeof(*metadata), 4, dp);
	if (numRead < 4 || !metadata[3]) {puts("Tree malformatted."); exit(0); }
	// Check file for compatibility with this compilation config
	#ifdef NO_COUNT
		#define CNT_SIZE 0
	#else
		#define CNT_SIZE sizeof(CNTTYPE)
	#endif
	size_t numNodes = metadata[3]; 
	if (metadata[0] != sizeof(WTYPE) || metadata[1] != CNT_SIZE ||
		metadata[2] != sizeof(IXTYPE)) {
			printf("ERROR. Input tree requires PACKSIZE=%u, CNTTYPE=%s, IXTYPE=%s\n", 
				metadata[0] << 2,TYPEARR[metadata[1]], TYPEARR[metadata[2]]);
		exit(0);
	}
	unsigned elemSize = sizeof(WTYPE) + CNT_SIZE + sizeof(IXTYPE);
	printf("Nodes in input tree: %llu (PACKSIZE=%u, CNTTYPE=%s, IXTYPE=%s, el=%u)\n",
		numNodes, metadata[0] << 2, TYPEARR[metadata[1]], TYPEARR[metadata[2]],elemSize);
	
	static const int DR_SZ = sizeof(WTYPE) + sizeof(IXTYPE);	
	#define DR_WORD_AT(raw,i) *(WTYPE*)((raw)+DR_SZ*(i))
	#define DR_IX_AT(raw,i) *(IXTYPE*)((raw)+DR_SZ*(i)+sizeof(WTYPE))
	#define DR_ADDR_AT(raw,i) ((raw) + DR_SZ*(i))
	//SZ
	#define DR_FIRST_WORD(raw) *(WTYPE*)(raw)
	#define DR_FIRST_IX(raw) *(IXTYPE*)((raw)+sizeof(WTYPE))
	#define DR_IXDIST(raw,new) ((new)-(raw))/DR_SZ
	
	char *Dump = malloc(numNodes * DR_SZ);
	fread(Dump,DR_SZ,numNodes,dp);
	UTree *utree = calloc(1,sizeof *utree);
	utree->sampStringSz = 10;
	int annotated = readSamplesFPdelim(utree, dp, 0);
	if (!annotated) puts("No annotation found in tree file."); 
	fclose(dp);
	
	// build N-character (2N-bit?) bin designators
	#define DR_BITS 24
	static const uint64_t bitsLeft = PACKSIZE*2 - DR_BITS, DR_NUMBINS = (1 << DR_BITS) + 1;
	#define DR_PREFIX(x) ((x) >> bitsLeft)
	if (numNodes < UINT32_MAX) puts("Using 32-bit counters");
	else puts("Holy smokes. Looks like we have over 4 billion k-mers here."),
		puts("Please complain to the developer.\nTrying something anyway...");
	uint64_t *BinIx = calloc(DR_NUMBINS,sizeof(*BinIx));
	for (size_t i = 0; i < numNodes; ++i) {
		uint32_t v = DR_PREFIX(DR_WORD_AT(Dump,i));
		if (!BinIx[v]) BinIx[v] = i;
		//printf("prefix=%u\n",(uint32_t)v);
	}
	BinIx[DR_NUMBINS-1] = numNodes;
	size_t u = 0; for (; !BinIx[u]; ++u); BinIx[u] = 0;
	for (size_t i = DR_NUMBINS-2; i > u; --i) if (!BinIx[i]) BinIx[i] = BinIx[i+1];
	
	// Compress and store the tree
	int cmpStr = 3; // curSz; //2
	char *cmask = malloc(sizeof(WTYPE)); // round up to whole word size of WTYPE
	int curSz = sizeof(WTYPE), cmpSz = sizeof(WTYPE) - cmpStr; //= sizeof(Darkness);
	for (int i = 0; i < cmpSz; ++i) cmask[i] = 0xff;
	for (int i = cmpSz; i < curSz; ++i) cmask[i] = 0x00;
	WTYPE mask = *(WTYPE*)cmask;
	FILE *prefixPtrs = fopen(outfile,"wb");
	if (!prefixPtrs) { puts("Invalid output filename"); exit(0); }
	
	fwrite(metadata,sizeof *metadata,4,prefixPtrs);
	//fwrite(BinIx,sizeof *BinIx, DR_NUMBINS, prefixPtrs);
	int ixSize = numNodes < UINT32_MAX ? sizeof(uint32_t) : sizeof(uint64_t);
	for (size_t i = 0; i < DR_NUMBINS; ++i) fwrite(&BinIx[i],ixSize,1,prefixPtrs);
	//for (size_t i = 0; i < DR_NUMBINS; ++i) printf("%llu\t%llu\n",i,BinIx[i]);
	for (size_t i = 0; i < numNodes; ++i) {
		fwrite(DR_ADDR_AT(Dump,i),1,cmpSz,prefixPtrs);
		fwrite(DR_ADDR_AT(Dump,i)+sizeof(WTYPE),sizeof(IXTYPE),1,prefixPtrs);
	}
	uint64_t total = 0;
	for (size_t i = 0; i < utree->sampIX+1; ++i) 
		total += utree->SampCnts[i], 
		fprintf(prefixPtrs,"%s\t%llu\n",utree->SampStrings[i],(uint64_t)utree->SampCnts[i]);
	printf("Total nodes in tree: %llu [%llu labels]\n",total,(uint64_t)utree->sampIX+1);
}

int UT_writeTreeBinary(UTree *utree, char* filename) {
	FILE *of = fopen(filename,"wb");
	if (!of) { puts("Invalid output filename"); return 0; }
	uint64_t metadata[4] = {sizeof(WTYPE), 0, sizeof(IXTYPE),0};
	fwrite(metadata,sizeof(*metadata),4,of);
	if (!utree->SampCnts) {
		utree->SampCnts = calloc(utree->sampIX + 1,sizeof(*utree->SampCnts));
		for (size_t i = 0; i < KHASH_SIZE; ++i) {
			traceTreeBUf(i,utree->Roots[i], utree->SampCnts, of);
		}
	}
	else { 
		for (size_t i = 0; i < KHASH_SIZE; ++i) {
			traceTreeBU(i,utree->Roots[i], of);
		}
	}
	// Embed sample info into tree
	uint64_t total = 0;
	for (size_t i = 0; i < utree->sampIX+1; ++i) 
		total += utree->SampCnts[i], 
		fprintf(of,"%s\t%llu\n",utree->SampStrings[i],(uint64_t)utree->SampCnts[i]);
	printf("Total nodes in tree: %llu [%llu labels]\n",total,(uint64_t)utree->sampIX+1);
	fseek(of,3*sizeof(uint64_t),SEEK_SET);
	fwrite(&total,sizeof(total),1,of);
	fclose(of);
	return 1;
}
#if defined BUILD_GG || defined SEARCH_GG
	#define DO_GG 1
#else
	#define DO_GG 0
#endif
// usage: RadixalTriecrobium [see module usage]
#define VER "[v2.0RF SigNature Edition]"
int main(int argc, char *argv[]) { 
	#ifdef COMPRESS
	if (argc != 3) {puts(VER " usage: xtree-compress preTree.ubt compTree.ctr"); exit(1); }
	XT_cmp32(argv[1],argv[2]); //("outNEW.cbt","compTre.ctr");
	exit(0);
	#endif
	#if defined SEARCH || defined SEARCH_GG
	if (argc < 4) {
		printf(VER " usage: xtree-search%s compTree.ctr fastaToSearch.fa output.txt [threads] [SPEED <X>] [RC]\n",
			DO_GG ? "GG" : ""); exit(1); }
	printf("This is UTree " VER "\n");
	int doRC = !strcmp(argv[argc-1],"RC"), threads = 1;
	argc -= doRC; // get rid of last commandline
	int speed = 0; if (!strcmp(argv[argc-2],"SPEED")) speed = atoi(argv[argc-1]), argc -= 2;
	printf("Reverse complement consideration is %sabled.\n",doRC? "en" : "dis");
	printf("Searching at speed %d.\n",speed);
	#ifdef _OPENMP
		threads = argc >= 5 ? atoi(argv[4]) : omp_get_max_threads();
		omp_set_num_threads(threads);
		printf("Using up to %d threads.\n",threads);
	#else
		puts("Multi-threading is disabled in this build.");
	#endif
	UTree *xtr = XT_read32(argv[1], DO_GG ? ';' : 0); //("compTre.ctr");
	printf("Searched %llu queries\n",
		XT_doSearch32(xtr,argv[2],argv[3], DO_GG ? 8 : 0, speed, doRC));
	exit(0);
	#endif
	#if defined BUILD || defined BUILD_GG
	if (argc < 5) {
		printf(VER " usage: utree-build%s input_fasta.fa labels.map output.ubt threads{0=auto} [complevel]\n",
			DO_GG ? "GG" : ""); exit(1); }
	printf("This is UTree " VER "\n");
	char *filename = argv[1], *dbname = argv[2], *outname = argv[3];
	int threads = 1;
	#ifdef _OPENMP
		threads = atoi(argv[4]) ?: omp_get_max_threads();
		omp_set_num_threads(threads);
		printf("Using up to %d threads.\n",threads);
	#else
		puts("Multi-threading is disabled in this build.");
	#endif
	
	UTree *myU = UT_createTree(threads);
	puts("Tree initialized.");
	uint32_t cl = 1;
	if (argc > 5) cl = atoi(argv[5]);
	printf("Setting compression level to %u\n",cl);
	UT_parseSampFastaExternOSFA(myU,filename,dbname,0,1,cl,DO_GG); // 0,2 for ixcol, ixlab in img
	puts("File parsed.");
	UT_writeTreeBinary(myU, outname);
	puts("Tree written.");
	
	char *logFN = calloc(4097,1);
	sprintf(logFN,"%s%s.log",outname, DO_GG ? ".gg" : "");
	UT_writeSamples(myU, logFN);
	exit(0);
	#endif
	
	puts("Nothing to see here.");
	puts("Compile with one of: -D BUILD, BUILD_GG, SEARCH, SEARCH_GG, COMPRESS");
	return 1;
}
