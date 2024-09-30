#ifndef HTS_HEADER
#define HTS_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <stdbool.h>

// Information on the bases
#define N_BASES 4
#define BASES "ACGT"
// Infomration on the bases, including deletions
#define N_DELBASES 5
#define DELBASES "ACGT-"

#define BAM_BASES "=ACMGRSVTWYHKDBN";
// For the main counting function
#define QUERY_SKIPPED 1
#define QUERY_PROCESSED 0
// Indices of the bases in BAM files
#define BAM_A 1 // A
#define BAM_C 2 // C
#define BAM_G 4 // G
#define BAM_T 8 // U
#define BAM_N 15 // A,C,G,U
// Indices of the bases in the output arrays
#define IX_A 0
#define IX_C 1
#define IX_G 2
#define IX_T 3
#define IX_U 3
#define IX_DEL 4
// Default counting parameters
#define DEFAULT_MIN_MAPPING_QUALITY 10
#define DEFAULT_MIN_NUCLEOTIDE_QUALITY 20
#define DEFAULT_MIN_QUERY_LENGTH 0
#define DEFAULT_MAX_DELETION_LENGTH 10
#define DEFAULT_COLLAPSE_MUTATIONS 0
#define DEFAULT_SPREAD_DELETIONS false
#define DEFAULT_NEIGHBOURS_TO_CHECK 0



typedef struct {
    const char *filename;
    faidx_t *fai;
    int numSequences; 
} IndexedFASTA;

IndexedFASTA openIndexedFASTA(const char *filename);

int getSequenceLength(IndexedFASTA ixFASTA, int ix);

char *getSequence(IndexedFASTA ixFASTA, int ix);

bool verifyHomogeneous(IndexedFASTA ixFASTA);

void closeIndexedFASTA(IndexedFASTA ixFASTA);



typedef struct {
    const char *filename;
    htsFile *file;
    bam1_t *alignment;
    hts_idx_t *index;
    bam_hdr_t *header;
    int numReferences;
    int numAlignments;
    int numUnmappedReads;
    int maxReferenceLength;
} IndexedBAM;

IndexedBAM openIndexedBAM(const char* filename);

void closeIndexedBAM(IndexedBAM bam);



typedef struct {
    int minMapQ;
    int minBaseQ;
    int maxDelLength;
    int minQueryLength;
    int collapseMutations;
    bool spreadDeletions;
    int numNeighboursToCheck;
} CountingFlags;

CountingFlags defaultCountingFlags();

bool verifyCountingFlags(CountingFlags cFlags);

static int accumulateMutations(
    float *mutations,
    float *insertions,
    char *referenceSequence,
    IndexedBAM ixBAM,
    CountingFlags cFlags
);

typedef struct {
    uint64_t numAlignments;
    uint64_t numSkipped;
} AlignmentCount;

AlignmentCount countMutations(
    float *mutations,
    float *insertions,
    char *referenceSequence,
    IndexedBAM ixBAM,
    int ix, 
    CountingFlags cFlags
);

void embedSequence(
    const char *sequence,
    float *embeddedSequence
);

#endif
