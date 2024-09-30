#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
// OMP
#include <omp.h>
// getopt
#include <getopt.h>
#include <errno.h>
// For printing integers with commas
#include <locale.h>
// Processing indexed BAM and FASTA files
#include "htsutils.h"
// HDF5
#include "h5utils.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))


void embedFASTA(
    IndexedFASTA ixFASTA,
    Dataspace dataspace,
    Memspace memspace
) {

    for (int ix_b = 0; ix_b < ixFASTA.numSequences; ix_b+=CHUNK_SIZE) {

        int loopSize = MIN(CHUNK_SIZE, ixFASTA.numSequences - ix_b);
        resizeMemspace(&memspace, loopSize);
        fillMemspace(memspace, 0.0);
        float *embeddedSequences = memspace.data;
        
        for (int ix = ix_b; ix < ix_b + loopSize; ix++) {

            char *sequence = getSequence(ixFASTA, ix);
            embedSequence(sequence, embeddedSequences);
            
            embeddedSequences += memspace.numElementsPerUnit;

            free(sequence);

        }

        writeSlab(dataspace, memspace, ix_b);

    }

}


bool checkAccess(const char *filename) {
    if (access(filename, F_OK) == 0) {
        return true;
    }
    return false;
}


int main(int argc, char **argv) {

    if (argc != 3) {
        printf("Error.\n");
        exit(EXIT_FAILURE);
    }

    const char *inFile = argv[1];
    const char *outFile = argv[2];

    // Open the FASTA file
    IndexedFASTA ixFASTA = openIndexedFASTA(inFile);

    if (!verifyHomogeneous(ixFASTA)) {
        fprintf(stderr, "The FASTA file is not homogeneous.\n");
        exit(EXIT_FAILURE);
    }



    bool overwrite = true;
    DataspaceFlags dFlags = defaultDataspaceFlags();
    dFlags.group = "grp";

    // Check if the output file exists; if it does, check if the
    // user wants it deleted
    bool outFileExists = false;
    if (checkAccess(outFile)) {
        if (overwrite) {
            int status = remove(outFile);
            if (status != 0) {
                printf("There was an error deleting the existing output file \"%s\".\n", outFile);
                exit(EXIT_FAILURE);
            }
        } else {
            // If the group is not specified, we cannot proceed
            if (!dFlags.group) {
                printf("The file \"%s\" already exists, and a group has not been specified. Please use the --overwrite flag to overwrite the file, or specify a group to write to the existing file with the --group option.\n", outFile);
                exit(EXIT_FAILURE);
            } 
            // Else, we check if the specified group exists in the file
            else {
                // Open the filespace
                Filespace filespace = openFilespace(outFile, H5F_ACC_RDONLY);
                if (groupExists(filespace, dFlags.group)) {
                    printf("The group \"%s\" in the file \"%s\" already exists. Please use the --overwrite flag to overwrite the file, or choose another group name to write to the existing file with the --group option.\n", dFlags.group, outFile);
                    closeFilespace(filespace);
                    exit(EXIT_FAILURE);
                } 
                // Else, we need not do anything
                closeFilespace(filespace);
                outFileExists = true;
            }
        }
    }

    // Get the first sequence in the file
    char *sequence = getSequence(ixFASTA, 0);
    size_t seqLength = strlen(sequence);
    // Create an HDF5 file to store the output
    Filespace filespace = makeFilespace(outFile);

    hsize_t dims[2] = {ixFASTA.numSequences, seqLength};
    Dataspace dataspace = makeDataspace(
        filespace, 
        "sequence", 
        2, 
        dims, 
        dFlags
    );
    Memspace memspace = initMemspace(dataspace, CHUNK_SIZE);

    embedFASTA(ixFASTA, dataspace, memspace);

    return EXIT_SUCCESS;

}
