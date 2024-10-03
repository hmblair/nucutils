#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv) {

    if (argc > 4 || argc < 3) {
        printf("  usage: %s STRING POS1 [POS2]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int p1 = atoi(argv[2]);
    int length = strlen(argv[1]);

    // Verify the index
    if (p1 < -length || p1 >= length) {
        fprintf(stderr, "  The indices must be between %d and %d inclusive.\n", -length, length-1);
        return EXIT_FAILURE;
    }
    // Allow for negative indexing
    if (p1 < 0) {
        p1 = length + p1;
    }

    if (argc == 4) {
        int p2 = atoi(argv[3]);

        // Verify the index
        if (p2 < -length || p2 > length) {
            fprintf(stderr, "  The indices must be between %d and %d inclusive.\n", -length, length-1);
            return EXIT_FAILURE;
        }
        // Allow for negative indexing
        if (p2 < 0) {
            p2 = length + p2;
        }

        for (int i = p1; i < p2; i++) {
            printf("%c", argv[1][i]);
        }
        printf("\n");
    } else {
        printf("%c\n", argv[1][p1]);
    }

    return EXIT_SUCCESS;

}
