#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    
    if (argc == 1) {
        printf("0\n");
        return EXIT_SUCCESS;
    } else {
        for (int i = 1; i < argc; i++) {
            size_t length = strlen(argv[i]);
            printf("%zu", length);
            if (i != argc -1) {
                printf(",");
            }
        }
    }
    printf("\n");

    return EXIT_SUCCESS;

}
