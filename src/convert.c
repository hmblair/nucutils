
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE 4096
#define HEADER_START '>'

int main(int argc, char *argv[]) {

    if (argc != 4) {
        fprintf(stderr, "  usage: %s <filename> <origBase> <newBase>\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (strlen(argv[2]) != 1 || strlen(argv[3]) != 1) {
        fprintf(stderr, "  Both bases must be single characters.\n");
        return EXIT_FAILURE;
    }

    char origBase = argv[2][0];
    char newBase = argv[3][0];

    FILE *file = fopen(argv[1], "rb+");
    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    char buffer[BUFFER_SIZE];
    size_t bytesRead;
    long pos;

    while ((bytesRead = fread(buffer, 1, BUFFER_SIZE, file)) > 0) {

        // Store the current file position
        pos = ftell(file);

        // Process each byte in the buffer
        for (size_t i = 0; i < bytesRead; i++) {
            // Check if the current position is the start of a line
            if (i == 0 || buffer[i - 1] == '\n') {
                // If the line starts with '>', skip to the next line
                if (buffer[i] == HEADER_START) {
                    // Find the end of the current line and then continue
                    while (i < bytesRead && buffer[i] != '\n') {
                        i++;
                    }
                    continue;
                }
            }

            // Replace occurrences of origBase with newBase
            if (buffer[i] == origBase) {
                buffer[i] = newBase;
            }
        }

        // Move the file pointer back to the start of the current buffer
        fseek(file, pos - bytesRead, SEEK_SET);
        // Write the modified buffer back to the file
        fwrite(buffer, 1, bytesRead, file);
        fflush(file);
    }

    fclose(file);
    return EXIT_SUCCESS;

}

