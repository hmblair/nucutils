#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define BUFFER_SIZE 4096

int process_file(const char *filename) {
    FILE *input = fopen(filename, "r");
    if (input == NULL) {
        perror("  Error opening input file");
        return EXIT_FAILURE;
    }

    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", filename);
    FILE *output = fopen(temp_filename, "w");
    if (output == NULL) {
        perror("  Error opening temporary file");
        fclose(input);
        return EXIT_FAILURE;
    }

    char buffer[BUFFER_SIZE];
    char prev_buffer[BUFFER_SIZE] = "";
    int is_first_line = 1;

    while (fgets(buffer, sizeof(buffer), input)) {
        if (!is_first_line && buffer[0] != '>' && prev_buffer[0] != '>') {
            // Remove newline from previous line
            size_t len = strlen(prev_buffer);
            if (len > 0 && prev_buffer[len - 1] == '\n') {
                prev_buffer[len - 1] = '\0';
            }
            fputs(prev_buffer, output);
        } else if (!is_first_line) {
            fputs(prev_buffer, output);
        }

        strcpy(prev_buffer, buffer);
        is_first_line = 0;
    }

    // Write the last line
    fputs(prev_buffer, output);

    fclose(input);
    fclose(output);

    // Replace original file with processed file
    if (rename(temp_filename, filename) != 0) {
        perror("  Error renaming temporary file");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
        fprintf(stderr, "  usage: %s FASTA\n", argv[0]);
        return EXIT_FAILURE;
    }

    return process_file(argv[1]);
}
