/*********

File        mafftBorderLocator.c
Author      EPW <estebanpw@uma.es> MSR <marcossr@uma.es>
Description Plots a mafft alignment with and without gapps with colors (even for gigabyte-sized files) and their own conservation plots (4 png images altogether)

Requires: 	libpng to be installed. You can obtain it by running "sudo apt-get install libpng-dev"

USAGE       ./mafft2plot <mafft alignment file>

MAKEFILE	g++ -g main.cpp pngwriter/src/pngwriter.h -o mafft2plot `freetype-config --cflags` -l:libPNGwriter.a -lm -lpng -lz -lfreetype

**********/

/* PREPROCESSOR COMMANDS */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <string.h>
#include <ctype.h>

using namespace std;

/* FUNCTIONS */
uint64_t get_seq_len(FILE * f);
uint64_t get_number_of_seqs(FILE * f);
void terror(const char * s);
char * get_basename(char * a, int file_length);

/* VARIABLES */
//#define PERCENTAGE 90
//#define MAX_HIT 5
// Output file name
char name[1024];

/* STATEMENTS & EXPRESSIONS */
int main(int argc, char ** av) {
    FILE * database = NULL;
    name[0] = '\0';

    if (argc != 4) {
        terror("USE: mafftBorderLocator <conservation percentage> <max consecutive hits> <mafft alignment file>\n");
    }

    int PERCENTAGE = atoi(av[1]);
    int MAX_HIT = atoi(av[2]);
    printf("--> mafftBorderLocator: %d%%, %d hits\n", PERCENTAGE, MAX_HIT);

    database = fopen(av[3], "rt");
    int n_seqs = (int) get_number_of_seqs(database);

    //fprintf(stdout, "Output file name is %s.png\n", get_basename(av[1], strlen(av[1])));
    //strcpy(name, get_basename(av[1], strlen(av[1])));
    //strcat(name, ".png");

    char c = 'N'; //Char to read character

    uint64_t file_length = get_seq_len(database);

    fprintf(stdout, "There are %d seqs of length %" PRIu64 "\n", n_seqs, file_length);

    char ** mat_nucl = (char ** ) malloc(n_seqs * sizeof(char * ));

    int i, j;
    for (i = 0; i < n_seqs; i++) {
        mat_nucl[i] = (char * ) malloc(file_length * sizeof(char));
	for (j = 0; j < file_length; j++) {
            mat_nucl[i][j] = '\0';
        }
    }

    // cargar los datos en una matriz
    printf("Reading file...\n");
    int x = -1, y = 0;
    c = fgetc(database);
    while ((!feof(database))) {
        if (c == '>') {
            ++x;
            // mostrar porcentaje restante
            if ((x * 100 / n_seqs) % 5 == 0) {
                printf("%d%%\r", (x * 100) / n_seqs);
            }
            y = 0;

            //fprintf(stdout, "\n");

            while (c != '\n') {
                //printf("%c", c);
                c = fgetc(database);
            }

            while (c != '>' && (!feof(database))) { //Until next id
                c = fgetc(database);
                c = toupper(c);
                //printf("read %c\n", c);
                //getchar();
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' || c == '-') {
                    mat_nucl[x][y] = c;
                    ++y;
                }
            }
        } else {
            c = fgetc(database);
        }
    }
	printf("Done\n");

	printf("Getting first cut...\n");
    // Recorrer la maldita matriz desde el principio
    float conservation_start = 0.0;
    float conservation_end = 0.0;
    int w, h, position_start, position_end;
    int hit = 0;
    int times = 0;
    for (w = 0; w < file_length; w++) {
        // contar de que tipo de nucleotido hay
        int count[4] = {0, 0, 0, 0};
        //char nucl[5] = "ACGT";
        for (h = 0; h < n_seqs; h++) {
            switch (mat_nucl[h][w]) {
                case 'A':
                    count[0]++;
                    break;
                case 'C':
                    count[1]++;
                    break;
                case 'G':
                    count[2]++;
                    break;
                case 'T':
                    count[3]++;
                    break;
                default:
                    break;
            }
        }
        // Get max
        int i = 0;
        float maximum = count[i];
        // int index = 0;
        for (i = 0; i < 4; i++) {
            if (count[i] > maximum) {
                maximum = count[i];
                // index = i;
            }
        }

        conservation_start = (maximum * 100 / n_seqs);

        if(conservation_start >= PERCENTAGE){
	    times++;
            //printf("1st#%d: %d\n", times, w);
            position_start = w;
	    hit++;
        }else{
	    hit=0;
	}

	if(hit == MAX_HIT){
	    hit = 0;
	    times = 0;
	    break;
	}
    }

    // Recorrer la maldita matriz desde el final
    printf("Getting second cut...\n");
    for (w = file_length; w > 0; w--) {
        // contar de que tipo de nucleotido hay
        int count[4] = {0, 0, 0, 0};
        //char nucl[5] = "ACGT";
        for (h = 0; h < n_seqs; h++) {
            switch (mat_nucl[h][w]) {
                case 'A':
                    count[0]++;
                    break;
                case 'C':
                    count[1]++;
                    break;
                case 'G':
                    count[2]++;
                    break;
                case 'T':
                    count[3]++;
                    break;
                default:
                    break;
            }
        }
        // Get max
        int i = 0;
        float maximum = count[i];
        // int index = 0;
        for (i = 0; i < 4; i++) {
            if (count[i] > maximum) {
                maximum = count[i];
                // index = i;
            }
        }

        conservation_end = (maximum * 100 / n_seqs);

        if(conservation_end >= PERCENTAGE){
	    times++;
            //printf("2nd#%d: %d\n", times, w);
            position_end = w;
            hit++;
        }else{
	    hit=0;
	}

        if(hit == MAX_HIT){ 
            hit = 0;
	    times = 0;
	    break;
        }

    }
    float percentage_start = (position_start*100)/file_length;
    float percentage_end = (position_end*100)/file_length;

    printf("1st %d [%.0f%%] \n2nd %d [%.0f%%(%.0f%%)]\n",
        position_start, percentage_start,
        position_end, percentage_end, 100-percentage_end);

    /*printf("%d %d %.0f%% %.0f%%(%.0f%%) %d\n", position_start, position_end, 
	percentage_start, percentage_end, 100-percentage_end, 
	position_end-position_start);*/

    fclose(database);

    /* VARIABLES ROCKING IN A FREE WORLD */
    int f;
    for(f=0; f<n_seqs; f++){
	free(mat_nucl[f]);
    }
    free(mat_nucl);

    return 0;
}

uint64_t get_number_of_seqs(FILE * f) {
    char c = '\0';
    uint64_t n = 0;
    while (!feof(f)) {
        c = fgetc(f);
        if (c == '>') n++;
    }

    rewind(f);
    return n;
}

uint64_t get_seq_len(FILE * f) {
    char c = '\0';
    uint64_t file_length = 0;

    while (c != '>') c = fgetc(f);
    c = '\0';

    while (c != '\n') c = fgetc(f);

    while (c != '>' && !feof(f)) {
        c = getc(f);
        if (c == '>') {
            break;
        }
        c = toupper(c);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N' || c == '-') {
            ++file_length;
        }
    }

    rewind(f);
    return file_length;
}

void terror(const char * s) {
    printf("%s", s);
    exit(-1);
}

char * get_basename(char * a, int file_length) {
    int p = file_length - 1;
    while (a[p] != '/' && p >= 0) p--;
    return &a[p + 1];
}
