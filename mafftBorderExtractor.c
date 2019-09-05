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
    FILE * file_beginning = NULL;
    FILE * file_ending = NULL;

    name[0] = '\0';

    if (argc != 4) {
        terror("USE: mafftBorderExtractor <first cut> <second cut> <mafft alignment file>\n");
    }

    int first_cut = atoi(av[1]);
    int second_cut = atoi(av[2]);

    printf("--> mafftBorderExtractor: %d, %d\n", first_cut, second_cut);

    database = fopen(av[3], "rt");
    file_beginning = fopen("border_beg.mafft", "wt");
    file_ending = fopen("border_end.mafft", "wt");

    int n_seqs = (int) get_number_of_seqs(database);

    //fprintf(stdout, "Output file name is %s.png\n", get_basename(av[1], strlen(av[1])));
    //strcpy(name, get_basename(av[1], strlen(av[1])));
    //strcat(name, ".png");

    char c = 'N'; //Char to read character

    uint64_t file_length = get_seq_len(database);

    fprintf(stdout, "There are %d seqs of length %" PRIu64 "\n", n_seqs, file_length);

    // cargar los datos en una matriz
    printf("Loading left border...\n");
    int x = -1;
    int count = 0;
    c = fgetc(database);
    while ((!feof(database))) {
        if (c == '>') {
            ++x;
            // mostrar porcentaje restante
            if ((x * 100 / n_seqs) % 5 == 0) {
                printf("%d%%\r", (x * 100) / n_seqs);
            }
            
            // copiar el nombre de la secuencia
            while (c != '\n') {
		fprintf(file_beginning, "%c", c);
		c = fgetc(database);
		//fprintf(file_beginning, "%c", c);
            }

	    // copiar el salto de linea
	    if(c == '\n'){
		fprintf(file_beginning, "%c", c);
	    }
            
	    // copiar la secuencia hasta el corte
	    while ((!feof(database)) && (c != '>') && (count < first_cut)) { //Until next id
                c = fgetc(database);
		if((c != '\n') && (c != '-')){
		    fprintf(file_beginning, "%c", c);
		}
		count++;
            }
	  	
            // si ya se ha recortado la primera secuencia se reinicia
	    if(count == first_cut){
		count = 0;
		fprintf(file_beginning, "%c", '\n');
	    }
        } else {
	    // en caso de que algun subnormal haya metido un espacio/salto de linea indeseado
            c = fgetc(database);

        }
    }
    printf("Done\n");
    rewind(database);

    // cargar los datos en una matriz
    printf("Loading right border...\n");
    x = -1;
    count = 0;
    c = fgetc(database);
    while ((!feof(database))) {
        if (c == '>') {
            ++x;
            // mostrar porcentaje restante
            if ((x * 100 / n_seqs) % 5 == 0) {
                printf("%d%%\r", (x * 100) / n_seqs);
            }

            // copiar el nombre de la secuencia
            while (c != '\n') {
                fprintf(file_ending, "%c", c);
                c = fgetc(database);
                //fprintf(file_beginning, "%c", c);
            }

            // copiar el salto de linea
            if(c == '\n'){
                fprintf(file_ending, "%c", c);
            }
            // saltar la secuencia hasta el corte
	    while ((!feof(database)) && (c != '>') && (count < second_cut)){
		c = fgetc(database);
		count++;
	    }

	    // copiar a partir del corte
            while ((!feof(database)) && (c != '>') && (count >= second_cut)) { //Until next id
                c = fgetc(database);
                if((c != '\n') && (c != '-') && (c != '>')){
                    fprintf(file_ending, "%c", c);
                }
		if(c == '>'){
		    fprintf(file_ending, "%c", '\n');
		}
                count++;
            }

            // si ya se ha recortado la primera secuencia se reinicia
            if(count == file_length){
                count = 0;
                fprintf(file_ending, "%c", '\n');
            }
        } else {
            // en caso de que algun subnormal haya metido un espacio/salto de linea indeseado
            c = fgetc(database);

        }
    }
    printf("Done\n");

    fclose(file_beginning);
    fclose(file_ending);
    fclose(database);

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
