/*********

File        mafft2plotComplete.c
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
#include "pngwriter/src/pngwriter.h"

using namespace std;

/* FUNCTIONS */
uint64_t get_seq_len(FILE * f);
uint64_t get_number_of_seqs(FILE * f);
int remove_gaps(float ** mat, float ** mat_ungapped, int wide, int height, int percentage, float * conservation, float * conservation_ungapped);
void convolve(float * altura, float * alturasmooth, int width, int window_size);
void plot(float ** mat, int wide, int height, char * name);

void terror(const char * s);
inline void get_color(char a, float * c1, float * c2, float * c3);
char * get_basename(char * a, int file_length);
void free_mat(void *** mat, int width);

/* VARIABLES */
// Minimum vertical size of plotting area
#define MIN_HEIGHT 3000
#define MAX_PERCENTAGE 100
#define TRACE_WINDOW_SIZE 135

// Output file name
char name[1024];

/* STATEMENTS & EXPRESSIONS */
int main(int argc, char ** av) {
    FILE * database = NULL;
    name[0] = '\0';

    if (argc != 4) {
        terror("USE: mafft2plot <mafft alignment file> <ungapped percentage> <window size>\n");
    }

    database = fopen(av[1], "rt");
    int n_seqs = (int) get_number_of_seqs(database);

    /*fprintf(stdout, "Output file name is %s.png\n", get_basename(av[1], strlen(av[1])));
    strcpy(name, get_basename(av[1], strlen(av[1])));
    strcat(name, ".png");
    */
    int ungapped_per = atoi(av[2]);
    int window_size = atoi(av[3]);

    char c = 'N'; //Char to read character

    uint64_t file_length = get_seq_len(database);

    fprintf(stdout, "There are %d seqs of length %"
    PRIu64 "\n", n_seqs, file_length);

    float ** mat = (float ** ) malloc(n_seqs * sizeof(float * ));
    float ** mat_ungapped = (float ** ) malloc(n_seqs * sizeof(float * ));

    char ** mat_nucl = (char ** ) malloc(n_seqs * sizeof(char * ));
    float ** mat_nucl_canvas = (float ** ) malloc((MAX_PERCENTAGE + 1) * sizeof(float * ));

    int k, j;
    for (k = 0; k < n_seqs; k++) {
        mat[k] = (float * ) malloc(3 * file_length * sizeof(float));
        mat_nucl[k] = (char * ) malloc(file_length * sizeof(char));
        mat_ungapped[k] = (float * ) malloc(3 * file_length * sizeof(float));
        for (j = 0; j < 3 * file_length; j++) {
            mat[k][j] = 0.0;
            mat_ungapped[k][j] = 0.0;
        }
        for (j = 0; j < file_length; j++) {
            mat_nucl[k][j] = '\0';
        }
    }

    // Rellenamos la matriz de color blanco
    for (k = 0; k < MAX_PERCENTAGE + 1; k++) {
        mat_nucl_canvas[k] = (float * ) malloc(3 * file_length * sizeof(float));
        for (j = 0; j < 3 * file_length; j++) {
            mat_nucl_canvas[k][j] = 1.0;
        }
    }

    int x = -1, y = 0;
    c = fgetc(database);
    while ((!feof(database))) {
        if (c == '>') {
            ++x;
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
                    float c1, c2, c3;
                    get_color(c, & c1, & c2, & c3);
                    mat[x][3 * y] = c1;
                    mat[x][3 * y + 1] = c2;
                    mat[x][3 * y + 2] = c3;

                    // para la conservacion
                    mat_nucl[x][y] = c;
                    //fprintf(stdout, ",%d,%d,%d", c1, c2, c3);
                    ++y;
                }
            }
        } else {
            c = fgetc(database);
        }
    }

    // Vector de conservacion
    float * conservation = (float * ) malloc(file_length * sizeof(float));
    float * conservation_ungapped = (float * ) malloc(file_length * sizeof(float));
    float * conservation_smooth = (float * ) malloc(file_length * sizeof(float));
    float * conservation_trace = (float * ) malloc(file_length * sizeof(float));

    uint64_t s, w;
    for (w = 0; w < file_length; w++) {
        int count[4] = {0, 0, 0, 0};
        char nucl[5] = "ACGT";
        for (s = 0; s < n_seqs; s++) {
            switch (mat_nucl[s][w]) {
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

        // Escalado
        conservation[w] = (maximum * 100 / n_seqs);
    }

    convolve(conservation, conservation_smooth, file_length, window_size);
    convolve(conservation_smooth, conservation_trace, file_length, TRACE_WINDOW_SIZE);

    // Rellenar canvas
    int c_s = 0, c_t = 0;
    for (w = 0; w < file_length; w++) {
        c_s = (int) conservation_smooth[w];
        mat_nucl_canvas[c_s][3 * w] = 1.0;
        mat_nucl_canvas[c_s][3 * w + 1] = 0.0;
        mat_nucl_canvas[c_s][3 * w + 2] = 0.0;

        c_t = (int) conservation_trace[w];
        mat_nucl_canvas[c_t][3 * w] = 0.741;
        mat_nucl_canvas[c_t][3 * w + 1] = 0.741;
        mat_nucl_canvas[c_t][3 * w + 2] = 0.741;
    }

    sprintf(name, "%s_%s_%s_alignment.png", get_basename(av[1], strlen(av[1])), av[2], av[3]);
    fprintf(stdout, "Output file name is %s\n", name);
    plot(mat, file_length, n_seqs, name);
	
    sprintf(name, "%s_%s_%s_conservation.png", get_basename(av[1], strlen(av[1])), av[2], av[3]);
    fprintf(stdout, "Output file name is %s\n", name);
    plot(mat_nucl_canvas, file_length, MAX_PERCENTAGE, name);

    for (k = 0; k < MAX_PERCENTAGE + 1; k++) {
        for (j = 0; j < 3 * file_length; j++) {
            mat_nucl_canvas[k][j] = 1.0;
        }
    }

    /* SMOOTHING BABY */
    int new_length = remove_gaps(mat, mat_ungapped, file_length, n_seqs, ungapped_per, conservation, conservation_ungapped);

    convolve(conservation_ungapped, conservation_smooth, new_length, window_size);
    convolve(conservation_smooth, conservation_trace, new_length, TRACE_WINDOW_SIZE);

    // Rellenar canvas
    c_s = 0, c_t = 0;
    for (w = 0; w < new_length; w++) {
        c_s = (int) conservation_smooth[w];
        mat_nucl_canvas[c_s][3 * w] = 1.0;
        mat_nucl_canvas[c_s][3 * w + 1] = 0.0;
        mat_nucl_canvas[c_s][3 * w + 2] = 0.0;

        c_t = (int) conservation_trace[w];
        mat_nucl_canvas[c_t][3 * w] = 0.741;
        mat_nucl_canvas[c_t][3 * w + 1] = 0.741;
        mat_nucl_canvas[c_t][3 * w + 2] = 0.741;
    }

    sprintf(name, "%s_%s_%s_alignment_ungapped.png", get_basename(av[1], strlen(av[1])), av[2], av[3]);
    fprintf(stdout, "Output file name is %s\n", name);
    plot(mat_ungapped, new_length, n_seqs, name);

    sprintf(name, "%s_%s_%s_conservation_ungapped.png", get_basename(av[1], strlen(av[1])), av[2], av[3]);
    fprintf(stdout, "Output file name is %s\n", name);
    plot(mat_nucl_canvas, new_length, MAX_PERCENTAGE, name);

    //end = clock();

    // data_database.total_len = pos_in_database;

    //fprintf(stdout, "[INFO] Database loaded and of length %"PRIu64". Hash table building took %e seconds\n", data_database.total_len, (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(database);


    /*free_mat((void ***) &mat, 3 * file_length);
    free_mat((void ***) &mat_ungapped, 3 * file_length);
    free_mat((void ***) &mat_nucl, file_length);
    free_mat((void ***) &mat_nucl_canvas, MAX_PERCENTAGE + 1);
    */
    /* VARIABLES ROCKING IN A FREE WORLD */
	int i;
	for(i=0; i<n_seqs; i++){
		free(mat[i]);
		free(mat_ungapped[i]);
		free(mat_nucl[i]);
	}
	free(mat);
	free(mat_ungapped);
	free(mat_nucl);

	for(i=0; i<MAX_PERCENTAGE+1; i++){
		free(mat_nucl_canvas[i]);
	}
	free(mat_nucl_canvas);

    free(conservation);
    free(conservation_smooth);
    free(conservation_ungapped);
    free(conservation_trace);

    return 0;
}

void free_mat(void *** mat, int width){
	int i;
	for(i=0; i<width; i++){
		free((*mat)[i]);
	}
	free(mat);
}

void convolve(float * conservation, float * conservation_smooth, int width, int window_size) {
    printf("Convolving conservation with %d of windows size\n", window_size);
    int i, j;
    for (i = 0; i < width; i++) {
        if ((i < window_size) || (i > width - window_size)) {
            conservation_smooth[i] = conservation[i];
        } else {
            float x = 0.0;
            for (j = -window_size; j < window_size; j++) {
                x += conservation[i + j];
            }
            conservation_smooth[i] = x / (window_size * 2);
        }
        //printf(" i:%d n:%f s:%f |",i , altura[i], alturasmooth[i]);
    }
}

int remove_gaps(float ** mat, float ** mat_ungapped, int wide, int height, int percentage, float * conservation, float * conservation_ungapped) {
    printf("Ungapping at %d%%\n", percentage);
    int u_width = 0;
    int gapps_count = 0;

    for (int w = 0; w < wide; w++) {
        // get gapps count
        gapps_count = 0;
        for (int h = 0; h < height; h++) {
            if ((mat[h][3 * w] == 0.0) && (mat[h][3 * w + 1] == 0.0) && (mat[h][3 * w + 2] == 0.0)) {
                gapps_count++;
            }
        }
        // Check if 90% are gapps, if it isn't then add it to the new mat and increase the new width
        if ((gapps_count * 100) / height < percentage) {
            for (int hu = 0; hu < height; hu++) {
                mat_ungapped[hu][3 * u_width] = mat[hu][3 * w];
                mat_ungapped[hu][3 * u_width + 1] = mat[hu][3 * w + 1];
                mat_ungapped[hu][3 * u_width + 2] = mat[hu][3 * w + 2];
            }
            conservation_ungapped[u_width] = conservation[w];
            u_width++;
        }
    }

    return u_width;
}

void plot(float ** mat, int wide, int height, char * name) {
    int size_h = wide;
    int size_v = height;
    int multiplier_v = 1;
    char stmp[256];
    int borde_h = 50;
    int borde_v = 50;

    // Add some pixels if the height is too small
    while (multiplier_v * size_v < MIN_HEIGHT) {
        multiplier_v++;
    }

    pngwriter png(size_h + 2 * borde_h, multiplier_v * size_v + 2 * borde_v, 0, name);

    fprintf(stdout, "Plotting using %d for (%d,%d)\n", multiplier_v, size_h, size_v);

    int i, j, k;
    for (i = 0; i < size_v; i++) {
        for (j = 0; j < size_h; j++) {
            if ((i * 100 / size_v) % 5 == 0) printf("%d%%\r", (i * 100 / size_v));
            for (k = 0; k < multiplier_v; k++) {
                png.plot(j + borde_h, i * multiplier_v + k + borde_v, mat[i][3 * j], mat[i][3 * j + 1], mat[i][3 * j + 2]);
            }
        }
    }
    png.close();
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

inline void get_color(char a, float * c1, float * c2, float * c3) {
    if (a == 'A' || a == 'a') {
        * c1 = 1.0;* c2 = 0.0;* c3 = 0.0;
        return;
    }
    if (a == 'C' || a == 'c') {
        * c1 = 0.0;* c2 = 1.0;* c3 = 0.0;
        return;
    }
    if (a == 'G' || a == 'g') {
        * c1 = 0.0;* c2 = 0.0;* c3 = 1.0;
        return;
    }
    if (a == 'T' || a == 't') {
        * c1 = 1.0;* c2 = 1.0;* c3 = 0.0;
        return;
    }
    if (a == 'N' || a == 'n') {
        * c1 = 1.0;* c2 = 1.0;* c3 = 1.0;
        return;
    }
    * c1 = 0.0;* c2 = 0.0;* c3 = 0.0;
}

char * get_basename(char * a, int file_length) {
    int p = file_length - 1;
    while (a[p] != '/' && p >= 0) p--;
    return &a[p + 1];
}