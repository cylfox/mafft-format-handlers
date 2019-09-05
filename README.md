# Mafft Format Handlers
Some scripts to handle mafft format files in C.
## Author
- EPW <estebanpw@uma.es>
- MSR <marcossr@uma.es>
## mafft2plotComplete
### Description
Plots a mafft alignment with and without gapps with colors (even for gigabyte-sized files) and their own conservation plots (4 png images altogether)
### Requirements
You'll need libpng to be installed. You can obtain it by running "sudo apt-get install libpng-dev"
### Makefile
```c
g++ -g mafft2plotComplete.c pngwriter/src/pngwriter.h -o mafft2plot `freetype-config --cflags` -l:libPNGwriter.a -lm -lpng -lz -lfreetype
```
### Usage
```c
./mafft2plot <mafft alignment file>
```
---
## mafftBorderLocator
### Description
Locates the borders of a alignment file in mafft format using as parameters the conservation percentage and the maximum number of consecutive hits to identify each border.
### Makefile
```c
g++ mafftBorderLocator.c -o mafftBorderLocator
```
### Usage
```c
./mafftBorderLocator <conservation percentage> <max consecutive hits> <mafft alignment file>
```
---
## mafftBorderExtractor
### Description
Extracts the borders of a alignment file in mafft format using as parameters the first and second positions in the sequence.
### Makefile
```c
g++ mafftBorderExtractor.c -o mafftBorderExtractor
```
### Usage
```c
./mafftBorderExtractor <first cut> <second cut> <mafft alignment file>
```
---
