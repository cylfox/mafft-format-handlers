all:
	g++ -g mafft2plotComplete.c pngwriter/src/pngwriter.h -o mafft2plot `freetype-config --cflags` -l:libPNGwriter.a -lm -lpng -lz -lfreetype
	g++ mafftBorderLocator.c -o mafftBorderLocator
	g++ mafftBorderExtractor.c -o mafftBorderExtractor
