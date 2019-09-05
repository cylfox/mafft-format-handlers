all:
	g++ -g main.cpp pngwriter/src/pngwriter.h -o mafft2plot `freetype-config --cflags` -l:libPNGwriter.a -lm -lpng -lz -lfreetype
	g++ mafftBorderLocator.cpp -o mafftBorderLocator
	g++ mafftBorderExtractor.cpp -o mafftBorderExtractor