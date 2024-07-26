CC = gcc
CPP = g++
CFLAGS = -Wall -g3 -Ofast

all: k2rasterProg
	
k2rasterProg:
	$(CPP) -std=c++11 -I ./libs/libcds2/includes -I ./include k2r.c utils.c ./no_restrictions/dacs.c ./no_restrictions/bitrankw32int.c ./no_restrictions/basics.c -o k2r $(CFLAGS)

clean:
	rm k2r utils ./no_restrictions/dacs ./no_restrictions/bitrankw32int ./no_restrictions/basics
