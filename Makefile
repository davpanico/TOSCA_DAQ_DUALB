CC = g++
CFLAGS = -g -O -pedantic -Wall
DEPS = keyb.h keyb.c Functions.c Functions.h
EXT_LIBS = $(shell root-config --libs) -lCAENDigitizer
INCLUDES = $(shell root-config --cflags)
FLAGS = -D_DEBUG_

ReadOut: ReadoutTest_Digitizer.cxx $(DEPS)
	$(CC) $^ -o $@ $(CFLAGS) $(FLAGS) $(INCLUDES) $(EXT_LIBS) 

clean: 
	rm -f ReadOut
