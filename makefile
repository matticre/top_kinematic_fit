CC = g++
CFLAGS = `root-config --cflags`
CLIBS  = `root-config --glibs`
PDG    = -lEG
MIN	   = -lMinuit

all: main

main: main.cpp
	$(CC) -o main main.cpp $(CFLAGS) $(CLIBS) $(PDG) $(MIN)
clean:
	rm *.o
