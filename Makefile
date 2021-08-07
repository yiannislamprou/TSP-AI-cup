# Makefile
# AI cup 2012

HEADER = TSP.h
SOURCE = TSP.cpp main.cpp
OBJ = TSP.o main.o
OUT = tsp

all: $(OBJ)
	g++ -o $(OUT) $(OBJ)

TSP.o: TSP.cpp
		g++ -c TSP.cpp

main.o: main.cpp
		g++ -c main.cpp		

count: 
	wc -wl $(SOURCE) $(HEADER)

clean: 
	rm -f $(OBJ) $(OUT)		