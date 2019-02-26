main: main.o RacahFunctions.o
	g++-8 -L/usr/local/lib main.o -fopenmp RacahFunctions.o -o main -lgsl -lblas  

main.o:
	g++-8 -I/usr/local/include -c -fopenmp -std=c++11 main.cpp  -O3

RacahFunctions.o:
	gfortran RacahFunctions.f -c

clean:
	rm *.o main
