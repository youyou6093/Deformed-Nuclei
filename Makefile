main: main.o RacahFunctions.o
	g++ main.o -fopenmp RacahFunctions.o -o main -lgsl -lblas  

main.o:
	g++ -c -fopenmp -std=c++11 main.cpp  -O3

RacahFunctions.o:
	gfortran RacahFunctions.f -c

clean:
	rm *.o main
