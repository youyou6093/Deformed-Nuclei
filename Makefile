main: main.o RacahFunctions.o
	g++ main.o RacahFunctions.o -o main -lgsl -lblas  

main.o:
	g++ -c -std=c++11 main.cpp -O3

clean:
	rm main.o main