mode := debug

all: compile link run

compile:
	g++ -c src/*.cpp -std=c++17 -g -Wall -m64 -I include

link:
	g++ -o bin/$(mode)/main *.o -lmingw32

run:
	.\bin\$(mode)\main