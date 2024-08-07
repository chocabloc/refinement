all: test

test: test.cpp objgradfun.h
	g++ -I. -O3 test.cpp -o test

run: test
	./test