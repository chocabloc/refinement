all: test

test: test.cpp
	g++ -I. -O3 test.cpp -o test

run: test
	./test