
CXXFLAGS +=  --std=c++0x -O2

.PHONY : all test clean

all : test

test : 
	g++ $(CXXFLAGS) *.cc test/*.cc -I. -lboost_unit_test_framework -o test_runner 
	./test_runner

clean : 
	rm -rf test_runner
