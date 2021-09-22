test: posit.cc test_posit.cc ../../../test/gtest/gtest-all.cc
	g++ -std=c++11 -Wall -g -pthread -I../../../include -I../../../build/include -I../../../test $^ -o $@ 

debug: posit.cc test_posit.cc ../../../test/gtest/gtest-all.cc
	g++ -std=c++11 -Wall -g -pthread -I../../../include -I../../../build/include -I../../../test -D PDEBUG $^ -o $@ 