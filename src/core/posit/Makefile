test: posit.cc test_posit.cc ../../../test/gtest/gtest-all.cc posit_helper.cc
	g++ -std=c++11 -Wall -g -pthread -I../../../include -I../../../build/include -I../../../test $^ -o $@ 

debug: posit.cc test_posit.cc ../../../test/gtest/gtest-all.cc posit_helper.cc
	g++ -std=c++11 -Wall -g -pthread -I../../../include -I../../../build/include -I../../../test -D PDEBUG $^ -o $@ 


posit2: posit2.cc test_posit2.cc ../../../test/gtest/gtest-all.cc
	g++ -std=c++11 -Wall -g -pthread -I../../../build/include -I../../../test -D PDEBUG $^ -o $@ 