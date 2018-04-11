FLAGS:=-std=c++14 -O3 # put -O3 here when ready

INCLUDE_PATH:=-I/usr/include/eigen3 -I./json/single_include -I./fifo_map/src
LIBS:=-lm -lgsl -lboost_filesystem -lboost_system

CC:=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK:=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

.PHONY: all clean

utils.o: utils.hpp utils.cpp 
	$(CC) -o utils.o utils.cpp

env_deng2016lin.o: env_deng2016lin.cpp env_deng2016lin.hpp
	$(CC) -o env_deng2016lin.o env_deng2016lin.cpp

env_deng2016scr.o: env_deng2016scr.cpp env_deng2016scr.hpp
	$(CC) -o env_deng2016scr.o env_deng2016scr.cpp

observers.o: observers.cpp observers.hpp
	$(CC) -o observers.o observers.cpp

cubature/hcubature.o: cubature/cubature.h cubature/hcubature.c
	$(CC) -o cubature/hcubature.o cubature/hcubature.c

clean:
	rm -f *.o
	rm -f cubature/hcubature.o
