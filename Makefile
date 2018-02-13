FLAGS=-std=c++14 # put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3
LIBS=-lgsl -lm

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

.PHONY: all clean 

clebsh.o: clebsh.cpp 
	$(CC) -o clebsh.o clebsh.cpp

env_deng2016lin.o: env_deng2016lin.cpp env_deng2016lin.hpp
	$(CC) -o env_deng2016lin.o env_deng2016lin.cpp

env_deng2016scr.o: env_deng2016scr.cpp env_deng2016scr.hpp
	$(CC) -o env_deng2016scr.o env_deng2016scr.cpp

observers.o: observers.cpp observers.hpp
	$(CC) -o observers.o observers.cpp

clean:
	rm -rf $(OUTPUT_DIR)
	rm -f *.o
