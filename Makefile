FLAGS=-std=c++14 # put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3
LIBS=-lgsl -lm

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

OUTPUT_DIR=build

.PHONY: all clebsh clean 

all: clebsh
	mkdir -p $(OUTPUT_DIR)
	$(LINK) -o $(OUTPUT_DIR)/main\
			clebsh.o main.cpp

clebsh: clebsh.o
clebsh.o: clebsh.cpp 
	$(CC) -o clebsh.o clebsh.cpp

clean:
	rm -rf $(OUTPUT_DIR)
	rm -f *.o
