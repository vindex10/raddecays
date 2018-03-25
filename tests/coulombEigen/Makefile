NAME=coulombEigen

SRC=../..
FLAGS=-std=c++14 -O3# put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3 -I$(SRC) -I$(SRC)/json/single_include -I$(SRC)/fifo_map/src
LIBS=-lm -lgsl -lboost_system -lboost_filesystem

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

EXTMAKE=make -C $(SRC)
DEPS=env_deng2016lin.o

.PHONY: all clean

all: $(DEPS)
	$(LINK) $(addprefix $(SRC)/, $(DEPS)) -o $(NAME).exe $(NAME).cpp

$(DEPS):
	$(EXTMAKE) $@
	
clean:
	rm -rf output/
	rm -f $(NAME).exe
