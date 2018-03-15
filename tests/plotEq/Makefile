NAME=plotEq

SRC=../..
FLAGS=-std=c++14 -O3 # put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3 -I$(SRC)
LIBS=-lm -lgsl

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

EXTMAKE=make -C $(SRC)
DEPS=env_deng2016scr.o

.PHONY: all clean wipe

all: $(DEPS)
	$(LINK) $(addprefix $(SRC)/, $(DEPS)) -o $(NAME).exe $(NAME).cpp

$(DEPS):
	$(EXTMAKE) $@
	
clean:
	rm -f *.o
	rm -f *.dat
	rm -f $(NAME).exe
