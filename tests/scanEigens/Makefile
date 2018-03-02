NAME=scanEigens

SRC=../..
FLAGS=-std=c++14 -O3 # put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3 -I$(SRC)
LIBS=-lm

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

EXTMAKE=make -C $(SRC)
DEPS=env_deng2016lin.o 

.PHONY: all clean wipe

all: $(DEPS)
	$(LINK) $(addprefix $(SRC)/, $(DEPS)) -o $(NAME) $(NAME).cpp

run:
	./$(NAME)

$(DEPS):
	$(EXTMAKE) $@
	
clean:
	rm -f *.o

wipe: clean
	rm -f *.dat
	rm -f $(NAME)
