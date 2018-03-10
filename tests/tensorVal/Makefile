NAME=tensorVal

SRC=../..
FLAGS=-std=c++14 -O3 # put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3 -I$(SRC)
LIBS=-lm -lgsl

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

EXTMAKE=make -C $(SRC)
DEPS=env_deng2016lin.o

.PHONY: all clean wipe

all: $(DEPS)
	$(LINK) $(addprefix $(SRC)/, $(DEPS)) -o $(NAME).exe $(NAME).cpp

$(DEPS):
	$(EXTMAKE) $@
	
run:
	./$(NAME)

clean:
	rm -f *.o

wipe: clean
	rm -f *.dat
	rm -f $(NAME).exe
