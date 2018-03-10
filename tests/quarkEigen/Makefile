NAME=quarkEigen

SRC=../..
FLAGS=-std=c++14 -O3# put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3 -I$(SRC) -I$(SRC)/json/single_include
LIBS=-lm -lgsl -lnlopt_cxx

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

EXTMAKE=make -C $(SRC)
DEPS=env_deng2016lin.o env_deng2016scr.o

.PHONY: all clean wipe

all: $(NAME)Lin.exe $(NAME)Scr.exe

$(DEPS):
	$(EXTMAKE) $@

$(NAME)Lin.exe:$(NAME).cpp env_deng2016lin.o
	$(LINK) -include $(NAME)Lin.hpp $(SRC)/env_deng2016lin.o -o $(NAME)Lin.exe $(NAME).cpp

$(NAME)Scr.exe:$(NAME).cpp env_deng2016scr.o
	$(LINK) -include $(NAME)Scr.hpp $(SRC)/env_deng2016scr.o -o $(NAME)Scr.exe $(NAME).cpp
	
clean:
	rm -f *.dat
	rm -rf output
	rm -f *.exe
