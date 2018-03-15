NAME=quarkEigen

SRC=../..
FLAGS=-std=c++14 -O3# put -O3 here when ready

INCLUDE_PATH=-I/usr/include/eigen3 -I$(SRC) -I$(SRC)/json/single_include
LIBS=-lm -lgsl -lnlopt_cxx

CC=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

EXTMAKE=make -C $(SRC)

.PHONY: all data clean

all: $(NAME)Lin.exe $(NAME)Scr.exe

generate: data $(NAME)Lin.exe $(NAME)Scr.exe
	./$(NAME)Scr.exe b-scr.cfg
	./$(NAME)Scr.exe c-scr.cfg
	./$(NAME)Lin.exe c-lin.cfg

data:
	../utils/json-merge.py states/b-scr-{S,P,D,F}.cfg > b-scr.cfg
	../utils/json-merge.py states/c-scr-{S,P,D}.cfg > c-scr.cfg
	../utils/json-merge.py states/c-lin-{S,P,D}.cfg > c-lin.cfg

$(NAME)Lin.exe:$(NAME).cpp env_deng2016lin.o
	$(LINK) -include $(NAME)Lin.hpp $(SRC)/env_deng2016lin.o -o $(NAME)Lin.exe $(NAME).cpp

$(NAME)Scr.exe:$(NAME).cpp env_deng2016scr.o
	$(LINK) -include $(NAME)Scr.hpp $(SRC)/env_deng2016scr.o -o $(NAME)Scr.exe $(NAME).cpp

env_deng2016lin.o env_deng2016scr.o:
	$(EXTMAKE) $@

clean:
	rm -rf output
	rm -f *.exe
	rm -f {b,c}-{scr,lin}.cfg
