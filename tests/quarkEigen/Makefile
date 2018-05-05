NAME:=quarkEigen

SRC:=../..
FLAGS:=-std=c++14 -O3# put -O3 here when ready

INCLUDE_PATH:=-I/usr/include/eigen3 -I$(SRC) -I$(SRC)/json/single_include -I$(SRC)/fifo_map/src
LIBS:=-lm -lgsl -lnlopt -lboost_system -lboost_filesystem

CC:=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK:=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

SRCMAKE:=make -C $(SRC)
DEPS:=env_deng2016lin.o env_deng2016scr.o
DEPS:=$(addprefix $(SRC)/,$(DEPS))

STATES:=$(basename $(notdir $(wildcard states/*)))

.PHONY: all generate clean $(STATES)

all: $(NAME)Lin.exe $(NAME)Scr.exe

generate: $(STATES)

$(filter %-lin, $(STATES)): $(NAME)Lin.exe
	./$(NAME)Lin.exe states/$@.cfg

$(filter %-scr, $(STATES)): $(NAME)Scr.exe
	./$(NAME)Scr.exe states/$@.cfg

$(NAME)Lin.exe:$(NAME).cpp $(SRC)/env_deng2016lin.o
	$(SRCMAKE) $(filter %.o, $(patsubst $(SRC)/%, %, $^))
	$(LINK) -include $(NAME)Lin.hpp $(SRC)/env_deng2016lin.o -o $(NAME)Lin.exe $(NAME).cpp

$(NAME)Scr.exe:$(NAME).cpp $(SRC)/env_deng2016scr.o
	$(SRCMAKE) $(filter %.o, $(patsubst $(SRC)/%, %, $^))
	$(LINK) -include $(NAME)Scr.hpp $(SRC)/env_deng2016scr.o -o $(NAME)Scr.exe $(NAME).cpp

clean:
	rm -rf ./output/
	rm -f $(NAME){Lin,Scr}.exe
