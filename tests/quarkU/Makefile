NAME:=quarkU

SRC:=../..
FLAGS:=-std=c++14 -O3 # put -O3 here when ready

INCLUDE_PATH:=-I/usr/include/eigen3 -I$(SRC) -I$(SRC)/json/single_include -I$(SRC)/fifo_map/src
LIBS:=-lm -lgsl -lboost_filesystem -lboost_system

CC:=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK:=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

SRCMAKE:=make -C $(SRC)
DEPS:=env_deng2016lin.o env_deng2016scr.o observers.o
DEPS:=$(addprefix $(SRC)/,$(DEPS))

STATES=$(basename $(notdir $(wildcard states/*)))

.PHONY: all generate prepare clean $(STATES)

all: $(NAME)Lin.exe $(NAME)Scr.exe

generate: prepare $(STATES)

$(filter %-lin, $(STATES)): $(NAME)Lin.exe
	./$(NAME)Lin.exe states/$@.cfg

$(filter %-scr, $(STATES)): $(NAME)Scr.exe
	./$(NAME)Scr.exe states/$@.cfg

prepare:
	make -C ../quarkEigen/ generate

$(NAME)Lin.exe:$(NAME).cpp $(addprefix $(SRC)/, env_deng2016lin.o observers.o)
	$(SRCMAKE) $(filter %.o, $(patsubst $(SRC)/%,%,$^))
	$(LINK) -include $(NAME)Lin.hpp $(SRC)/env_deng2016lin.o $(SRC)/observers.o -o $(NAME)Lin.exe $(NAME).cpp

$(NAME)Scr.exe:$(NAME).cpp $(SRC)/env_deng2016scr.o $(SRC)/observers.o
	$(SRCMAKE) $(filter %.o, $(patsubst $(SRC)/%,%,$^))
	$(LINK) -include $(NAME)Scr.hpp $(SRC)/env_deng2016scr.o $(SRC)/observers.o -o $(NAME)Scr.exe $(NAME).cpp

$(DEPS):
	$(SRCMAKE) $(patsubst $(SRC)/%,%,$<)
	
clean:
	rm -rf ./output/
	rm -f $(NAME){Lin,Scr}.exe
