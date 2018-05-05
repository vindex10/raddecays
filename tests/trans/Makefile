NAME:=trans

SRC:=../..
FLAGS:=-std=c++14 -O3 # put -O3 here when ready

INCLUDE_PATH:=-I/usr/include/eigen3 -I$(SRC) -I$(SRC)/json/single_include -I$(SRC)/cubature -I$(SRC)/fifo_map/src
LIBS:=-L$(HOME)/work/build/lib -lm -lgsl -lboost_filesystem -lboost_system

CC:=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK:=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

SRCMAKE:=make -C $(SRC)
DEPS:=env_deng2016lin.o env_deng2016scr.o cubature/hcubature.o utils.o
DEPS:=$(addprefix $(SRC)/,$(DEPS))

TRANSITIONS:=$(basename $(notdir $(wildcard transitions/*)))

.PHONY: all prepare generate clean $(TRANSITIONS)

all: $(NAME)Lin.exe $(NAME)Scr.exe

generate: prepare $(TRANSITIONS)

prepare:
	make -C ../quarkU/ generate

$(filter %-lin, $(TRANSITIONS)): $(NAME)Lin.exe
	./$(NAME)Lin.exe transitions/$@.cfg

$(filter %-scr, $(TRANSITIONS)): $(NAME)Scr.exe
	./$(NAME)Scr.exe transitions/$@.cfg

$(NAME)Lin.exe:$(NAME).cpp $(addprefix $(SRC)/, env_deng2016lin.o cubature/hcubature.o utils.o)
	$(SRCMAKE) $(filter %.o, $(patsubst $(SRC)/%,%,$^))
	$(LINK) -include $(NAME)Lin.hpp -o $(NAME)Lin.exe $^

$(NAME)Scr.exe:$(NAME).cpp $(addprefix $(SRC)/, env_deng2016scr.o cubature/hcubature.o utils.o)
	$(SRCMAKE) $(filter %.o, $(patsubst $(SRC)/%,%,$^))
	$(LINK) -include $(NAME)Scr.hpp  -o $(NAME)Scr.exe $^

$(DEPS):
	$(SRCMAKE) $(patsubst $(SRC)/%,%,$<)

clean:
	rm -rf ./output/
	rm -f $(NAME){Lin,Scr}.exe
