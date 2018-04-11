NAME:=clebsch

SRC:=../..
FLAGS:=-std=c++14 -g # put -O3 here when ready

INCLUDE_PATH:=-I$(SRC) -I$(SRC)/json/single_include -I$(SRC)/fifo_map/src
LIBS:=-lm -lgsl -lboost_filesystem -lboost_system

CC:=g++ -c $(FLAGS) $(INCLUDE_PATH)
LINK:=g++ $(FLAGS) $(INCLUDE_PATH) $(LIBS)

SRCMAKE:=make -C $(SRC)
DEPS:=utils.o
DEPS:=$(addprefix $(SRC)/,$(DEPS))

.PHONY: all clean $(DEPS)

all: $(NAME).exe

$(NAME).exe:$(NAME).cpp $(addprefix $(SRC)/, utils.o)
	$(LINK) -o $(NAME).exe $^

$(DEPS):
	$(SRCMAKE) $(patsubst $(SRC)/%,%,$@)
	
clean:
	rm -rf ./output/
	rm -f $(NAME).exe
