NAME:=trans

include ../../MakeVars
include $(SRC)/Makefile

LIBS:=-lm -lgsl -lgslcblas -lboost_filesystem -lboost_system
LINK:=$(LINK) $(LIBS)

TRANSITIONS:=$(basename $(notdir $(wildcard transitions/*)))

.PHONY: all prepare generate clean $(TRANSITIONS)

all_$(NAME): $(NAME)Lin.exe $(NAME)Scr.exe

generate: prepare $(TRANSITIONS)

prepare:
	make -C ../quarkU/ generate

$(filter %-lin, $(TRANSITIONS)): $(NAME)Lin.exe
	./$(NAME)Lin.exe transitions/$@.cfg

$(filter %-scr, $(TRANSITIONS)): $(NAME)Scr.exe
	./$(NAME)Scr.exe transitions/$@.cfg

$(NAME)Lin.exe:$(NAME).cpp $(addprefix $(BUILD)/, env_deng2016lin.o hcubature.o utils.o)
	$(LINK) -include $(NAME)Lin.hpp -o $@ $^

$(NAME)Scr.exe:$(NAME).cpp $(addprefix $(BUILD)/, env_deng2016scr.o hcubature.o utils.o)
	$(LINK) -include $(NAME)Scr.hpp  -o $@ $^

clean:
	$(RM) ./output/
	$(RM) $(NAME){Lin,Scr}.exe
