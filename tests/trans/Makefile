NAME:=trans

include ../../MakeVars
include $(SRC)/Makefile

LIBS:=-lm -lgsl -lgslcblas -lboost_filesystem -lboost_system

TRANSITIONS:=$(basename $(notdir $(wildcard transitions/*.cfg)))

.PHONY: all prepare generate export clean $(TRANSITIONS)

all_$(NAME): $(NAME)Lin.exe $(NAME)Scr.exe

report: generate
	jupyter nbconvert --exec $(NAME).ipynb --stdout > /dev/null

generate: prepare $(TRANSITIONS)

prepare:
	make -C ../quarkU/ generate

$(filter %-lin, $(TRANSITIONS)): $(NAME)Lin.exe
	./$(NAME)Lin.exe transitions/$@.cfg

$(filter %-scr, $(TRANSITIONS)): $(NAME)Scr.exe
	./$(NAME)Scr.exe transitions/$@.cfg

$(NAME)Lin.exe:$(NAME).cpp $(addprefix $(BUILD)/, env_deng2016lin.o hcubature.o utils.o) $$(DEPREQS)
	$(PRECOMPILE)
	$(LINK) -include $(NAME)Lin.hpp -o $@ $(NAME).cpp $(addprefix $(BUILD)/, env_deng2016lin.o hcubature.o utils.o)
	$(POSTCOMPILE)

$(NAME)Scr.exe:$(NAME).cpp $(addprefix $(BUILD)/, env_deng2016scr.o hcubature.o utils.o) $$(DEPREQS)
	$(PRECOMPILE)
	$(LINK) -include $(NAME)Scr.hpp  -o $@ $(NAME).cpp $(addprefix $(BUILD)/, env_deng2016scr.o hcubature.o utils.o)
	$(POSTCOMPILE)

export:
	jupyter nbconvert --to html --execute  $(NAME).ipynb
	mv $(NAME).html export/

clean:
	$(RM) ./output/
	$(RM) ./export/*
	$(RM) $(NAME){Lin,Scr}.exe
	$(RM) $(DEPDIRNAME)
