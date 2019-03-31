include ../../MakeVars
include $(SRC)/Makefile

NAME:=quarkU

LIBS:=-lm -lgsl -lgslcblas -lboost_filesystem -lboost_system

STATES=$(basename $(notdir $(wildcard states/*)))

.PHONY: all generate prepare export clean $(STATES)

all: $(NAME)Lin.exe $(NAME)Scr.exe

generate: prepare $(STATES)

$(filter %-lin, $(STATES)): $(NAME)Lin.exe
	./$(NAME)Lin.exe states/$@.cfg

$(filter %-scr, $(STATES)): $(NAME)Scr.exe
	./$(NAME)Scr.exe states/$@.cfg

prepare:
	make -C ../quarkEigen/ generate

$(NAME)Lin.exe:$(NAME).cpp $(addprefix $(BUILD)/, env_deng2016lin.o observers.o) $$(DEPREQS)
	$(PRECOMPILE)
	$(LINK) -include $(NAME)Lin.hpp $(BUILD)/env_deng2016lin.o $(BUILD)/observers.o -o $@ $(NAME).cpp
	$(POSTCOMPILE)

$(NAME)Scr.exe:$(NAME).cpp $(addprefix $(BUILD)/, env_deng2016scr.o observers.o) $$(DEPREQS)
	$(PRECOMPILE)
	$(LINK) -include $(NAME)Scr.hpp $(BUILD)/env_deng2016scr.o $(BUILD)/observers.o -o $@ $(NAME).cpp
	$(POSTCOMPILE)

export:
	jupyter nbconvert --to html --execute  $(NAME).ipynb
	mv $(NAME).html export/

clean:
	$(RM) ./output/
	$(RM) $(NAME){Lin,Scr}.exe
	$(RM) $(DEPDIRNAME)
