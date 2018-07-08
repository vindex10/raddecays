include ../../MakeVars

NAME:=quarkEigen

LIBS:=-lm -lgsl -lgslcblas -lnlopt -lboost_system -lboost_filesystem

STATES:=$(basename $(notdir $(wildcard states/*)))

.PHONY: all generate clean $(STATES)

all_$(NAME): $(NAME)Lin.exe $(NAME)Scr.exe

generate: $(STATES)

$(filter %-lin, $(STATES)): $(NAME)Lin.exe
	./$(NAME)Lin.exe states/$@.cfg

$(filter %-scr, $(STATES)): $(NAME)Scr.exe
	./$(NAME)Scr.exe states/$@.cfg

$(NAME)Lin.exe:$(NAME).cpp $(BUILD)/env_deng2016lin.o $$(DEPREQS)
	$(PRECOMPILE)
	$(LINK) -include $(NAME)Lin.hpp $(BUILD)/env_deng2016lin.o -o $(NAME)Lin.exe $(NAME).cpp
	$(POSTCOMPILE)

$(NAME)Scr.exe:$(NAME).cpp $(BUILD)/env_deng2016scr.o $$(DEPREQS)
	$(PRECOMPILE)
	$(LINK) -include $(NAME)Scr.hpp $(BUILD)/env_deng2016scr.o -o $(NAME)Scr.exe $(NAME).cpp
	$(POSTCOMPILE)

clean:
	$(RM) ./output/
	$(RM) $(NAME){Lin,Scr}.exe
	$(RM) $(DEPDIRNAME)

include $(SRC)/Makefile
