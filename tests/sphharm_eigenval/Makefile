NAME=sphharm_eigenval
PRJ_DIR=../..

INCLUDE_PATH=-I$(PRJ_DIR)
LIBS=-lm -lgsl

LINK=g++ -std=c++14 $(INCLUDE_PATH) $(LIBS)

OUTPUT_DIR=build

all:
	make -C $(PRJ_DIR) clebsh
	mkdir -p $(OUTPUT_DIR)
	$(LINK) -o $(OUTPUT_DIR)/$(NAME) $(PRJ_DIR)/clebsh.o $(NAME).cpp

clean:
	rm -rf $(OUTPUT_DIR)
