NAME=sphharm_eigenval
PRJ_DIR=../..

INCLUDE_PATH=-I$(PRJ_DIR)
LIBS=-lm -lgsl

LINK=g++ -std=c++14 $(INCLUDE_PATH) $(LIBS)

all:
	make -C $(PRJ_DIR) clebsh.o
	$(LINK) -o $(NAME).exe $(PRJ_DIR)/clebsh.o $(NAME).cpp

clean:
	rm -f $(NAME).exe
