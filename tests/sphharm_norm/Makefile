NAME=sphharm_norm
PRJ_DIR=../..

INCLUDE_PATH=-I$(PRJ_DIR) -I$(PRJ_DIR)/cubature
LIBS=-lm -lgsl

LINK=g++ -std=c++14 $(INCLUDE_PATH) $(LIBS)

all:
	make -C $(PRJ_DIR) clebsh.o
	$(LINK) -o $(NAME).exe\
		$(PRJ_DIR)/clebsh.o $(PRJ_DIR)/cubature/hcubature.c $(NAME).cpp

clean:
	rm -f $(NAME).exe
