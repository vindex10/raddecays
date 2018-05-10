$(BUILD)/utils.o: $(addprefix $(SRC)/,utils.hpp utils.cpp)
	$(CC) -o $@ $(SRC)/utils.cpp

$(BUILD)/env_deng2016lin.o: $(addprefix $(SRC)/,env_deng2016lin.cpp env_deng2016lin.hpp)
	$(CC) -o $@ $(SRC)/env_deng2016lin.cpp

$(BUILD)/env_deng2016scr.o: $(addprefix $(SRC)/,env_deng2016scr.cpp env_deng2016scr.hpp)
	$(CC) -o $@ $(SRC)/env_deng2016scr.cpp

$(BUILD)/observers.o: $(addprefix $(SRC)/,observers.cpp observers.hpp)
	$(CC) -o $@ $(SRC)/observers.cpp

$(BUILD)/hcubature.o:
	make -C $(SRC)/extsrc/ $@

