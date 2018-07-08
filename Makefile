all: all_$(NAME)

$(BUILD)/utils.o: $(addprefix $(SRC)/,utils.hpp utils.cpp) $$(DEPREQS)
	$(PRECOMPILE)
	$(CC) -o $@ $(SRC)/utils.cpp
	$(POSTCOMPILE)

$(BUILD)/env_deng2016lin.o: $(addprefix $(SRC)/,env_deng2016lin.cpp) $$(DEPREQS)
	$(PRECOMPILE)
	$(CC) -o $@ $(SRC)/env_deng2016lin.cpp
	$(POSTCOMPILE)

$(BUILD)/env_deng2016scr.o: $(addprefix $(SRC)/,env_deng2016scr.cpp) $$(DEPREQS)
	$(PRECOMPILE)
	$(CC) -o $@ $(SRC)/env_deng2016scr.cpp
	$(POSTCOMPILE)

$(BUILD)/observers.o: $(addprefix $(SRC)/,observers.cpp) $$(DEPREQS)
	$(PRECOMPILE)
	$(CC) -o $@ $(SRC)/observers.cpp
	$(POSTCOMPILE)

$(BUILD)/hcubature.o:
	make -C $(SRC)/extsrc/ $@

%.d: ;
.PRECIOUS: $(shell find $(SRC) -type f -name '*.d')

include $(shell find $(SRC) -name '*.d' -type f)
