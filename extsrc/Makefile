include ../MakeVars
EXTSRC:=$(SRC)/extsrc

.PHONY: all

all: $(addprefix $(BUILD)/, include/nlohmann/json.hpp include/fifo_map.hpp lib64/libnlopt.so hcubature.o)

$(BUILD)/include/nlohmann/json.hpp:
	$(MKDIR) $(BUILD)/include
	$(CP) $(EXTSRC)/json/single_include/nlohmann $(BUILD)/include

$(BUILD)/include/fifo_map.hpp:
	$(MKDIR) $(BUILD)/include
	$(CP) $(EXTSRC)/fifo_map/src/fifo_map.hpp $@

$(BUILD)/lib64/libnlopt.so:
	$(MKDIR) $(EXTSRC)/nlopt/build
	cmake -DCMAKE_INSTALL_PREFIX=$(BUILD) -B$(EXTSRC)/nlopt/build -H$(EXTSRC)/nlopt
	make -C $(EXTSRC)/nlopt/build
	make -C $(EXTSRC)/nlopt/build install

$(BUILD)/hcubature.o:
	$(CC) -o $@ $(EXTSRC)/cubature/hcubature.c
	$(CP) $(EXTSRC)/cubature/cubature.h $(BUILD)/include/cubature.h

