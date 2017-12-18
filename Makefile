OUTPUT_DIR=build

.PHONY: debug prod clean

debug:
	mkdir -p $(OUTPUT_DIR)
	g++ -std=c++14 -I /usr/include/eigen3 -o $(OUTPUT_DIR)/main  main.cpp
prod: clean
	mkdir -p $(OUTPUT_DIR)
	g++ -std=c++14 -O3 -I /usr/include/eigen3 -o $(OUTPUT_DIR)/main  main.cpp
clean:
	rm -rf $(OUTPUT_DIR)
