# Discussion 10
#   make setup   download the Eigen library (once)
#   make         compile main.cpp
#   make run     compile (if needed) and run
#   make clean   remove the compiled program

EIGEN_VERSION = 3.4.0
EIGEN_DIR     = eigen-$(EIGEN_VERSION)
EIGEN_URL     = https://gitlab.com/libeigen/eigen/-/archive/$(EIGEN_VERSION)/$(EIGEN_DIR).tar.bz2

CXX      = g++
CXXFLAGS = -std=c++17 -Wall -I $(EIGEN_DIR)

.PHONY: all setup run clean

all: main

setup: $(EIGEN_DIR)/Eigen/Dense

$(EIGEN_DIR)/Eigen/Dense:
	curl -fL -o $(EIGEN_DIR).tar.bz2 $(EIGEN_URL)
	tar xjf $(EIGEN_DIR).tar.bz2
	rm -f $(EIGEN_DIR).tar.bz2

main: main.cpp $(EIGEN_DIR)/Eigen/Dense
	$(CXX) $(CXXFLAGS) main.cpp -o main

run: main
	./main

clean:
	rm -f main
