CXX=g++
CXFLAGS=-O2 -std=c++11
LDFLAGS=
EXEC=simu

SRC=Simulator.cpp Pattern.cpp main.cpp
HEADER=Simulator.hpp Pattern.hpp
OBJ=$(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXFLAGS) -o $@ -c $< $(LDFLAGS)

.PHONY: clean cleanall

clean:
	rm -rf $(OBJ)

cleanall: clean
	rm -rf $(EXEC)
