CXX = g++
CXXFlags = -Wall -Wextra -pedantic -g -O2 -DNDEBUG

IDIR =$(HOME)/include
LIBDIR=$(HOME)/lib

EXEC =	RedshiftBiasEFT
SRC = $(wildcard *.cpp)
OBJ = $(SRC: .cpp=.o) 
STD = -std=c++14
GSL = -lm -lgsl -lgslcblas
CUBA = -I$(IDIR) -L$(LIBDIR) -lcuba -lm
FFTW = -lfftw3


all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $^ -o $@ $(STD) $(GSL) $(FFTW)

%.o: %.cpp
	$(CXX) $(CXXFlags) -o $@ -c $^ 

.PHONY: clean

clean:
	rm $(EXEC) *~

