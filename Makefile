MAIN = main
OBJS = part.o load.o FM.o

CXX = g++
CXXFLAGS = -O3 -Wall

all: $(MAIN)

$(MAIN): $(MAIN).o $(OBJS)
	$(CXX) -o $@ $^

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o $(MAIN)
