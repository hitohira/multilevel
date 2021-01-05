MAIN = edge_main
VERT = vert_main
OBJS = part.o load.o FM.o FMvert.o FMvert2.o

CXX = g++
CXXFLAGS = -O3 -Wall

all: $(MAIN) $(VERT)

$(MAIN): $(MAIN).o $(OBJS)
	$(CXX) -o $@ $^

$(VERT): $(VERT).o $(OBJS)
	$(CXX) -o $@ $^


.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o $(MAIN) $(VERT)
