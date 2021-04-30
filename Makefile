MAIN = edge_main
VERT = vert_main
SUB = sub_main
EXPPL = exp_powerlaw
OBJS = part.o load.o FM.o FMvert.o FMvert2.o clustermap.o ordering.o

CXX = g++
CXXFLAGS = -O3 -Wall

all: $(MAIN) $(VERT) $(SUB) $(EXPPL)

$(MAIN): $(MAIN).o $(OBJS)
	$(CXX) -o $@ $^

$(VERT): $(VERT).o $(OBJS)
	$(CXX) -o $@ $^

$(SUB): $(SUB).o $(OBJS)
	$(CXX) -o $@ $^

$(EXPPL): $(EXPPL).o $(OBJS)
	$(CXX) -o $@ $^

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o $(MAIN) $(VERT)
