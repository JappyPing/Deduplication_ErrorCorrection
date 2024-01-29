sources=Sequence.cpp Utils.cpp Options.cpp SeqFileParser.cpp PairedSeqFileParser.cpp Cluster.cpp ClusterPaired.cpp \
	SingleEnd.cpp PairedEnd.cpp main.cpp

root_dir = .
core_dir = .
objs_dir = objs
objs += $(patsubst %.cpp, $(objs_dir)/%.cpp.o, $(sources))
mylibs = -lpthread -lm -lz -lgomp
EXEC = ParDRe

#compile optionsK
CXXFLAGS = -I . -pthread -std=c++11 -O3 -funroll-loops -ffast-math -Wall
MYLIBS += -lm -lz
CXX = mpic++

all: dir $(objs)
	$(CXX) $(objs) -o $(EXEC) $(MYLIBS) $(CXXFLAGS)
	strip $(EXEC)

dir:
	mkdir -p $(objs_dir)

clean:
	-rm -rf $(objs_dir) $(EXEC)

$(objs_dir)/%.cpp.o: $(core_dir)/%.cpp
	$(CXX) $(MYLIBS) -o $@ -c $< $(CXXFLAGS)
	
test:
	prove t/*.t
