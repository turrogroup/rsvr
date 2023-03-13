CXX=g++
CXXFLAGS=-I. -std=c++17

rsvr : rsvr.o enc.o idtypes.o seqfx.o utils.o
	$(CXX) $(CXXFLAGS) -o rsvr rsvr.o enc.o idtypes.o seqfx.o utils.o

clean:
	rm *.o

install:
	cp rsvr /usr/local/bin
	cp scripts/* /usr/local/bin
