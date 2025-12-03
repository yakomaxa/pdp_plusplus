CXX	= g++-13
CXXFLAGS  =  -std=c++20 -Wall -I./gemmi/include/ -O3 -Wl,-ld_classic
#CXXFLAGS  =  -std=c++20 -Wall -I./gemmi/include/ -O0 -Wl,-ld_classic -g
# If you want to complile on linux, consider removing -Wl,-ld_classic; it's for MacOS things.
LIBS    = -lm -lz 
VPATH   = ./src/
SRC     = $(shell ls $(VPATH)/*.cpp)
HEADERS = $(shell ls $(VPATH)/*.hpp) 
OBJS    = $(SRC:.cpp=.o)
TARGET  = pdp++

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

$(OBJS): $(HEADERS)

clean:
	rm -f $(OBJS) $(TARGET)
