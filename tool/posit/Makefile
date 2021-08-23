CXX := g++
LD := g++

INCLUDE 	:= .
EXE		:= cnn

SRCS    	:= $(wildcard *.cc)
OBJS    	:= $(patsubst %.cc,%.o,$(SRCS))
DATA_OBJS 	:= images.o labels.o myfilesnap.o
OBJS		:= $(DATA_OBJS) $(OBJS)

CFLAGS  := -I$(INCLUDE) -DLITE_POSIT -g

all: $(EXE)

$(EXE): $(OBJS)
	$(LD) $^ -o $@
    
%.o: %.cc
	$(CXX) $(CFLAGS) -c $< -o $@
