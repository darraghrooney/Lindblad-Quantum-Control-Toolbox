IDIR =../include
CC=g++
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS= -l armadillo

_DEPS = 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = Central.o Menus.o LSystem2d.o Numerics.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: prep Lindblad2

prep:
	mkdir -p $(ODIR)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c $^ -o $@ $(CFLAGS)

Lindblad2: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -r -f $(ODIR)
	rm -f $(IDIR)/*~ *.dat *.gnu *~ core 

