CC      =g++ 
CXXFLAGS  =-g -std=c++11 -Wall -pedantic -O2
LDFLAGS =-lncurses -pthread -lz 

OBJ = main.o DNA.o fastq.o  fqless.o

fqless: $(OBJ)
		$(CC) $(CFLAGS) -o fqless $(OBJ) $(LDFLAGS)

%.o: %.c
		$(CC) $(CFLAGS) -c $<


