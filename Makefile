CC      =g++ 
CXXFLAGS  =-g -std=c++11 -Wall -pedantic
LDFLAGS =-lncurses

OBJ = main.o DNA.o fastq.o 

fqless: $(OBJ)
	$(CC) $(CFLAGS) -o fqless $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<


