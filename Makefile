CC = gcc
CFLAGS = -W -Wall -std=gnu99 -ffast-math -fgnu89-inline -O3
LDFLAGS = -lm -lgsl -lgslcblas 
OBJ = main.o io.o aux.o train.o

hmm: $(OBJ)
		$(CC) $(OBJ) $(LDFLAGS) -o hmm

aux.o: aux.c aux.h
io.o: io.c io.h
train.o: train.c train.h datastruct.h
main.o: main.c io.h datastruct.h aux.h train.h

clean: 
		rm -f $(OBJ) hmm













