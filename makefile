CC = gcc
CFLAGS = "-g"

all: diffseries winmax

diffseries: diffseries.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o diffseries diffseries.c simple_input.o 

winmax: winmax.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o winmax winmax.c simple_input.o 

simple_input.o: simple_input.h simple_input.c
	${CC} ${CFLAGS} -c simple_input.c

