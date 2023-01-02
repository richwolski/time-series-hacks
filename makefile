CC = gcc
CFLAGS = "-g"

all: diffseries winmax interarrival-run markov-series markov-seasonal-series

diffseries: diffseries.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o diffseries diffseries.c simple_input.o 

winmax: winmax.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o winmax winmax.c simple_input.o 

interarrival-run: interarrival-run.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o interarrival-run interarrival-run.c simple_input.o 

markov-series: markov-series.c simple_input.o simple_input.h red_black.h Hval.h red_black.o
	${CC} ${CFLAGS} -DTEST -o markov-series markov-series.c simple_input.o red_black.o

markov-seasonal-series: markov-seasonal-series.c simple_input.o simple_input.h red_black.h Hval.h red_black.o
	${CC} ${CFLAGS} -DTEST -o markov-seasonal-series markov-seasonal-series.c simple_input.o red_black.o

simple_input.o: simple_input.h simple_input.c
	${CC} ${CFLAGS} -c simple_input.c

red_black.o: red_black.h red_black.c Hval.h
	${CC} ${CFLAGS} -c red_black.c
