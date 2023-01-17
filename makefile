CC = gcc
CFLAGS = "-g"

ALL = diffseries winmax interarrival-run markov-series markov-seasonal-series inar-seasonal-series
all: ${ALL}

MIO=../mio
MIOLIB=${MIO}/mio.o ${MIO}/mymalloc.o
DIST=../distributions
DISTLIB=${DIST}/poisson.o

diffseries: diffseries.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o diffseries diffseries.c simple_input.o 

winmax: winmax.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o winmax winmax.c simple_input.o 

interarrival-run: interarrival-run.c simple_input.o simple_input.h
	${CC} ${CFLAGS} -DTEST -o interarrival-run interarrival-run.c simple_input.o 

markov-series: markov-series.c simple_input.o simple_input.h redblack.h hval.h redblack.o
	${CC} ${CFLAGS} -DTEST -o markov-series markov-series.c simple_input.o redblack.o

markov-seasonal-series: markov-seasonal-series.c simple_input.o simple_input.h redblack.h hval.h redblack.o
	${CC} ${CFLAGS} -DTEST -o markov-seasonal-series markov-seasonal-series.c simple_input.o redblack.o

inar-seasonal-series: inar-seasonal-series.c ${MIO}/mio.h ${MIOLIB} ${DIST}/poisson.h ${DISTLIB}
	${CC} ${CFLAGS} -I${MIO} -I${DIST} -o inar-seasonal-series inar-seasonal-series.c ${MIOLIB} ${DISTLIB}

simple_input.o: simple_input.h simple_input.c
	${CC} ${CFLAGS} -c simple_input.c

redblack.o: redblack.h redblack.c hval.h
	${CC} ${CFLAGS} -c redblack.c

clean:
	rm -f *.o ${ALL}
