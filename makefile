CC = gcc
CFLAGS = "-g"

ALL = diffseries winmax interarrival-run markov-series markov-seasonal-series inar-seasonal-series inar-series acp acp-gen acp-seasonal
all: ${ALL}

MIO=../mio
APATH=../matrix
MIOLIB=${MIO}/mio.o ${MIO}/mymalloc.o
DIST=../distributions
DISTLIB=${DIST}/poisson.o -lm
#OSX
#YWLIB=../meanvar/meanvar.o ../autoc/autoc.o ../autoc/yw-estimate.o -L/opt/homebrew/opt/lapack/lib -L /opt/homebrew/opt/openblas/lib ${APATH}/mioarray.o -llapacke -llapack -lopenblas
YWLIB=../meanvar/meanvar.o ../autoc/autoc.o ../autoc/yw-estimate.o ${APATH}/mioarray.o ${APATH}/lapack-3.8.0/liblapacke.a ${APATH}/lapack-3.8.0/liblapack.a -lgfortran -L/usr/lib64 -L/usr/lib64/atlas -lsatlas -lblas -lm

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

inar-series: inar-series.c ../autoc/yw-estimate.h ${MIO}/mio.h ${MIOLIB} ${DIST}/poisson.h ${DISTLIB} 
	${CC} ${CFLAGS} -I../autoc -I../meanvar -I${MIO} -I${DIST} -o inar-series inar-series.c ${MIOLIB} ${DISTLIB} ${YWLIB}

acp: acp.c ${MIO}/mio.h ${MIOLIB} ${DIST}/poisson.h ${DISTLIB} 
	${CC} ${CFLAGS} -I../meanvar -I${MIO} -I${DIST} -o acp acp.c ${MIOLIB} ${DISTLIB} ${YWLIB}

acp-seasonal: acp-seasonal.c ${MIO}/mio.h ${MIOLIB} ${DIST}/poisson.h ${DISTLIB} 
	${CC} ${CFLAGS} -I../meanvar -I${MIO} -I${DIST} -o acp-seasonal acp-seasonal.c ${MIOLIB} ${DISTLIB} ${YWLIB}

acp-gen: acp-gen.c ${DIST}/poisson.h ${DISTLIB} 
	${CC} ${CFLAGS} -I${DIST} -o acp-gen acp-gen.c ${DISTLIB}

simple_input.o: simple_input.h simple_input.c
	${CC} ${CFLAGS} -c simple_input.c

redblack.o: redblack.h redblack.c hval.h
	${CC} ${CFLAGS} -c redblack.c

clean:
	rm -f *.o ${ALL}
