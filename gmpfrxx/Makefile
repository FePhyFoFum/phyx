cc     = gcc
CC     = g++
CFLAGS = -g -Wall

INC    =
LIBS   = -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm

#ROOT   = /Users/wilken
#INC    = -I$(ROOT)/include
#LIBS   = -L. -L$(ROOT)/lib -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm

headers = $(shell ls *.h)

all : libgmpfrxx.a

absolute : libgmpfrxx.a example

libgmpfrxx.a : gmpfrxx.o mpfr_mul_d.o
	ar rs libgmpfrxx.a gmpfrxx.o mpfr_mul_d.o

example : example.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

%.o : %.f
	$(FC) -c $(FFLAGS) $<

%.o : %.cc $(headers)
	$(CC) -c $(CFLAGS) $(INC) $<

%.o : %.cpp $(headers)
	$(CC) -c $(CFLAGS) $(INC) $<

%.o : %.c $(headers)
	$(cc) -c $(CFLAGS) $(INC) $<

files = COPYING README example.cpp gmpfrxx.h gmpfrxx.cpp \
	mpfr_mul_d.h mpfr_mul_d.c Makefile

.PHONY : always
gmpfrxx.zip : always
	cp -p $(files) gmpfrxx
	zip -r gmpfrxx gmpfrxx

.PHONY : clean
clean:
	rm -f *.o example
