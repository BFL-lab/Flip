#
# Makefile for flip
#
# written by Alexander Nip

CC=gcc
CFLAGS=-O3 
LIBS=$(SYSLIB) -lm
OBJS=flip.o

# Rules to create .o files from .c files
.c.o:
	$(CC) -g -c $(CFLAGS) $(INCLUDES) $*.c

all: flip

flip:	$(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(LIBS) -o $@

clean:
	-rm -f $(OBJS)
	-rm -f core
	-rm -f flip

flip.o:	flip.c code.h flip.h
