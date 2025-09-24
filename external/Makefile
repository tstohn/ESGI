UNAME_S := $(shell uname -s)

CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
BINDIR=/usr/local/bin

WIN_SHIM = rand_win.c
ifeq ($(findstring MINGW, $(UNAME_S)),MINGW)
	WIN_SHIM = rand_win.c
endif

ifeq ($(findstring MSYS, $(UNAME_S)),MSYS)
	WIN_SHIM = rand_win.c
endif

all: seqtk

seqtk: seqtk.c $(WIN_SHIM) khash.h kseq.h
	$(CC) $(CFLAGS) seqtk.c $(WIN_SHIM) -o $@ -lz -lm

install: all
	install seqtk $(BINDIR)

clean:
	rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*