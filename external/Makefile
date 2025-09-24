UNAME_S := $(shell uname -s)

CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function
BINDIR=/usr/local/bin

WIN_SHIM = rand_win.c
#ifeq ($(findstring MINGW, $(UNAME_S)),MINGW)
#	WIN_SHIM = rand_win.c
#endif

#ifeq ($(findstring MSYS, $(UNAME_S)),MSYS)
#	WIN_SHIM = rand_win.c
#endif

all: seqtk

seqtk: seqtk.c rand_win.c khash.h kseq.h rand_win.h
	$(CC) $(CFLAGS) seqtk.c rand_win.c -o $@ -lz -lm

install: all
	install seqtk $(BINDIR)

clean:
	rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*