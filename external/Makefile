UNAME_S := $(shell uname -s)

CC:=gcc
CFLAGS:=-g -Wall -O2 -Wno-unused-function
BINDIR:=/usr/local/bin

all: seqtk

seqtk: seqtk.c rand_win.c khash.h kseq.h rand_win.h
	@if echo $(UNAME_S) | grep -q MINGW; then \
		echo "Patching seqtk.c to include rand_win.h..."; \
		sed -i '1i#ifdef _WIN32\n#include "rand_win.h"\n#endif' seqtk.c; \
	fi
	$(CC) $(CFLAGS) seqtk.c rand_win.c -o $@ -lz -lm

install: all
	install seqtk $(BINDIR)

clean:
	rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*