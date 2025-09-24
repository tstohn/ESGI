// lrand48_win.c – shim for Windows to support lrand48/drand48/srand48

#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "rand_win.h"

// internal state
static uint64_t lrand48_state = 0;

void srand48(long seedval) 
{
    lrand48_state = ((uint64_t)seedval << 16) | 0x330E;
}

long lrand48(void) 
{
    lrand48_state = (0x5DEECE66DULL * lrand48_state + 0xB) & ((1ULL << 48) - 1);
    return (long)(lrand48_state >> 17);
}

double drand48(void) 
{
    return lrand48() / (double)(1L << 31);
}
