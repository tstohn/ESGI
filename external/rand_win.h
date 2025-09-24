// rand_win.h
#ifndef _LRAND48_WIN_H
#define _LRAND48_WIN_H

#ifdef _WIN32
void srand48(long seedval);
long lrand48(void);
double drand48(void);
#endif

#endif