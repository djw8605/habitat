//---------------------------------------------------------------------------

#ifndef mt19937H
#define mt19937H
//---------------------------------------------------------------------------
void sgenrand(unsigned int seed); /* initialises the generator */
double genrand(void);   /* returns a single double precision real on [0..1] */
unsigned long igenrand(void);   /* returns a single long int on [0..MAXINT] */
#endif
