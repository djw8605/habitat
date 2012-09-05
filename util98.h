//---------------------------------------------------------------------------

#ifndef util98H
#define util98H
//---------------------------------------------------------------------------
#include <stdio.h>

/* nr macros */
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void error(char error_text[]);
/* parameter list structure */
#define NOTYPE 0
#define INTEGER 1
#define LONGINT 2
#define REAL 3
#define DOUBLE 4
#define STRING 5
#define VNAMELEN 15
#define VPRECISION 3


struct varlist
{
	char name[VNAMELEN];
	 void *address;
	 int type;
	 struct varlist *next;
};

/* table definitions */
#define TABMAX 50
struct rtable
{
	double xval;
	 double yval;
};

/* utility definitions */
#define FALSE 0
#define TRUE 1
#define LINELEN 80

/* parameter reading functions */
extern "C" void rdsvar(FILE *infile, char varname[], void *varvalue, int parmtype);

/* multiple run reading functions */
extern "C" int getnextrun(char rerunname[]);

/* variable list maintenance functions */
extern "C" int addsvar(struct varlist *list, char varname[], void *address, int parmtype);
extern "C" int dumpvars(struct varlist *list);

/* table interpolation functions */
extern "C" double lint(char tname[], struct rtable table[], double key);
extern "C" int tsearch(char tname[], struct rtable table[], double key);
extern "C" void rdtable(FILE *initfiile, char searchstr[], struct rtable table[]);

/* vector reading functions */
extern "C" void rdvector(FILE *initfile, char searchstr[], void *vector, int parmtype);

/* maths functions for rubecula model */
extern "C" double chop(double x, int a, int b);
/*extern "C" int round(double x);*/
extern "C" double con_negbin(int i, double k, double mean);
extern "C" double negbin(int i, double k, double mean);
extern "C" double gammln(double xx);
extern "C" double factrl(int n);
extern "C" float gasdev(void);
extern "C" int randompick(double prob[], int maxint);

/* nr functions */
extern "C" float **matrix(long nrl, long nrh, long ncl, long nch);
extern "C" float *vector(long nl, long nh);
extern "C" void free_vector(float *v, long nl, long nh);
extern "C" void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
extern "C" void hqr(float **a, int n, float wr[], float wi[]);
extern "C" float dfridr(float (*func)(float), float x, float h, float *err);
extern "C" float qromb(float (*func)(float), float a, float b);

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


/* randomize and random are not defined in ANSI C */
#define randomize() srand((unsigned)time(NULL))
#define rand0to1() genrand()
#define random(a) (igenrand() % a)

#ifndef __STDC__
#include <alloc.h>
extern int checkmemory(void);
#endif

#endif
