//---------------------------------------------------------------------------
#include "util98.h"
#include "mt19937.h"

//---------------------------------------------------------------------------

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define NR_END 1
#define FREE_ARG char*
#define ENDRUN '!'
#define COMMENT '*'
#define MAXOUTFILES 10
#define ALLFILES -1
#define MAKESTR(A,B,C) A ## B ## C

/* singly linked list of parameter structures */
/* used by the multiple run functions */
struct varlist plist = {"parmlist",NULL,NOTYPE,NULL};


void error(char error_text[])
/* Drew Tyre's standard error handler */
{
	fprintf(stderr,"Drew Tyre run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}  /* end function error */

int getnextrun(char rerunname[])
/* reads all meaningful lines in a rerun file, resets the named variables,
   and returns non-zero if another rerun is to be performed. It handles the
   opening of the rerun file itself. If the file does not exist, or if the
   end of the file has been reached, then the function returns zero */
{
	static long int oldfilepos = 0L;
    char buffer[LINELEN] = "\0", message[LINELEN] = "\0";
    unsigned int foundit = FALSE, check1 = 0, check2 = 0;
	 struct varlist *curparm;
    FILE *rerunfile;

    if ((rerunfile = fopen(rerunname,"r")) == NULL)
    {   /* no rerun file exists, so return 0 to calling program */
    	return 0;
    }

    /* move to last position read */
	 if (fseek(rerunfile,oldfilepos,SEEK_SET))
	{
		error("Bad oldfilepos in getnewrun");
    }

    /* 	read in file one line at a time, discarding lines beginning with
    	COMMENT, and stopping when a line starts with ENDRUN or EOF
		is reached */
    fgets(buffer,LINELEN,rerunfile);
	 while (!feof(rerunfile) && (buffer[0] != ENDRUN))
    {
        if (buffer[0] != COMMENT)
        {	/* line contains a variable, so find the variable in plist */
        	curparm = plist.next;
            foundit = FALSE;
            while (curparm != NULL)
            {
            	check1 = 0;
				check2 = 0;
            	while (buffer[check1++] == curparm->name[check2++]);
				if (check1 > strlen(curparm->name))
                {   /* found the variable in plist */
                	foundit = TRUE;
                	switch(curparm->type)
                    {
                    case REAL :
                    	sscanf(&buffer[check1]," %f",(float *) curparm->address);
								break;
                    case DOUBLE :
                    	sscanf(&buffer[check1]," %lf",(double *) curparm->address);
                        break;
                    case INTEGER :
                    	sscanf(&buffer[check1]," %d",(int *) curparm->address);
                        break;
                    case LONGINT :
                    	sscanf(&buffer[check1]," %ld",(long int *) curparm->address);
								break;
                    case STRING :
                    	sscanf(&buffer[check1]," %s",(char *) curparm->address);
                      break;
                    default :
                    	strcat(message,"Bad parmtype in plist, ");
                        strcat(message,curparm->name);
                        error(message);
                    }
						  /* exit loop over plist */
                    break;
                }
                /* skip to next entry in plist */
                curparm = curparm->next;
            }	/* end of loop over plist */

            if (!foundit)
            {
					strcat(message,"Bad variable name ");
				strcat(message,buffer);
				strcat(message," in rerunfile");
            	error(message);
            }
        }	/* end if buffer[0] != comment */
        /* get next line in rerunfile */
        fgets(buffer,LINELEN,rerunfile);

	 }	/* end of while !feof && buffer[0] != comment */

    /* remember where we reached in the file */
    oldfilepos = ftell(rerunfile);

    /* close the rerunfile, or the file handles overflow! */
    fclose(rerunfile);

    if (foundit)
	 {   /* found something, so do another run */
    	return 1;
    }
    else
    {	/* found nothing, so don't do another run */
		return 0;
	 }

}	/* end of getnewrun */

void rdsvar(FILE *infile, char varname[], void *address, int parmtype)
/* reads a variable of parmtype in from infile. Assumes that the next
	line is the one that has the variable in question but will wrap once to
   find the variable */
{
	int foundit, namelen, check1, check2, wrapped, gotit;
    char buffer[LINELEN] = "\0";	/* for storing the line found in the file */
    char message[LINELEN] = "\0";	/* for storing errormessages */

    namelen = strlen(varname);	/* figure out how long the string is */
    foundit = FALSE;
    wrapped = FALSE;

 	do
    {   /* loop through file looking for varname */
    	check1 = 0;
        check2 = 0;
        fgets(buffer,LINELEN,infile);
		if (feof(infile) && !(wrapped))
        {   /* reached the end of the file for the first time */
            wrapped = TRUE;
            rewind(infile);
        }

        while (buffer[check1++] == varname[check2++]);
        if (check1 > (namelen-1))
        {	/* varname matches upto namelen */
        	foundit = TRUE;

            switch(parmtype)
            {
            	case REAL :
					gotit = sscanf(&buffer[check1]," %f", (float *) address);
                    break;
            	case DOUBLE :
					gotit = sscanf(&buffer[check1]," %lf", (double *) address);
                    break;
                case INTEGER :
                	gotit = sscanf(&buffer[check1]," %d", (int *) address);
                    break;
                case LONGINT :
                	gotit = sscanf(&buffer[check1]," %ld", (long int *) address);
                    break;
                case STRING :
                	gotit = sscanf(&buffer[check1]," %s", (char *) address);
                  break;
                default :
                	strcat(message,"Bad parmtype in rdsvar, ");
                    strcat(message,varname);
                	error(message);
            }	/* end of switch(parmtype) */

            if (!gotit)
            {
            	strcat(message,"Unable to scan ");
                strcat(message,varname);
                error(message);
            }	/* end of if (!gotit) */

            if (!addsvar(&plist,varname,address,parmtype))
            {   /* unable to stick variable info into plist */
            	strcat(message,"Unable to put ");
                strcat(message,varname);
                strcat(message," in plist");
                error(message);
            }	/* end of if !addsvar */
        }	/* end of if (check1 > (namelen - 1)) */
    }	/* end of do-while */
    while (!foundit && !(feof(infile) && wrapped));
    if (!foundit)
    {
    	strcat(message, "Unable to find ");
        strcat(message,varname);
        error(message);
    }
    return;
}

int addsvar(struct varlist *list, char varname[],void *address,int parmtype)
/* adds a variable to a singly linked list plist. The list is used by the
   functions rdsets and rdrun to update variables for multiple runs */
{

    /* find the end of the list */
	while (list->next != NULL) list = list->next;

    /* allocate a new structure */
    list->next = (struct varlist *) malloc(sizeof(struct varlist));

    /* check for success */
    if (list->next == NULL) return 0;

    /* set list to new structure */
    list = list->next;

    /* fill in the blanks */
    strcpy(list->name,varname);
    list->address = address;
    list->type = parmtype;
    list->next = NULL;		/* VERY IMPORTANT!!! */

    return 1;	/* successfully added variable */

}

int dumpvars(struct varlist *list)
/* dumps the entire variable list */
{
	struct varlist *top, *next;

    /* can't deallocate the head of the list, so move to the next item */
    top = list->next;

    while (top != NULL)
    {
    	next = top->next;
        free(top);
        top = next;
    }

    return 0;	/* successfully dumped list */

}

double lint(char tname[], struct rtable table[], double  key)
/* This function performs a linear interpolation on a function stored in */
/* table. key gives the independent variable that is to be used to */
/* initiate the table lookup. lint assumes that key is sorted in */
/* order in table */
{
	double result, slope;
   int bottom, top;

   top = tsearch(tname, table, key);
   bottom = top - 1;

   slope = (table[top].yval - table[bottom].yval) /
   		  (table[top].xval - table[bottom].xval);

   result = (slope * (key - table[bottom].xval)) + table[bottom].yval;
   return result;
}

int tsearch(char tname[], struct rtable table[], double key)
/* This function searches table for the first value greater than or */
/* equal to key, and returns the array index of that value */
/* The function relies on the fact that the last key in the table */
/* will be much larger than any key that will be passed. The function also */
/* "remembers" which table it last searched, and where it reached in that */
/* table. The static variables lasttable and lastreached are used for this. */
/* This fact is used to speed up sequential searches on sorted tables. */
{
	int current, result;
   static char lasttable[VNAMELEN];
   static int lastreached;

   if((strcmp(tname,lasttable)==0) && (table[lastreached].xval < key))
	{    /* searching the same table as last time */
   	current = lastreached;
   }
   else
   {   /* searching a new table, or a previous value */
   	current = 0;
      strcpy(lasttable,tname);
   }

   /* loop through the table until the value is greater than or */
   /* equal the key */
   while (table[current++].xval < key);

	lastreached = current - 1; /* set the global variable lastreached */
   result = current - 1;     /* return a 'pointer' to the previous entry */
   return result;
}

void rdtable(FILE *initfile, char searchstr[], struct rtable table[])
/* Reads the x and y values from a table. Table has to have the */
/* following format : */
/* tablename tablelength */
/* xvalue yvalue         */
/*  ...    ...           */
{
	int i, arraysize;
	char foundstr[80];

	rdsvar(initfile,searchstr,&arraysize,INTEGER);

	if (arraysize > TABMAX)
	{
   		fprintf(stderr, "table '%s' has too many lines\n", searchstr);
	    exit(1);
	}

	for (i = 0; i <= (arraysize-1); i++)
	{
		fgets(foundstr,80,initfile);
	    sscanf(foundstr," %lf %lf",&table[i].xval,&table[i].yval);
	}

	return;
}

void rdvector(FILE *initfile, char searchstr[], void *vector, int parmtype)
/* Reads vector values from a list. List has to have the */
/* following format : */
/* vectorname # of entries */
/* value         */
/*  ...          */
/*********************************************************************/
/* NB!! assumes that sufficient space has been set aside in vector!! */
/*********************************************************************/
{
	int i, vectorsize;
	char foundstr[80];
  char message[LINELEN] = "\0";	/* for storing errormessages */

  /* use rdsvar to find the start of the list */
	rdsvar(initfile,searchstr,&vectorsize,INTEGER);

	for (i = 0; i <= (vectorsize-1); i++)
	{
		fgets(foundstr,80,initfile);
    switch(parmtype){
	    case REAL :
		    sscanf(foundstr," %f", &(((float *) vector)[i]));
        break;
      case INTEGER :
      	sscanf(foundstr," %d", &(((int *) vector)[i]));
         break;
      case DOUBLE :
      	sscanf(foundstr," %lf", &(((double *) vector)[i]));
        break;
      default :
      	strcat(message,"Bad parmtype in rdvector , ");
        strcat(message,searchstr);
        error(message);
      }	/* end switch */
	}

	return;
}

double chop(double x, int a, int b)
{
	if (x < a)
	{
		x = a;
	}
	else if (x > b)
	{
		x = b;
	}

	return x;
}   /* end of chop */



/*/*int round(double x)
{
	 int dum1, result;
	 double dum2;

	 dum1 = x;
	 dum2 = x - dum1;
	 if(dum2 < 0.5)
	 {
		result =  floor(x);   /* round x down /
	 }
	 else
	 {
		result = ceil(x);     /* round x up /
		}
	return result;
} /* end of round */

/* Returns the probability of finding an odourconcentration of */
/*	i individuals in a quadrat, */
/* given a mean of 'mean' and the aggregation coefficient 'k' */
/* Taken from Krebs 1989 Ecological Methodology, page 82 */
/* I assume that n hosts produce an odourconc of n */

double negbin(int i, double k, double mean)
{
	 double tmp1, tmp2, tmp3, result;

	 tmp1 = exp(gammln(k + i)) / (factrl(i) * exp(gammln(k)));
	 tmp2 = mean / (mean + k);
	 tmp3 = k / (k + mean);
	 if (tmp2 < 0) printf("Domain error: tmp2\n");
	 if (tmp3 < 0) printf("Domain error: tmp3\n");
	 result = tmp1 * pow(tmp2,(double)i) * pow(tmp3,k);
	 return result;
}

/* Computes the log of Lanczos approximation to gamma function */
/* From Numerical recepies in C. error < 2E-10!! */
/* Some modifications were also taken from the numerical recipes book */
/* in fortran  */
double gammln(double xx)
{
	 double czero = 1.000000000190015;
	 double roottwopi = 2.5066282746310005;
	 double cof[8] = {0, 76.1800917294716, -86.50532032941677,
							24.01409824083091, -1.231739572450155,
							0.1208650973866179E-2, -0.5395239384953E-5};
	 double x, series, tmp1, tmp2;
	 int j;
	 double result;

	 x = xx - 1;
	 tmp1 = x + 5.5;
	 tmp2 = (x + 0.5) * log(tmp1) - tmp1;

	 /* initialize the series with the first term */
	 series = czero;

	 /* loop summs each coefficient divided by x+j */
	 for(j = 1;j <= 6;j ++)
	 {
		series += cof[j] / (x + j);
	 }
	 result = tmp2 + log(roottwopi * series);    /* return this value */

	 return result;
}


/* Returns the factorial of n as a doubleing point number */
double factrl(int n)
{
	double result;
	int ntop = 4;    /* This is the top of the table */
	/* This is the table of results, fill up as required */
	double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};

	if(n < 0)
	{
		printf("Negative factorial in routine factl\n");
		exit(n);
	}

	if(n > 32)    /* larger value than is in table, probably will overflow */
	{
		result = exp(gammln(n + 1.0));
	}
	else
	{
		while(ntop < n)         /* fill in table up to required value */
		{
			ntop ++;
			a[ntop] = a[ntop - 1] * ntop;
		}
		result = a[n];        /* return table value */
	}

	return result;
}

/* This function calculates the probability that there are i hosts on
	 a plant, given the plant is occupied */
double con_negbin(int i, double k, double mean)
{
	double result, b, a;

	a = negbin(i,k,mean);
	b = 1-negbin(0,k,mean);
	result = a/b;

	return result;
}
/* Returns a normally distributed deviate with zero mean and unit variance,
	using ran1(idum) as the source of uniform deviates. Taken from
	numerical recipes book in c, p.289 */

float gasdev(void)
{
/*	float ran3(long *idum); */
	static int iset=0;
	static float gset;
	float fac, rsq, v1, v2, temp;

	/* We don't have an extra deviate handy, so pick 2 uniform numbers in
		the square extending from -1 to +1 in each direction, see if they
		are in the unit circle, and if they are not, try again.*/
	if (iset == 0){
		do{
			v1 = 2.0*rand0to1()-1.0;
			v2 = 2.0*rand0to1()-1.0;
			rsq=v1*v1+v2*v2;
		}while (rsq >= 1.0 ||rsq == 0.0);

      temp = -2.0*log(rsq)/rsq;
      if (temp >= 0.0) fac = sqrt(temp);
      else {
      	printf("error in gasdev\n");
		}
		/* Now make the Box-Muller transformation to get 2 normal deviates.
		Return one and save the other for the next time. */
		gset = v1*fac;
		iset=1; 				/*set flag */
		return v2*fac;
	} else{
		/* We have an extra deviate handy, so unset the flag, and return it.*/
		iset = 0;
		return gset;
	}
}

int randompick(double prob[], int maxint)
{
	double x, sumprob = 0;
	int z = 0, choose = 0, i = -1;

	/* x determines the target variable */
	x = genrand();

	do
	{
		i++;
		sumprob += prob[i];
		if (x <= sumprob)
		{
			choose = i;
			z = 1;
		}
	} while (!z);

	if (choose > maxint){
		fprintf(stderr,"choose: %d maxint: %d\n",choose,maxint);
		for (i = 0; i<=choose; i++) fprintf(stderr,"%d %f\n",i,prob[i]);
		fprintf(stderr,"bad probability sum in randompick!\n");
		exit(1);
	}

	return choose;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) error("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void hqr(float **a, int n, float wr[], float wi[])
{
	int nn,m,l,k,j,its,i,mmin;
	float z,y,x,w,v,u,t,s,r,q,p,anorm;

	anorm=fabs(a[1][1]);
	for (i=2;i<=n;i++)
		for (j=(i-1);j<=n;j++)
			anorm += fabs(a[i][j]);
	nn=n;
	t=0.0;
	while (nn >= 1) {
		its=0;
		do {
			for (l=nn;l>=2;l--) {
				s=fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s == 0.0) s=anorm;
				if ((float)(fabs(a[l][l-1]) + s) == s) break;
			}
			x=a[nn][nn];
			if (l == nn) {
				wr[nn]=x+t;
				wi[nn--]=0.0;
			} else {
				y=a[nn-1][nn-1];
				w=a[nn][nn-1]*a[nn-1][nn];
				if (l == (nn-1)) {
					p=0.5*(y-x);
					q=p*p+w;
					z=sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z=p+SIGN(z,p);
						wr[nn-1]=wr[nn]=x+z;
						if (z) wr[nn]=x-w/z;
						wi[nn-1]=wi[nn]=0.0;
					} else {
						wr[nn-1]=wr[nn]=x+p;
						wi[nn-1]= -(wi[nn]=z);
					}
					nn -= 2;
				} else {
					if (its == 30) error("Too many iterations in hqr");
					if (its == 10 || its == 20) {
						t += x;
						for (i=1;i<=nn;i++) a[i][i] -= x;
						s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y=x=0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2);m>=l;m--) {
						z=a[m][m];
						r=x-z;
						s=y-z;
						p=(r*s-w)/a[m+1][m]+a[m][m+1];
						q=a[m+1][m+1]-z-r-s;
						r=a[m+2][m+1];
						s=fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if ((float)(u+v) == v) break;
					}
					for (i=m+2;i<=nn;i++) {
						a[i][i-2]=0.0;
						if (i != (m+2)) a[i][i-3]=0.0;
					}
					for (k=m;k<=nn-1;k++) {
						if (k != m) {
							p=a[k][k-1];
							q=a[k+1][k-1];
							r=0.0;
							if (k != (nn-1)) r=a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							if (k == m) {
								if (l != m)
								a[k][k-1] = -a[k][k-1];
							} else
								a[k][k-1] = -s*x;
							p += s;
							x=p/s;
							y=q/s;
							z=r/s;
							q /= p;
							r /= p;
							for (j=k;j<=nn;j++) {
								p=a[k][j]+q*a[k+1][j];
								if (k != (nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn<k+3 ? nn : k+3;
							for (i=l;i<=mmin;i++) {
								p=x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l < nn-1);
	}
}

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

float dfridr(float (*func)(float), float x, float h, float *err)
{
	int i,j;
	float errt,fac,hh,**a,ans;

	if (h == 0.0) error("h must be nonzero in dfridr.");
	a=matrix(1,NTAB,1,NTAB);
	hh=h;
	a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) {
			free_matrix(a,1,NTAB,1,NTAB);
			return ans;
		}
	}
	free_matrix(a,1,NTAB,1,NTAB);
	return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE

/* following 3 routines do romberg integration on closed intervals */
void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) error("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}

#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

float qromb(float (*func)(float), float a, float b)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	float trapzd(float (*func)(float), float a, float b, int n);
	float ss,dss;
	float s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	error("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

