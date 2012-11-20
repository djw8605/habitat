//---------------------------------------------------------------------------
#include "landscpe.h"

//---------------------------------------------------------------------------
/*************************************************************************
This file implements a landscape with variable geometry
Last Modified: 18/03/97
Written by: Drew Tyre, Dept. of Environmental Science and Management,
University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
Ph: (08) 303 7931 Fax: (08) 303 7956
email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "util98.h"
#include <math.h>
#include <limits.h>
#include <string.h>
#include "habitat.h"


int MakeFractalLandscape(LandScape *L, LParameters *LPm);
int MidPointFM2D(float **X, int maxlevel, float sigma, float H,
    int addition);

int quantile(int *H, int N, int q);
int hqCompare(const void *e1,  const void *e2);


#define MAXBOUNDARYTYPE 12
#define MAXDIRECTION FULLSQR
#define CENTER  0
#define EVENROW 0
#define TOP 1
#define BOTTOM 2
#define LEFT 4
#define RIGHT 8
#define ODD -1
#define ODDROW 11
#define LEFTODD (LEFT + ODD)
#define RIGHTODD (RIGHT + ODD)
#define TOPLEFT (TOP + LEFT)
#define TOPRIGHT (TOP + RIGHT)
#define BOTLEFT (BOTTOM + LEFT)
#define BOTRIGHT (BOTTOM + RIGHT)

/* next two constants used to calculate the distance between two territories
on a hex grid */
#define HEXROW (0.866025) /* distance between two rows in a hex grid */
#define HEXCOL (-0.5)         /* distance by which odd rows are offset horizontally */
int OffSetArray[MAXBOUNDARYTYPE][MAXDIRECTION] = {{0}};

int FullSqrTemplate[][MAXDIRECTION] =
{{0,0,0,0,0,0,0,0}, /* center is empty */
{0,0,1,1,1,1,1,0},  /* top edge */
{1,1,1,0,0,0,1,1},  /* bottom edge */
{0,0,0,0,0,0,0,0},   /* spacer row */
{1,1,1,1,1,0,0,0},  /* left edge */
{0,0,1,1,1,0,0,0},   /* Top Left corner */
{1,1,1,0,0,0,0,0},   /* Bottom Left */
{0,0,0,0,0,0,0,0},   /* spacer row */
{1,0,0,0,1,1,1,1},   /* right edge */
{0,0,0,0,1,1,1,0},   /* top right corner */
{1,0,0,0,0,0,1,1}};  /* bottom right corner */

int OrthoTemplate[][MAXDIRECTION] =
{{0,0,0,0}, /* center is empty */
{0,1,1,1},  /* top edge */
{1,1,0,1},   /* bottom edge */
{0,0,0,0},   /* spacer row */
{1,1,1,0},   /* left edge */
{0,1,1,0},   /* top left corner */
{1,1,0,0},   /* bottom left corner */
{0,0,0,0},  /* spacer row */
{1,0,1,1},   /* right edge */
{0,0,1,1},  /* top right corner */
{1,0,0,1}};  /* bottom right corner */

int HexTemplate[][MAXDIRECTION] =
{{0,0,0,0,0,0}, /* center is empty (even rows) */
{0,1,1,1,1,0},  /* top edge */
{1,1,0,0,1,1}, /* bottom edge */
{1,1,1,0,0,0}, /* left edge, odd row */
{1,1,1,1,0,1},   /* left edge, even row */
{0,1,1,1,0,0},  /* top left corner */
{1,1,0,0,0,0},   /* bottom left corner */
{1,0,1,1,1,1},   /* right edge, odd row */
{0,0,0,1,1,1}, /* right edge, even row */
{0,0,0,1,1,0}, /* top right corner */
{1,0,0,0,1,1}, /* bottom left corner */
{0,0,0,0,0,0}};  /* odd rows in center */

void InitOffSet(Parameters *Pm)
{
    int i,j;
    switch (Pm->Sides)
    {
    case FULLSQR :
        OffSetArray[CENTER][0] = -Pm->Xsize;        /* Straight up */
        OffSetArray[CENTER][1] = -Pm->Xsize + 1;   /* up and right */
        OffSetArray[CENTER][2] = 1;                 /* right */
        OffSetArray[CENTER][3] = Pm->Xsize + 1;    /* down and right */
        OffSetArray[CENTER][4] = Pm->Xsize;         /* down */
        OffSetArray[CENTER][5] = Pm->Xsize - 1;    /* down and left */
        OffSetArray[CENTER][6] = -1;                   /* left */
        OffSetArray[CENTER][7] = -Pm->Xsize - 1;   /* up and left */
        for(i = 1; i < MAXBOUNDARYTYPE; i++)
            for(j = 0; j < Pm->Sides; j++)
                OffSetArray[i][j] = FullSqrTemplate[i][j]*OffSetArray[CENTER][j];
        break;
    case ORTHO :
        OffSetArray[CENTER][0] = -Pm->Xsize;        /* Straight up */
        OffSetArray[CENTER][1] = 1;                    /* right */
        OffSetArray[CENTER][2] = Pm->Xsize;            /* down */
        OffSetArray[CENTER][3] = -1;                   /* left */
        for(i = 1; i < MAXBOUNDARYTYPE; i++)
            for(j = 0; j < Pm->Sides; j++)
                OffSetArray[i][j] = OrthoTemplate[i][j]*OffSetArray[CENTER][j];
        break;
    case HEX :
        OffSetArray[EVENROW][0] = -Pm->Xsize + 1;   /* up and right */
        OffSetArray[EVENROW][1] = 1;                   /* right */
        OffSetArray[EVENROW][2] = Pm->Xsize + 1;   /* down and right */
        OffSetArray[EVENROW][3] = Pm->Xsize;           /* down and left */
        OffSetArray[EVENROW][4] = -1;                  /* left */
        OffSetArray[EVENROW][5] = -Pm->Xsize;       /* up and left */
        OffSetArray[ODDROW][0]  = -Pm->Xsize;      /* up and right */
        OffSetArray[ODDROW][1]  = 1;                    /* right */
        OffSetArray[ODDROW][2]  = Pm->Xsize;           /* down and right */
        OffSetArray[ODDROW][3]  = Pm->Xsize - 1;    /* down and left */
        OffSetArray[ODDROW][4]  = -1;                  /* left */
        OffSetArray[ODDROW][5]  = -Pm->Xsize - 1;  /* up and left */
        for(j = 0; j < Pm->Sides; j++){
            OffSetArray[TOP][j] = HexTemplate[TOP][j]*OffSetArray[EVENROW][j];
            OffSetArray[BOTTOM][j] = HexTemplate[BOTTOM][j]*OffSetArray[ODDROW][j];
            OffSetArray[LEFTODD][j] = HexTemplate[LEFTODD][j]*OffSetArray[ODDROW][j];
            OffSetArray[LEFT][j] = HexTemplate[LEFT][j]*OffSetArray[EVENROW][j];
            OffSetArray[TOPLEFT][j] = HexTemplate[TOPLEFT][j]*OffSetArray[EVENROW][j];
            OffSetArray[BOTLEFT][j] = HexTemplate[BOTLEFT][j]*OffSetArray[ODDROW][j];
            OffSetArray[RIGHTODD][j] = HexTemplate[RIGHTODD][j]*OffSetArray[ODDROW][j];
            OffSetArray[RIGHT][j] = HexTemplate[RIGHT][j]*OffSetArray[EVENROW][j];
            OffSetArray[TOPRIGHT][j] = HexTemplate[TOPRIGHT][j]*OffSetArray[EVENROW][j];
            OffSetArray[BOTRIGHT][j] = HexTemplate[BOTRIGHT][j]*OffSetArray[ODDROW][j];
        }
        break;
    default :
        error("Bad sides parameter in InitOffSet");
        break;
    }

    return;
}

int CheckForBoundary(LandScape *L, Location Where){
    /* returns an integer that uniquely identifies the edge position */
    int i,j,Boundary = 0;

    i = Where / L->Xsize;
    j = Where % L->Xsize;

    if (j == 0) Boundary += LEFT;                                /* on left edge */
    else if (j == (L->Xsize - 1)) Boundary += RIGHT;   /* on right edge */

    if (i == 0) Boundary += TOP;                             /* on top edge */
    else if (i == (L->Ysize-1)) Boundary += BOTTOM;        /* on bottom edge */
    else if (L->Sides == HEX){
        if ((i % 2) > 0) Boundary += ODD;            /* in odd row; hex grid only */
        if (Boundary < 0) Boundary = ODDROW;
    }

    return Boundary;
}


LandScape *InitLandScape(Parameters *Pm, LParameters *LPm){
    /* this function creates a square empty landscape */
    int i,j;
    LandScape *L = NULL;
    char Message[80];
    //    unsigned int temp;
    void *tempptr;
    HQuality Value;

    if ((L = (LandScape *) malloc((size_t)sizeof(LandScape))) == NULL){
        error("Insufficient Memory for Landscape");
    }

    /* initialize landscape parameters */
    L->NumLayers = Pm->NumLayers;
    L->Sides = Pm->Sides;
    L->Xsize = Pm->Xsize;
    L->Ysize = Pm->Ysize;
    L->Map = NULL;

    /* Init Locks for each cell */
    if((L->CellLocks = (omp_lock_t *) malloc(sizeof(omp_lock_t) * (L->Xsize * L->Ysize + 1))) == NULL)
    {
        error("Insufficient Memory for Landscape Locks");
    }
    for(i=0;i<=(Pm->Xsize*Pm->Ysize);i++)
    {
        omp_init_lock(&L->CellLocks[i]);
    }

    /* Initialize habitat layers */
    for (i=0; i<Pm->NumLayers; i++){
        if ((L->H[i] = (HQuality *)
            malloc((size_t)sizeof(HQuality)*(Pm->Xsize*Pm->Ysize + 1))) == NULL){
                error("Insufficient Memory for Habitat Layer");
        }
        /* set habitat values to MeanVal */
        ClearLand(L->H[i],LPm->meanval[i],Pm);
    }   /* end for i < NumLayers */

    /* create fractal landscape in first layer */
    if (LPm->Hdim > 0){
        (void) MakeFractalLandscape(L, LPm);
    }
    else{   /* fill landscape completely randomly */
        for (i=0;i<Pm->Xsize*Pm->Ysize;i++)
            L->H[0][i] = (HQuality) chop(gasdev()*17.1667 + 50,0,100);/* random numbers btw. 0 and 99 */
    }

    /* Allocate Map Layer */
    //    temp = Pm->Xsize*Pm->Ysize+1;
    if ((L->Map = (LandCell *) malloc((size_t)sizeof(LandCell)*(Pm->Xsize*Pm->Ysize+1))) == NULL){
        error("Insufficient Memory for Map Layer");
    }

    /* loop over map and set all pointers to zero */
    for(i=0;i<=(Pm->Xsize*Pm->Ysize);i++){
        L->Map[i] = NULL;
    }

    /* initialize offset array */
    InitOffSet(Pm);

    /* exit(0);   /* exit here for debugging new landscape functions */

    return L;
}   /* end function InitLandScape */

int DestroyLandscape(LandScape **L){
    /* deallocate landscape storage and set pointer to NULL */
    int i;
    for (i=0;i<(*L)->NumLayers;i++) free((*L)->H[i]);
    free((*L)->Map);
    free(*L);
    *L = NULL;
    return 0;
}

int GetDirection(LandScape *L){
    /* return a random direction from the landscape */
    return (int) floor(rand0to1()*L->Sides);
}

HQuality GetValue(LandScape *L, int Layer, Location Where){
    /* return the value of the habitat variable in Layer at Where */
    return L->H[Layer][(int)Where];
}

Location GetNextCell(LandScape *L, Location Where){
    /* return a pointer to a cell in a random direction from Where */
    Location New = 0;
    int a,b,c;
    do
    {
        a = CheckForBoundary(L,Where);
        b = GetDirection(L);
        c = OffSetArray[a][b];
        New = c;
    }while (New == 0);   /* keep looking for a direction that is valid */
    New += Where;
    if (New > L->Xsize*L->Ysize){
        printf("oops\n");
    }

    return New;
}   /* end function GetNextCell() */

int GetCellinDir(LandScape *L, Location Where, int direction){
    /* get address of cell at direction from Where */
    int a, b, c;
    int NewCell;

    a = CheckForBoundary(L,Where);   /* compute boundary position */
    c = OffSetArray[a][direction];               /* get offset value */
    if (c != 0)                             /* if valid direction */
        NewCell = Where + c;
    else                                     /* else can't go that way */
        NewCell = -1;

    return NewCell;
}

int GetNeighbours(LandScape *L, Location Where, HQuality *hq, int Force){
    /* fill array hq with the qualities of neighbouring cells
    invalid cells get quality -1, occupied cells are decremented by the
    value of Force (usually the number of steps taken to reach the current
    point assumes that hq has at least FULLSQR elements for safety */
    int i;
    int a, c;

    for (i=0; i<L->Sides; i++){
        a = CheckForBoundary(L,Where);  /* compute boundary position */
        c = OffSetArray[a][i];              /* get offset value */
        if (c != 0)                             /* if valid direction */
            hq[i] = L->H[0][Where+c];
        else                                        /* else can't go that way */
            hq[i] = -1;
        if (Force)                              /* if only want unoccupied cells */
            hq[i] = L->Map[Where+c] ? hq[i]-Force : hq[i];  /* if occupied set to 0 */
    }

    return L->Sides;
}

Location GetCell(LandScape *L, int force){
    /* return the location of a cell from the landscape at random */
    /* force == 0 causes only a single grab. force <> 0 looks for an
    unoccupied cell */
    /* could be problems if landscape has more than 64 K cells */
    Location New;
    float a,d;
    do{
        a = rand0to1();
        d = floor( L->Xsize*L->Ysize * a );
        New = (unsigned int) d;
    }while(L->Map[New] && force);   /* force == 0 causes only a single grab */
    return New;
}


void printmap(LandScape *L){
    int i,j,k;

    printf("____________________\n");
    for(j=0;j<L->Ysize;j++){
        for(i=0;i<L->Xsize;i++){
            k = j * L->Xsize + i;
            switch(((Individual*)L->Map[k])->AgeClass){
            case NewBorn :
                printf("N ");
                break;
            case Juvenile :
                printf("J ");
                break;
            case Adult :
                printf("A ");
                break;
            default :
                printf("  ");
                /* no break */
            }
        }
        printf("\n");
    }
    return;
}

int LayerPrint(int run, int rep, HQuality *H, int Xsize, int Ysize, FILE *out){
    /* prints the layer values out to a file */
    int i,j,k;
    for(j=0;j<Ysize;j++){
        for(k=0;k<Xsize;k++){
            i = j * Xsize + k;
            fprintf(out,"%-3d %-3d %-3d %-3d %-3d\n",run,rep,j,k,H[i]);
        }
    }
    return 0;
}

float LayerStats(int run, int rep, int returntype, HQuality *H, int Xsize, int Ysize, FILE *out){
    /* calculates a variety of statistics about the landscape layer */
    int i,j;
    long sum=0, sum2=0, N = Xsize * Ysize;
    int min = 101, max = 0, t;
    float mean, stddev, returnvalue;

    if (N <= 0) return 1;       /* error */

    /* loop for sum, sum2, min, max */
    for(i=0;i<(Xsize*Ysize);i++){
        t = H[i];
        sum += t;
        sum2 += t*t;
        if (t > max) max = t;
        if (t < min) min = t;
    }
    mean = (float) sum / N;
    stddev = (float) sum2 - (sum * sum)/N;
    stddev /= (N - 1);          /* calculates the variance */
    stddev = (float) sqrt((double) stddev);

    //    fprintf(out,"%d %d %f %f\n",min,max,mean,stddev);
    if (!returntype) returnvalue = mean;
    else returnvalue = quantile(H, N, returntype);
    return returnvalue;
}

int quantile(int *H, int N, int q)
{
    /* returns the qth rank value from the array H, which has N values unsorted */
    int *iHsorted;
    int returnvalue;
    int i;
    float iMax;

    iHsorted = (int *) malloc(sizeof(int)*N+1);

    for (i=0;i<N;i++) iHsorted[i] = H[i];

    qsort((void *)iHsorted,(size_t) N, sizeof(int), hqCompare);

    iMax = ((float)q/100)*N;

    for (i=0;i<=iMax;i++);  /* skip up to the qth value */

    returnvalue =  iHsorted[i];

    free(iHsorted);

    return returnvalue;
}

int hqCompare(const void *e1,  const void *e2)
{
    if (*((int *) e1) < *((int *) e2)) return -1;
    else if (*((int *) e1) == *((int *) e2)) return 0;
    else if (*((int *) e1) > *((int *) e2)) return 1;
    return -9999;

}



#ifdef __BORLANDC__
#pragma argsused
#endif
double LayerCorr(HQuality *H1, HQuality *H2, int Xsize, int Ysize, FILE *out){
    /* computes the aspatial correlation between two layers */
    unsigned int i;
    unsigned long productsum = 0, sum1=0, sum2=0, N = Xsize * Ysize;
    unsigned int t1, t2;
    double r;

    if (N <= 0) return 1;       /* error */

    /* loop for sum, sum2, min, max */
    for(i=0;i<N;i++){
        t1 = H1[i];
        t2 = H2[i];
        productsum += t1*t2;
        sum1 += t1*t1;
        sum2 += t2*t2;
    }

    r = productsum / sqrt((double) sum1*sum2);

    return r;
}

void spatial_corr(LandScape *L, float *r){
    /* *r is a pointer to a vector Xsize*Xsize that will be filled with
    correlation coefficients for that distance from L->H[] */

    return;
}

float Distance(LandScape *L, Location X1, Location X2){
    float x1, x2, y1, y2;

    y1 = (float) (X1 / L->Xsize);
    x1 = (float) (X1 % L->Xsize);
    y2 = (float) (X2 / L->Xsize);
    x2 = (float) (X2 % L->Xsize);

    if (L->Sides == HEX){
        if ((int)y1 % 2) x1 += HEXCOL;
        if ((int)y2 % 2) x2 += HEXCOL;
        y1 *= HEXROW;
        y2 *= HEXROW;
    }

    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

void TruncLand(HQuality **L, HQuality low[], HQuality high[], Parameters *Pm){
    int i,j;
    for(i=0;i<=Pm->NumLayers-1;i++){
        for(j=1;j<=(Pm->Xsize*Pm->Ysize - 1);j++){
            if (L[i][j] < low[i]) L[i][j] = low[i];
            if (L[i][j] > high[i]) L[i][j] = high[i];
        }
    }
    return;
}

void ClearLand(HQuality L[], HQuality value, Parameters *Pm){
    int i;
    //for(i=0;i<=(Pm->Xsize*Pm->Ysize - 1);i++){
    for(i=0;i<=(Pm->Xsize*Pm->Ysize);i++){
        L[i] = value;
    }
    return;
}

void WriteLand(HQuality **L, FILE *outfile, Parameters *Pm){
    int i,j;
    for(i=0;i<=Pm->NumLayers-1;i++){
        for(j=0;j<=(Pm->Xsize*Pm->Ysize-1);j++){
            if ((fwrite(&L[i][j],sizeof(HQuality),1, outfile)) != 1){
                error("Error writing land to output");
            }
        }
    }
    return;
}

int MakeFractalLandscape(LandScape *L, LParameters *LPm)
{
    float **X = NULL;
    int N, i, j, maxlevel;
    float maximum = 0;
    float minimum = 0;
    float slope;
    float constant;
    int temp;
    HQuality *H = NULL;

    /* first figure out how many powers of two there are */
    maxlevel = (int) floor(log((double)L->Ysize) / log((double)2));

    /* then add one, to allow for extra space to accomadate landscape sizes other than
    powers of two */
    maxlevel += 1;
    N = (int) pow((double)2,maxlevel);

    /* now allocate the needed space for the fractal algorithm */
    X = matrix(0,N,0,N);
    N = MidPointFM2D(X,maxlevel,LPm->fracvar,LPm->Hdim,0);

    /* LOOP for maximum, minimum, and average, if L->Xsize isn't 2^n + 1 then this is
    taking only a portion of the upper left hand corner of the landscape */
    for(i = 0; i < L->Xsize; i++)
    {
        for(j = 0; j < L->Ysize; j++)
        {
            if (X[i][j] < minimum) minimum = X[i][j];
            if (X[i][j] > maximum) maximum = X[i][j];
        }
    }

    slope = (float) MAXHAB / (maximum - minimum);
    constant = (float) MAXHAB - (slope * maximum);

    H = L->H[0];

    /* initialize landscape layer with transformed values */
    for(i = 0; i < L->Xsize; i++)
    {
        for (j = 0; j < L->Ysize; j++)
        {
            temp = i*L->Xsize + j;
            H[temp] = (HQuality) floor(X[i][j]*slope + constant + 0.01);
        }
    }

    free_matrix(X,0,N,0,N);

    return N;
}

int MidPointFM2D(float **X, int maxlevel, float sigma, float H,
    int addition)
    /* assumes that random number generator already seeded
    Parameters:
    X[][] double indexed array of real numbers to hold landscape
    maxlevel = maximal # of recursions (N = 2^maxlevel)
    sigma = initial standard deviation
    H = determines fractal dimension D = 3 - H
    addition = boolean parameter to turn random additions on/off

    algorithm taken from Peitgen HO, Saupe D 1988 The science of fractal images
    Springer Verlag, New York. Chapter 2.

    */
{
    int i;
    int N = (int) pow((double)2,maxlevel);  /* size of array being filled */
    int stage;
    float temp;
    const float onehalf = 0.5f;

    float delta = sigma;    /* standard deviation */

    int x, y, y0, D, d; /* array indexing variables */

    /* declare interpolation functions */
    float f3(float delta, float x0, float x1, float x2);
    float f4(float delta, float x0, float x1, float x2, float x3);

    /*    if (N > MAXN) error("Attempted fractal landscape too big.");*/

    X[0][0] = delta * gasdev();
    X[0][N] = delta * gasdev();
    X[N][0] = delta * gasdev();
    X[N][N] = delta * gasdev();

    D = N;
    d = N/2;

    for (stage = 1; stage <= maxlevel; stage++)
    {
        /* going from grid type 1 to grid type 2 */

        temp = H * onehalf;
        temp = pow(onehalf, temp);
        delta *= temp;

        /* interpolate and offset points */
        for (x=d; x<=N-d; x += D)
        {
            for (y=d; y<=N-d; y += D)
            {
                X[x][y] = f4(delta,X[x+d][y+d],X[x+d][y-d],X[x-d][y+d],X[x-d][y-d]);
            }
        }

        /* going from grid type 2 to grid type 1 */
        delta = delta * pow(onehalf, onehalf*H);

        /* interpolate and offset boundary points */
        for (x = d; x <= N-d; x += D)
        {
            X[x][0] = f3(delta,X[x+d][0],X[x-d][0],X[x][d]);
            X[x][N] = f3(delta,X[x+d][N],X[x-d][N],X[x][N-d]);
            X[0][x] = f3(delta,X[0][x+d],X[0][x-d],X[d][x]);
            X[N][x] = f3(delta,X[N][x+d],X[N][x-d],X[N-d][x]);
        }

        /* interpolate and offset interior points */
        for (x=d; x<=N-d; x += D)
        {
            for (y=D; y<=N-d; y += D)
            {
                X[x][y] = f4(delta,X[x][y+d],X[x][y-d],X[x+d][y],X[x-d][y]);
            }
        }
        for (x=D; x<=N-d; x += D)
        {
            for (y=d; y<=N-d; y += D)
            {
                X[x][y] = f4(delta,X[x][y+d],X[x][y-d],X[x+d][y],X[x-d][y]);
            }
        }

        /* scale down resolution */
        D /= 2;
        d /= 2;

    }
    return N;
}

float f3(float delta, float x0, float x1, float x2)
{
    const float three = 3.0f;
    return (x0 + x1 + x2)/three + delta * gasdev();
}

float f4(float delta, float x0, float x1, float x2, float x3)
{
    const float four = 4.0f;
    return (x0 + x1 + x2 + x3)/four + delta * gasdev();
}

