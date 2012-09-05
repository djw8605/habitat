//---------------------------------------------------------------------------

#include "glider.h"

//---------------------------------------------------------------------------

/*************************************************************************
This file implements greater gliders
Last Modified: 09/05/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
                        University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
                Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "util98.h"
#include "habitat.h"
#include "landscpe.h"
#include "linked.h"
#include "queue.h"
#include "glider.h"

Individual TheGodGlider =
    {NULL,NULL,0,0,Female,0,0,AVGHAB,AVGHAB,1.0,0,0,0,0,Alive,0,0,0,0,0,0,{0}};

void fillturndist(Individual *it, Parameters *Pm);
int maxvector(float *list, int N);

LinkList *InitPop(LinkList *L, Parameters *Pm){
/* creates a new population with the distribution given in AgeStruc */
  int j;
  int dt = 0, af = 0;
  int a;
  Individual *I;
  for(a=0;a<=Pm->MaxAgeClass;a++){
    for(j=0;j<(Pm->AgeStruc)[a];j++){
        I = Birth(&TheGodGlider, Pm);
        I->Age = a;
        I->AgeClass = a;
        if (Pm->competition)
        {
            if (dt > 2) dt = 0;
            I->dispersaltype = dt++;
            if (af > 1) af = 0;
            I->selectortype = af++;
        }
        else
        {
            I->dispersaltype = Pm->dispersaltype;
            I->selectortype = Pm->acceptfunction;
        }
        InsertList(L,I);
     }
  }
  return L;
}

LinkList *InitInvader(LinkList *L, LandScape *Land, Parameters *Pm, int adults){
/* adds a new batch of invaders to the population */
    int j;
    Individual *I;
    for(j=0;j<adults;j++){
        I = Birth(&TheGodGlider,Pm);
        I->Age = Pm->MaxAgeClass;
        I->AgeClass = Pm->MaxAgeClass;
        I->dispersaltype = Pm->invaderdt;
        I->selectortype = Pm->invaderaf;
        InsertList(L,I);
        I->Where = GetCell(Land,1);
      I->Quality = Land->H[0][I->Where];
      Land->Map[I->Where] = I;
    }
  return L;
}

Individual *Birth(Individual *Parent, Parameters *Pm){
/* Gives birth to a newborn and returns a pointer to it */
    Individual *NewInd;
    int i;

    if ((NewInd = (Individual *) malloc(sizeof(Individual))) != NULL){
        NewInd->PrevInd = NULL;
        NewInd->NextInd = NULL;
        NewInd->Age = 0;
        NewInd->AgeClass = 0;
        NewInd->Sex = Female;
        NewInd->Where = Parent->Where;
        NewInd->BornHere = Parent->Where;
        NewInd->Quality = Parent->Quality;
        NewInd->avghabquality = NewInd->Quality;
        NewInd->avgstepstounoccupied = 0;
        NewInd->numtries = 0;
        NewInd->numsteps = 0;
        NewInd->totalsteps = 0;
        NewInd->numbabies = 0;
        NewInd->died = Alive;
        NewInd->dispersaltype = Parent->dispersaltype;
        NewInd->selectortype = Parent->selectortype;
        NewInd->marked = 0;
        NewInd->numbumps = 0;
        NewInd->numbumped = 0;
        NewInd->queuepos = 0;
        if ((Parent->tdist[0] == 0.0) && (Pm->dispersaltype==3))
            fillturndist(NewInd,Pm);
        else
            for(i=0;i<=Pm->Sides;i++) NewInd->tdist[i] = Parent->tdist[i];
    }
    else{
        error("Insufficient Memory in Birth!");
    }
    return NewInd;
}

int Death(Individual **I, float rate, LinkList *Population, LandScape *L){
/* kills an individual glider returns 1 */

    LinkList *TempInd;

    TempInd = (*I)->PrevInd;    /* remember the previous individual */
    DeleteList(Population,*I);  /* delete the individual from the population */
    L->Map[(*I)->Where] = NULL;             /* empty map location */
    free((void *) *I);  /* destroy the individual */
    *I = TempInd;           /* set the current individual to be the previous one */

    return 1;
}

float IGR(int N, float *br, float *dr)
/* uses hqr algorithm to calculate the population growth rate for a given
    set of birth and death rates. N is the number of age classes. br is a
    vector of birth rates, dr is a vector of death rates */
{
    int i, j, biggest;
    float **a, *wr, *wi;
    float igr;

    a = matrix(1,N,1,N);

    /* construct the matrix, remembering that the input vectors are indexed
        from 0 */
    for (i=1;i<=N;i++) {
        for (j=1;j<=N;j++) a[i][j] = 0.0;
    }
    for (j=1;j<=N;j++) a[1][j] = br[j-1];
    for (i=2;i<=N;i++) a[i][i-1] = dr[i-2];
    a[N][N] = dr[N-1];

    wr=vector(1,N);
    wi=vector(1,N);

/*  printf("matrix:\n");
    for (i=1;i<=N;i++) {
        for (j=1;j<=N;j++) printf("%8.5f",a[i][j]);
        printf("\n");
    }*/

    hqr(a,N,wr,wi);

/*  printf("eigenvalues:\n");
    printf("%11s %16s \n","real","imag.");
    for (i=1;i<=N;i++) printf("%15f %14f\n",wr[i],wi[i]);*/

    biggest = maxvector(wr,N);
    if (wi[biggest] != 0.0){
        printf("eigenvalues:\n");
        printf("%11s %16s \n","real","imag.");
        for (i=1;i<=N;i++) printf("%15f %14f\n",wr[i],wi[i]);
        error("complex growth rate in IGR");
    }
    igr = wr[biggest];

    free_matrix(a,1,N,1,N);
    free_vector(wi,1,N);
    free_vector(wr,1,N);

    return igr;
}

#define VNEG    -3.37E+38F

/* operating on NR base 1 vectors!! */
int maxvector(float *list, int N){
/* scans 1 to N entries in the vector list and returns the index of the
    largest. returns -1 on error (no entry larger than MINFLOAT) */
    float max = VNEG;
    int maxi = -1;
    int i;

    for (i = 1; i<=N; i++) {
        if (list[i] > max){
            max = list[i];
            maxi = i;
        }
    }

    return maxi;
}

#undef VNEG
/* redeclare these global variables with extern storage directive in modules
    that need to call betafunc */
float alpha = 1.0;

float betafunc(float x){
/* returns the value of eqn 2. at x. BetaDist(a) in excel is the integral
    of  this function from 0 to a. The parameters alpha and beta are passed
    thru global variables; this version assumes alpha = beta*/
    static float t1, t2;
    float t3;
    static float oldalpha = 0.0;

    if (alpha <= 1) error("Bad alpha value in betafunc, glider.c");

    /* the following lines only have to be done once for a given value of
        alpha! many different x's can use the same values of t1 and t2 */
    if (alpha != oldalpha){
        t1 = gammln(2 * alpha);
        t2 = gammln(alpha);
        t2 *= 2;
        t1 -= t2;
        t2 = exp(t1);
        oldalpha = alpha;
    }
    t3 = t2 *  pow(x,(alpha - 1));
    t3 *= pow(1-x,(alpha - 1));

    return t3;
}

float betadist(float a, float x1, float x2){

    float t1;

    alpha = a;
    t1 = qromb(betafunc,x1,x2);
    return t1;

}

void fillturndist(Individual *it, Parameters *Pm){

    float x, y = 0;
    int index = 0;
    float sum = 0.0;

    for (index = 0, x = (float) 1 / (Pm->Sides * 2); x <= 1.0; index++, x += (float) 1 / Pm->Sides) {
         sum += betadist(Pm->alpha,y,x);
         it->tdist[index] = sum;
/*       printf("%d %f %f %f\n",index,y,x,it->tdist[index]);*/
         y = x;
    }
    it->tdist[index] = 1.0;
/*  printf("%d %f %f %f\n",index,y,x,it->tdist[index]);      */

    return;
}

