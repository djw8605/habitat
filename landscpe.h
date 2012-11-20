//---------------------------------------------------------------------------

#ifndef landscpeH
#define landscpeH
//---------------------------------------------------------------------------
/*************************************************************************
This file implements a landscape with variable geometry
Last Modified: 15/3/98
Written by: Drew Tyre, Dept. of Environmental Science and Management,
                        University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "util98.h"
#include "habitat.h"

LandScape *InitLandScape(Parameters *Pm, LParameters *LPm);
int DestroyLandscape(LandScape **L);
int GetDirection(LandScape *L);
HQuality GetValue(LandScape *L, int Layer, Location Where);
Location GetNextCell(LandScape *L, Location Where);
int GetCellinDir(LandScape *L, Location Where, int direction);
Location GetCell(LandScape *L, int force);
int GetNeighbours(LandScape *L, Location Where, HQuality *hq, int Force);
void printmap(LandScape *L);
int LayerPrint(int run, int rep, HQuality *H, int Xsize, int Ysize, FILE *out);
float LayerStats(int run, int rep, int returntype, HQuality *H, int Xsize, int Ysize, FILE *out);
double LayerCorr(HQuality *H1, HQuality *H2, int Xsize, int Ysize, FILE *out);
float Distance(LandScape *L, Location X1, Location X2);
void TruncLand(HQuality **L, HQuality low[], HQuality high[], Parameters *Pm);
void ClearLand(HQuality L[], HQuality value, Parameters *Pm);
void WriteLand(HQuality **L, FILE *outfile, Parameters *Pm);

#define A 0.5774
#define TWOA 1.1547
#define B (float) (2/3)

/* maximum size of fractal landscape */
#define MAXN 128
#define MAXMAXLEVEL 7

#define cAVERAGE 0x00

// Parameters for offset (had to be moved to the .h because it is going to be used in the habitat.c)
#define MAXBOUNDARYTYPE 12
#define MAXDIRECTION FULLSQR

#endif
 
