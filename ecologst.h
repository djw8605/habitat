//---------------------------------------------------------------------------

#ifndef ecologstH
#define ecologstH
//---------------------------------------------------------------------------
/*************************************************************************
This file contains the declarations for the ecologist!
Last Modified: 10/05/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "util98.h"
#include "linked.h"
#include "queue.h"
#include "glider.h"
#include "habitat.h"
#include "landscpe.h"


#define COHORTOUTPUT 0x01
#define COHORTSUMMARY 0x13
#define DEATHRECORD 0x02
#define BIRTHRECORD 0x03
#define MAXAGERECORDED 20
#define MARK 0x04
#define RESETCOHORT 0x05

#define RANDOMSAMPLE 0x06
#define STRATIFIED 0x07
#define REGULAR 0x08

#define RESETTRANSITIONS 0x10
#define OUTPUTTRANSITIONS 0x11
#define SUMTRANSITIONS 0x12

int StratifySites(LandScape *L, Parameters *Pm);
int RegularSites(LandScape *L, Parameters *Pm);
int RandomSites(LandScape *L, Parameters *Pm);
int esttransitions(int flag, LandScape *L, Parameters *Pm);
int census(LinkList *P,int AgeStruc[], int Competitors[], Parameters *Pm);
int LocateAll(int run, int rep, int year, int gridsize, LinkList *P, LandScape *L);
int MapAll(int run, int rep, int year, int gridsize, LinkList *P, LandScape *L, Parameters *Pm);
int DeathWatch(int run, int rep, int year, LinkList *CurInd, LandScape *L);
int TransectOccupancy(int run, int rep, int year, int gridsize, LandScape *L);
int SpatialCV(int run, int rep, int year, int gridsize, LandScape *L);
int report(Parameters *Pm, float AvgHabitat, float DispDeaths, float AgeDeaths, int TotalSize,
	int AgeStruc[], int Competitors[]);
int freport(Parameters *Pm, float AvgHabitat, float DispDeaths, float AgeDeaths, int TotalSize,
	int AgeStruc[], int Competitors[]);
int census_habitat(LandScape *L, Parameters *Pm);
int sample_habitat(LandScape *L, Parameters *Pm);
int summary_habitat(LandScape *L, Parameters *Pm);
int reportvariance(int flag, int AgeStruc[], Parameters *Pm);
int sample_cohort(int Flag, LinkList *I, Parameters *Pm);
int GISexpert(LandScape *L, Parameters *Pm);
#endif
