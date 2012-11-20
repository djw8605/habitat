//---------------------------------------------------------------------------

#include "ecologst.h"

//---------------------------------------------------------------------------

/*************************************************************************
This file contains the implementation of the ecologist!
Last Modified: 18/03/97
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
						Ph: (08) 303 7931 Fax: (08) 303 7956
						email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "util98.h"
#include "linked.h"
#include "queue.h"
#include "glider.h"
#include "habitat.h"
#include "ecologst.h"
#include "landscpe.h"

/* following extern declarations of functions and variables in habitat.c are
	for calculating sensitivities of growth matrix entries */
float tGrowthRate(float x);
float br[MAXAGECLASS+1] = {0};
float dr[MAXAGECLASS+1] = {0};
int Index, MSize;
float *Vector = NULL;

/* StratifiedSample is an array of territory locations chosen to meet some
	stratification criterion */
#define MAXSAMPLE 1090
#define MAXCATEGORIES 10

int SampleSites[MAXSAMPLE + 1];

int RegularSites(LandScape *L, Parameters *Pm)
{
/* picks Pm->nSamples sites on a regular grid off the landscape, and stores the
	addresses in SampleSites array */
	int i = 0,j,k;
   int totalj,totalk;
	Location Where;

	totalj = Pm->Ysize / Pm->SampleSep;
	totalk = Pm->Xsize / Pm->SampleSep;

	for(j = Pm->SampleSep;j < Pm->Ysize; j += Pm->SampleSep){
   for(k = Pm->SampleSep;k < Pm->Xsize; k += Pm->SampleSep){
   	Where = j*Pm->Ysize + k;
      SampleSites[i++] = Where;
	}
	}

	return 0;
}

int RandomSites(LandScape *L, Parameters *Pm)
{
/* picks Pm->nSamples sites randomly off the landscape, and stores the
	addresses in SampleSites array */

	int i = 0,j/*,k*/;
	/*int totalj,totalk;*/
	int Lsize;
	int *picked = NULL;
	/*Location Where;*/

	Lsize = Pm->Ysize*Pm->Xsize + 1;
	picked = (int *) calloc(Lsize, sizeof(int));

	for (i=0;i<Pm->nSamples;i++){
		while(picked[j=random(Lsize)]);
		picked[j] = 1;
		SampleSites[i] = j;
	}
    free(picked);
	return 0;
}

int StratifySites(LandScape *L, Parameters *Pm)
{
/* picks (Xsize/Samplesep - 1)^2 points stratified across the landscape
	and stores locations in SampleSites array */
	int i,j,k,C;
	int desired, unassigned;
	int Categories[MAXCATEGORIES + 1][4] = {{0}}; /* identify which categories each territory in */
	/* Categories [][0] is the number of territories in the category
		Categories [][1] is the number of desired samples
		Categories [][2] is the actual samples assigned to date
		Categories [][3] is an index of how many territories have been added to that
		category */
	int *TerritoryByCat[MAXCATEGORIES + 1] = {NULL};
	int pick;
	static int fileflag = 1;
	FILE *logfile;

	/* if this is the first run/rep, erase any previous file of that name */
	if (fileflag){
		fileflag = 0;
		remove(Pm->StratifyFile);
	}

	logfile = fopen(Pm->StratifyFile,"a+");

	if (Pm->nSamples > MAXSAMPLE) error("too many samples to stratify");

	/* this loop counts the territories in each of 10 evenly spaced categories */
	for (k = 0; k<(Pm->Xsize*Pm->Ysize); k++)
	{
		/* C should range from 0 to MAXCATEGORIES */
		C = GetValue(L,Pm->BirthLayer,k) / MAXCATEGORIES;
		/* force Q = 100 into the top category */
		if (C == MAXCATEGORIES) C = MAXCATEGORIES - 1;
		Categories[C][0]++;
	}

	/* calculate desired samples in each category */
	/* and allocate memory in TerritoryByCat */
	desired = Pm->nSamples / MAXCATEGORIES;
	unassigned = 0;
	for (k=0; k<MAXCATEGORIES; k++)
	{
		if (Categories[k][0] < desired)
		{
			unassigned += (desired - Categories[k][0]);
			Categories[k][1] = Categories[k][0];
		}
		else
		{
			Categories[k][1] = desired;
		}
		TerritoryByCat[k] = (int *) malloc((size_t)(Categories[k][0]+1)*sizeof(int));
	}
	/* loop over categories again, redistributing unassigned samples evenly */
	/* this will oversample lower categories by no more than 1 per category */
	do
	{
		for (k=0; k<MAXCATEGORIES; k++)
		{
			/* if category not full yet and samples remain */
			if ((Categories[k][0] > Categories[k][1]) && unassigned)
			{
				unassigned--;
				Categories[k][1]++;
			}
		}
	}while (unassigned > 0);

	/* loop over map, putting individual territories into the appropriate cat. */
	for (k = 0; k<(Pm->Xsize*Pm->Ysize); k++)
	{
		/* C should range from 0 to MAXCATEGORIES */
		C = GetValue(L,Pm->BirthLayer,k) / MAXCATEGORIES;
		/* force Q = 100 into the top category */
		if (C == MAXCATEGORIES) C = MAXCATEGORIES - 1;
		TerritoryByCat[C][Categories[C][3]++] = k;
	}

	j = 0;
	for (k=0; k<MAXCATEGORIES; k++)
	{
		if (Categories[k][0] == Categories[k][1])
		{	/* as many territories as samples */
			for (i = 0; i<Categories[k][1]; i++)
			{
				SampleSites[j++] = TerritoryByCat[k][i];
			}
		}
		else if (Categories[k][0] > Categories[k][1] )/* more territories than samples */
		{
			for (i = 0; i<Categories[k][1]; i++)
			{
				do
				{
					pick = random(Categories[k][3]);
				}while (TerritoryByCat[k][pick] < 0);
				SampleSites[j++] = TerritoryByCat[k][pick];
				TerritoryByCat[k][pick] = -1;
			}

		}
		else /* buggered up! less sites than samples */
		{
			error("fewer sites than samples in StratifySites");
		}

		free(TerritoryByCat[k]);
		TerritoryByCat[k] = NULL;
	}

	if (j > Pm->nSamples) error("too many samples selected in StratifySites");

	/* print out stratification results */
	for (k=0; k<MAXCATEGORIES; k++)
	{
		fprintf(logfile,"%-5d%-5d%-5d%-5d%-5d%-5d\n",Pm->Run,Pm->Rep,
							 k,Categories[k][0],Categories[k][1],Categories[k][3]);
	}

	fclose(logfile);

	return j;
}

#define N00 0
#define N01 1
#define N10 2
#define N11 3

int esttransitions(int flag, LandScape *L, Parameters *Pm)
{
	/* estimates transition probabilities on the landscape over the last 100 years or so */
	int i,j,k;
	static int counts[MAXSAMPLE + 1][4] = {{0}};
	static int prevobs[MAXSAMPLE+1] = {0};
	int occupied,habitat;
	static int fileflag = 1;
	FILE *logfile;

	/* if this is the first run/rep, erase any previous file of that name */
	if (fileflag){
		fileflag = 0;
		remove(Pm->TransitFile);
		logfile = fopen(Pm->TransitFile,"a+");
		fprintf(logfile,"Run Rep hq x y occupied N00 N01 N10 N11\n");
      fclose(logfile);
	}

	switch(flag)
	{
		case RESETTRANSITIONS :
		{
			for (i=0;i<MAXSAMPLE;i++)
			{
				counts[i][N00] = 0;
				counts[i][N01] = 0;
				counts[i][N10] = 0;
				counts[i][N11] = 0;
			}
			for (i=0;i<Pm->nSamples;i++)
			{
				prevobs[i] = (L->Map[SampleSites[i]]) ? 1 : 0;
			}
			break;
		}
		case SUMTRANSITIONS :
		{
			for (i=0;i<Pm->nSamples;i++)
			{
				occupied = (L->Map[SampleSites[i]]) ? 1 : 0;
				/*C = GetValue(L,0,SampleSites[i]); / MAXCATEGORIES;
				if (C == MAXCATEGORIES) C = MAXCATEGORIES - 1;      */
				if (prevobs[i])
				{
					if (occupied) counts[i][N11]++;
					else counts[i][N10]++;
				}
				else
				{
					if (occupied) counts[i][N01]++;
					else counts[i][N00]++;
				}
				prevobs[i] = occupied;
			}
			break;
		}
		case OUTPUTTRANSITIONS :
		{
			logfile = fopen(Pm->TransitFile,"a+");
			for (i = 0; i < Pm->nSamples; i++)
			{
				habitat = GetValue(L,Pm->SurvLayer,SampleSites[i]);
				occupied = (L->Map[SampleSites[i]] ? 1 : 0);
				k = SampleSites[i] % Pm->Xsize;
				j = SampleSites[i] / Pm->Xsize;
				fprintf(logfile,"%-5d%-5d%-5d%-5d%-5d%-5d%-6d%-6d%-6d%-6d\n",
					Pm->Run,Pm->Rep,habitat,k,j,occupied,
					counts[i][N00],counts[i][N01],counts[i][N10],counts[i][N11]);
			}

			fclose(logfile);
			break;
		}
	}
	return 0;
}

int census(LinkList *P, int AgeStruc[], int Competitors[], Parameters *Pm){
/* censuses the population */
/* returns the total number of individuals in the population */
	LinkList *CurInd;
	int i=0,TotalSize=0;

  /* clear age structure array */
  for(i=0;i<=Pm->MaxAgeClass;i++) AgeStruc[i] = 0;
  for(i=0;i<=7;i++) Competitors[i] = 0;

	/* loop over population counting individuals */
	for(CurInd = P->NextInd; CurInd->Age >= 0; CurInd = CurInd->NextInd){
		if (CurInd->Age < Pm->MaxAgeClass) AgeStruc[CurInd->Age]++;
   	else AgeStruc[Pm->MaxAgeClass]++;
      Competitors[CurInd->dispersaltype + (CurInd->selectortype * 3)]++;
		TotalSize++;
  }
	return TotalSize;
}	/* end function census */

int LocateAll(int run, int rep, int year, int gridsize, LinkList *P, LandScape *L){
/* prints out the location of each individual the population */
/* returns the total number of individuals in the population */
	LinkList *CurInd;
	static int fileflag = 0;
	static FILE *mapfile;

	if (!fileflag){
		mapfile = fopen("mapfile.txt","w");
		fileflag = 1;
	}

	/* loop over population counting individuals */
	for(CurInd = P->NextInd; CurInd->Age >= 0; CurInd = CurInd->NextInd){
		fprintf(mapfile,"%-5d%-5d%-5d%-5d%-5d%-5d%-5d\n",run, rep, year, CurInd->Age,
			CurInd->Where / gridsize, CurInd->Where % gridsize, L->H[0][CurInd->Where]);
	}
	return 0;
}	/* end function locateall */

int MapAll(int run, int rep, int year, int gridsize, LinkList *P, LandScape *L, Parameters *Pm){
/* prints out the habitat quality and occupancy of every cell on the landscape */
/* returns the total number of individuals in the population */
	LinkList *CurInd;
	static int fileflag = 0;
	static FILE *mapfile;
	int i, habitat, occupied, age;

	if (!fileflag){
		mapfile = fopen("mapfile2.txt","w");
		fileflag = 1;
	}

	/* loop over population counting individuals */
	for (i = 0; i<(Pm->Xsize*Pm->Xsize); i++){
		habitat = GetValue(L,Pm->SurvLayer,i);
		occupied = (L->Map[i] ? 1 : 0);
		if (occupied){
			age = L->Map[i]->Age;
		}
		else {
			age = 0;
		}
		fprintf(mapfile,"%-5d%-5d%-5d%-5d%-5d%-5d%-5d%-5d\n",run, rep, year, age,occupied,
		i / gridsize, i % gridsize, habitat);
	}
	return 0;
}	/* end function locateall */

int DeathWatch(int run, int rep, int year, LinkList *CurInd, LandScape *L){
/* prints out the history of each indivdual at death */
	static int fileflag = 0;
	static FILE *deathfile;
	
	if (!fileflag){ 	/* this code executed once the first time called */
		deathfile = fopen("deadfile.txt","w");
		fileflag = 1;
	}

	fprintf(deathfile,
		"%-5d%-5d%-5d%-5d%-5d%-5d%-5d%-5d%-7.2f%-7.2f%-7.2f%-5d%-5d%-5d%-5d%-5d%-5d%-5d\n",
		run, rep, year, L->H[0][CurInd->Where], CurInd->AgeClass, CurInd->Age,
		CurInd->Where, CurInd->BornHere, Distance(L,CurInd->Where,CurInd->BornHere),
		CurInd->avghabquality, CurInd->avgstepstounoccupied, CurInd->queuepos,
		CurInd->numbumps, CurInd->numbumped, CurInd->numtries,
		CurInd->numsteps + CurInd->totalsteps,
		CurInd->numbabies, CurInd->died);

	return 0;
}	/* end function deathwatch */

int TransectOccupancy(int run, int rep, int year, int gridsize, LandScape *L){
/* prints out the location of each individual the population */
/* returns the total number of individuals in the population */
	int i,j;
	int Where;
	int Occupied;
	static int fileflag = 0;
	static FILE *mapfile;

	if (!fileflag){
		mapfile = fopen("transect.txt","w");
		fileflag = 1;
	}

	/* loop over Landscape */
	for(i=0;i<L->Xsize;i+=5){
	for(j=0;j<L->Ysize;j++){
		Where = i*L->Xsize + j;
		if (L->Map[Where] != NULL) Occupied = 1;
		else Occupied = 0;
		fprintf(mapfile,"%-5d%-5d%-5d%-5d%-5d%-5d\n",run, rep, year,
						i,j,Occupied);
	}
	}
	return 0;
}	/* end function transectoccupancy */

int SpatialCV(int run, int rep, int year, int gridsize, LandScape *L){
/* prints out the location of each individual the population */
/* returns the total number of individuals in the population */
	int i;
	int Where;
	static int fileflag = 0;
	static FILE *mapfile;
	float density;
	int start, span = 0;
	int tempsum = 0;
	int midgrid;
	int lastrow, currow;

	if (!fileflag){
		mapfile = fopen("sptialcv.txt","w");
		fileflag = 1;
	}

	midgrid = (gridsize/2)-1;	/* put the transect in the middle */

	fprintf(mapfile,"%-5d%-5d%-5d",run, rep, year);

	/* loop over transect counting individuals */
	for(start=(gridsize/2)-1;start>=0;start--){
		span += 2;							/* increment the scale */
		Where = start*gridsize + start;	/* position of the first cell */
		lastrow = (start+span-1)*gridsize + start;
		if (L->Map[Where] != NULL) tempsum++;	/* add in first cell */
		currow = Where;
		for(i=1;i<span-1;i++){
			currow += gridsize;
			if (L->Map[Where+i] != NULL) tempsum++;
			if (L->Map[lastrow+i] != NULL) tempsum++;
			if (L->Map[currow] != NULL) tempsum++;
			if (L->Map[currow+span-1] != NULL) tempsum++;
		}
		if (L->Map[Where+span-1] != NULL) tempsum++;	/* add in last cell */
		if (L->Map[lastrow] != NULL) tempsum++;
		if (L->Map[lastrow+span-1] != NULL) tempsum++;
		density = (float) tempsum / (span*span);
		fprintf(mapfile,"%-6.3f",density);
	}
	fprintf(mapfile,"\n");
	return 0;
}	/* end function SpatialCV */

int report(Parameters *Pm, float AvgHabitat, float DispDeaths, float AgeDeaths, int TotalSize,
	int AgeStruc[], int Competitors[])
{
/* reports total size and agestruc to stdout */
	int i;
	int resident, invader;

	if (Pm->invasion)
	{
   	resident = Pm->dispersaltype + Pm->acceptfunction * 4;
		invader = Pm->invaderdt + Pm->invaderaf * 4;
		printf("%-5d%-5d",resident, invader);
	}

	printf("%-5d%-5d%-6d%-7.3f%-7.3f%-7.3f%-7d",
		Pm->Run, Pm->Rep, Pm->Year, AvgHabitat,
		DispDeaths, AgeDeaths, TotalSize);

	for(i=0;i<=Pm->MaxAgeClass;i++) printf("%-7d",AgeStruc[i]);

	if (Pm->competition || Pm->invasion)
	   for(i=0;i<6;i++) printf("%-7d",Competitors[i]);

	if (Pm->randomalpha) printf("%-8.5f",Pm->alpha);

	printf("\n");

  return 0;
}	/* end function freport() */

int freport(Parameters *Pm, float AvgHabitat, float DispDeaths, float AgeDeaths, int TotalSize,
	int AgeStruc[], int Competitors[])
{
/* reports total size and agestruc to stdout */
	int i;
	static int fileflag=0;
	static FILE *reportfile;
	int resident, invader;

	if (!fileflag){
		reportfile = fopen("report.txt","w");
		fileflag = 1;
	}

	if (Pm->invasion)
	{
		resident = Pm->dispersaltype + Pm->acceptfunction * 4;
      invader = Pm->invaderdt + Pm->invaderaf * 4;
      fprintf(reportfile,"%-7d%-7d",resident, invader);
   }

	fprintf(reportfile,"%-5d%-5d%-6d%-7.3f%-7.3f%-7.3f%-7d",
   	Pm->Run, Pm->Rep, Pm->Year, AvgHabitat,
		DispDeaths, AgeDeaths, TotalSize);

	for(i=0;i<=Pm->MaxAgeClass;i++) fprintf(reportfile,"%-7d",AgeStruc[i]);

	if (Pm->competition || Pm->invasion)
		for(i=0;i<6;i++) fprintf(reportfile,"%-7d",Competitors[i]);

	if (Pm->randomalpha) fprintf(reportfile,"%-8.5f",Pm->alpha);

	fprintf(reportfile,"\n");

  return 0;
}	/* end function freport() */

int reportvariance(int flag, int AgeStruc[],Parameters *Pm){
/* reports average and variance of age structure to variance.txt */
	int i;
	float variance;
	static int fileflag=0;
	static FILE *reportfile;
	static float averages[MAXAGECLASS + 1] = {0};	/* sums until the end */
	static float sumsqs[MAXAGECLASS + 1] = {0}; /* sum of squares */
	static int n = 0;	/* sample size */

	if (!fileflag){
		reportfile = fopen("variance.txt","w");
		fileflag = 1;
	}

	if (flag)
	{
		n++;
		for(i=0;i<=Pm->MaxAgeClass;i++)
		{
			averages[i] += AgeStruc[i];
			sumsqs[i] += AgeStruc[i]*AgeStruc[i];
		}
	}
	else	/* print out results */
	{
		fprintf(reportfile,"%-5d%-5d%-5d%-5d", Pm->Run, Pm->Rep, Pm->Year, n);
		for(i=0;i<=Pm->MaxAgeClass;i++)
		{
			variance = (sumsqs[i] - (averages[i]*averages[i])/n)/n;
			averages[i] /= n;
			fprintf(reportfile,"%-9.2f%-9.2f",averages[i],variance);
			averages[i] = 0;	/* reset static variables for next run */
			sumsqs[i] = 0;
		}
		n = 0;
		fprintf(reportfile,"\n");
	}
	return 0;
}	/* end function reportvariance() */

int sample_habitat(LandScape *L, Parameters *Pm){
/* samples the habitat at a particular grid size and writes output to a file */
	int i,j,k;
	int habitat, occupied, count = 0;
    int true_habitat, true_occupied;
	static int fileflag = 1;
	FILE *logfile;

	/* if this is the first run/rep, erase any previous file of that name */
	if (fileflag){
		fileflag = 0;
		remove(Pm->LogisticFile);
	}

	logfile = fopen(Pm->LogisticFile,"a+");

	/* if not stratifying, choose regular grid of sites */
	/* assumes that stratification already done elsewhere */
	/*if (!Pm->Stratify) RegularSites(L,Pm);*/
	/**** NOW DONE ALL IN HABITAT.C ******/

	for (i = 0; i<Pm->nSamples; i++)
	{
		true_habitat = GetValue(L,Pm->SurvLayer,SampleSites[i]);
        habitat = true_habitat + (int) floor(gasdev()*Pm->HabitatErr*true_habitat);

		occupied = true_occupied = (L->Map[SampleSites[i]] ? 1 : 0);
        if ((Pm->OccupancyErr > 0)&&(true_occupied))
        {
            if (rand0to1() < Pm->OccupancyErr)
                occupied = 0;
        }

		k = SampleSites[i] % Pm->Xsize;
		j = SampleSites[i] / Pm->Xsize;

		if (occupied) count++;

		fprintf(logfile,"%d %d %d %d %d %d %d %d %d\n",Pm->Run,Pm->Rep,
			(Pm->dispersaltype + Pm->acceptfunction*3),k,j,
            occupied, true_occupied,habitat, true_habitat);
	}

	fclose(logfile);

	return count;
}

int census_habitat(LandScape *L, Parameters *Pm){
/* writes the frequency of occupancy out for all habitat qualities */
/* also the average intrinsic growth rate for that quality */
	int i;
	int habitat, occupied;
	static int fileflag = 1;
	int count[MAXHAB+1] = {0};
	int occupants[MAXHAB+1] = {0};
	FILE *logfile;

	/* if this is the first run/rep, erase any previous file of that name */
	if (fileflag){
		fileflag = 0;
		remove(Pm->CHabitatFile);
	}

	logfile = fopen(Pm->CHabitatFile,"a+");

	for (i = 0; i<(Pm->Xsize*Pm->Xsize); i++)
	{
		habitat = GetValue(L,Pm->SurvLayer,i);
		occupied = (L->Map[i] ? 1 : 0);

		count[habitat]++;
		if (occupied) occupants[habitat]++;

	}

	for (i = 0; i<=MAXHAB; i++)
	{
		fprintf(logfile,"%d %d %d %d %d %f\n",Pm->Run,Pm->Rep,i,
			count[i],occupants[i],L->r[i]);
	}
	fclose(logfile);

	return 0;
}

int summary_habitat(LandScape *L, Parameters *Pm){
/* writes one line summary of habitat characteristics */
	int i;
	int habitat, occupied;
	static int fileflag = 1;
	int count[MAXHAB+1] = {0};
	int occupants[MAXHAB+1] = {0};
	int median_occ = 0, median_hab = 0, total_hab = 0;
	int sum_occ = 0, sum_hab = 0, total_occ = 0;
	FILE *logfile;

	/* if this is the first run/rep, erase any previous file of that name */
	if (fileflag){
		fileflag = 0;
		remove(Pm->CHabitatFile);
	}

	logfile = fopen(Pm->CHabitatFile,"a+");

	total_hab = (Pm->Xsize*Pm->Xsize);
	for (i = 0; i< total_hab; i++)
	{
		habitat = GetValue(L,Pm->SurvLayer,i);
		occupied = (L->Map[i] ? 1 : 0);

		count[habitat]++;
		if (occupied) {
			occupants[habitat]++;
			total_occ++;
		}

	}

	for (i = 0; (i<=MAXHAB)&&(total_occ); i++)
	{
		sum_occ += occupants[i];
		sum_hab += count[i];
		if ((((float) sum_occ / total_occ) > 0.5) && !median_occ) median_occ = i;
		if ((((float) sum_hab / total_hab) > 0.5) && !median_hab) median_hab = i;
		if (median_occ && median_hab) break;
	}

	fprintf(logfile,"%d %d %d %d %d %6.3f\n",Pm->Run, Pm->Rep, median_hab, median_occ,
		Pm->Qpivot, Pm->SurvSlope);
	fclose(logfile);

	return 0;
}

int sample_cohort(int Flag, LinkList *I, Parameters *Pm)
{
	static int total_cohort = 0;
	static int lx[MAXAGERECORDED+1] = {0};
	static int dx[MAXAGERECORDED+1] = {0};
	static int dispersal = 0;
	float DispRate = 0;
	float SurvRate[MAXAGECLASS + 1] = {0};
	float FecRate[MAXAGECLASS + 1] = {0};
	int nAdults = 0;
	int Age, i, j;
	int mx;
	float Remember, err, temppcDD;

	static int fileflag = 1;
	FILE *cohortfile;

	/* if this is the first run/rep, erase any previous file of that name */
	if (fileflag){
		fileflag = 0;
		remove("cohort.txt");
	}

	if (I)
	{
		if (I->Age > MAXAGERECORDED) Age = MAXAGERECORDED;
		else Age = I->Age;
	}
	else if ((Flag == MARK) || (Flag == BIRTHRECORD) || (Flag == DEATHRECORD))
	{
		error("bad individual pointer in sample_cohort");
	}

	if (total_cohort == INT_MAX) error("cohort overflow!");

	switch (Flag)
	{
		case RESETCOHORT :
			total_cohort = 0;
			dispersal = 0;
			for (i = 0; i <= MAXAGERECORDED; i++)
			{
				lx[i] = 0;
				dx[i] = 0;
			}
			break;
		case MARK :
			total_cohort++;
			I->marked = 1;
			break;
		case  BIRTHRECORD :
			lx[Age]++; break;
		case  DEATHRECORD :
			if (I->died == Aging) dx[Age]++;
			else dispersal++; break;
		case  COHORTOUTPUT :
			if ((cohortfile = fopen("cohort.txt","a+"))==NULL)
				error("Unable to open cohortfile in sample_cohort");
			/* mx is the number alive at the beginning of the interval */
			mx = total_cohort;
			/* print out dispersal oriented deaths */
			fprintf(cohortfile,"%d %d %d %d %d %d\n",
				Pm->Run, Pm->Rep, -1, 0, dispersal, mx);
			mx -= dispersal;
			for (i = 0; i <= MAXAGERECORDED; i++)
			{
				fprintf(cohortfile,"%d %d %d %d %d %d\n",
							Pm->Run, Pm->Rep, i, lx[i], dx[i], mx);
				mx -= dx[i];
			}
			fclose(cohortfile);

			break;
		case COHORTSUMMARY :
			/* prints a summary of actual lifehistory parameters for this run */
			if ((cohortfile = fopen("cohort.txt","a+"))==NULL)
				error("Unable to open cohortfile in sample_cohort");

			/* mx is the number alive at the beginning of the interval */
			mx = total_cohort;

			/* mortality rate during dispersal */
			DispRate = 1 - (float) dispersal / mx;
			mx -= dispersal;

			/* first loop over juvenile age classes */
			for (i = 0; (i < Pm->MaxAgeClass)&&(mx>0); i++)
			{
				SurvRate[i] = 1 - (float) dx[i] / mx;
				mx -= dx[i];
			}
			/* total up fecundity and survival over all years, throwing away last
				category */
			for (i = Pm->MaxAgeClass; (i < MAXAGERECORDED) && (mx > 0); i++)
			{
				SurvRate[Pm->MaxAgeClass] += 1 - (float) dx[i] / mx;
				FecRate[Pm->MaxAgeClass] += (float) lx[i];
				nAdults += mx;
				mx -= dx[i];
			}

			SurvRate[Pm->MaxAgeClass] /= (MAXAGERECORDED - Pm->MaxAgeClass);
			if (nAdults > 0) FecRate[Pm->MaxAgeClass] /= nAdults;

			/* print out rates and sample sizes on which they are based for variance
				calculations */
			fprintf(cohortfile,"%d %d %f %d %f %d", Pm->Run, Pm->Rep,
				FecRate[Pm->MaxAgeClass], nAdults, DispRate, total_cohort);

			mx = total_cohort - dispersal;
			for (i = 0; i <= Pm->MaxAgeClass; i++) {
				fprintf(cohortfile," %f %d",SurvRate[i],mx);
				mx -= dx[i];
				}

			/* toss dispersal deaths in with 1st years for growth rate calc.*/
			SurvRate[0] = 1 - (float) (dx[0] + dispersal) / total_cohort;

			/* calculate growth rate and sensitivities for all entries */
			if (nAdults > 0){
				fprintf(cohortfile," %f",IGR(Pm->MaxAgeClass+1,FecRate,SurvRate));

				/* br and dr are globals in habitat.c */
				for (j=0; j<=Pm->MaxAgeClass; j++){
					br[j] = FecRate[j];
					dr[j] = SurvRate[j];
				}
				Index = 2;
				Vector = dr;
				MSize = Pm->MaxAgeClass+1;

				for (i=0; i<=Pm->MaxAgeClass; i++){
					Index = i;
					Remember = Vector[Index];
					temppcDD = dfridr(tGrowthRate, Vector[Index], 0.07, &err);
					Vector[Index] = Remember;
					fprintf(cohortfile," %f",temppcDD);
				}
				Vector = br;
				Index = Pm->MaxAgeClass;
				Remember = Vector[Index];
				temppcDD = dfridr(tGrowthRate, Vector[Index], 0.07, &err);
				Vector[Index] = Remember;
				fprintf(cohortfile," %f\n",temppcDD);
			}
			else
				fprintf(cohortfile,"\n");

			fclose(cohortfile);
			break;
		default : error("bad flag value in sample_cohort");
		/* no break */
	}

	return total_cohort;
}

float tGrowthRate(float x){
/* changes index I in br or dr (all globals), then calls IGR with those
	vectors. Called by dfridr */

	Vector[Index] = x;
	return IGR(MSize,br,dr);
}
 
int GISexpert(LandScape *L, Parameters *Pm){
	/* outputs the habitat quality grid as a ESRI grid raster in .asc format */
	FILE *esrifile;
	/* coordinate system for sample_habitat has k = x and j = y*/
	int i,x,y;
	HQuality habitat;
	
	if ((esrifile = fopen(Pm->ESRIFile,"w"))==NULL)
		error("Unable to open ESRIFile in GISexpert : ecologst.cpp");
	fprintf(esrifile,"ncols %i\n",Pm->Xsize);
	fprintf(esrifile,"nrows %i\n",Pm->Ysize);
	fprintf(esrifile,"xllcorner 0\n");
	fprintf(esrifile,"yllcorner 0\n");
	fprintf(esrifile,"cellsize 1\n");
	fprintf(esrifile,"NODATA_value  -9999\n");
	for (y = Pm->Ysize-1; y>=0; y--){
		for (x = 0; x<Pm->Xsize; x++)
		{
			i = y * Pm->Xsize + x;
			habitat = GetValue(L,Pm->SurvLayer,i);
			fprintf(esrifile,"%i ",habitat);
		}
		fprintf(esrifile,"\n");
	}
	
	fclose(esrifile);
	return 0;
}
