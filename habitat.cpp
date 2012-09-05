//---------------------------------------------------------------------------

/****f*********************************************************************
Modified 9/18/2006 Drew Tyre
School of Natural Resources
University of Nebraska-Lincoln
416 Hardin Hall, East Campus
Lincoln, NE 68583-0974
phone: +1 402 472 4054 fax: +1 402 472 2946
email: atyre2@unl.edu
http://snr.unl.edu/tyre

This version should compile under the MinGW developer suite 
http:://www.mingw.org/



This file implements an individual based model of an age-structured
vertebrate. The animal is based on the Greater Glider, Petauroides volans.
Parameters also available for Mountain Brushtail possum/

In this version there is/are:
1) Dispersal is sequential
2) Habitat quality scaled to quantile (provided as parameter)
3) No Catastrophes
4) No Males
5) began adding new dispersal strategies.

Last Modified: 15/3/98
Written by: Drew Tyre, Dept. of Environmental Science and Management,
                        University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
        Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "util98.h"
#include "habitat.h"
#include "queue.h"
#include "linked.h"
#include "glider.h"
#include "ecologst.h"
#include "landscpe.h"

#define MAXTRIES 2
#define DEBUG_IT 0x01
#undef DEBUG_IT /* comment out this line to activate debug printing*/

int natality(LandScape *L, LinkList *Population, Queue *Dispersers, Parameters *Pm);
int dispersal(LandScape *L, LinkList *Population, Queue *Dispersers, Parameters *Pm);
int accept0(LandScape *L, Location NewCell, Individual *newbie, Parameters *Pm);
int accept1(LandScape *L, Location NewCell, Individual *newbie, Parameters *Pm);
Location disperse0(LandScape *L, Individual *newbie, Parameters *Pm);
Location disperse1(LandScape *L, Individual *newbie, Parameters *Pm);
Location disperse2(LandScape *L, Individual *newbie, Parameters *Pm);
Location disperse3(LandScape *L, Individual *newbie, Parameters *Pm);
Location disperse4(LandScape *L, Individual *newbie, Parameters *Pm);
int initfuzzytable(Parameters *Pm);
int fuzzydirection(int bd, Individual *p);
int mortality(LandScape *L, LinkList *Population, Parameters *Pm);
float fchop(float a, float min, float max);
float habitatvalue(HQuality sq, HQuality fq, int AgeClass);
/*float tGrowthRate(float x);*/

/* global array of deviations for "fuzzydirection" functions */
int deviations[MAXSIDES+1][MAXSIDES+1] = {{0}};

/* arrays of pointers to habitat acceptance functions */
int (*accept[])(LandScape *L, Location NewCell, Individual *newbie, Parameters *Pm)
    = {accept0, accept1, NULL};
Location (*disperse[])(LandScape *L, Individual *newbie, Parameters *Pm)
    = {disperse0, disperse1, disperse2, disperse3, disperse4, NULL};


/* survival transformation functions and matrices */
void MakeDTrnsfrmTable(Parameters *Pm, float AvgHabitat);
void MakeBTrnsfrmTable(Parameters *Pm, float AvgHabitat);
void OutputSurvivalProb(Parameters *Pm, LandScape *L);

int PatchLayer;
int PrintSurvivalProb;
int SamplePopulation =0;
int SampleLocation =0;
int SampleReps[] = {0,95,75,94,92,74,59};   /* reps closest to average for each run */

float DTable[MAXAGECLASS+1][MAXHAB+1];
float BTable[MAXAGECLASS+1][MAXHAB+1];

/* variables previously in model() function, now global in scope for this module */
Queue *Dispersers = NULL;
LinkList *Population = NULL;
LinkList *CurInd = NULL;
static LandScape *AdaBlock = NULL;
LandCell *test = NULL;
int year;
int TotalSize;
int NewBorns,DispDeaths,AgeDeaths;
int PopStruc[MAXAGECLASS+1] = {0};
int Competitors[7] = {0};
/* setting numtypes > 1 ensures that runs with no competition will complete */
int numtypes = 6;
int iExtinct = 0;
float pcDDeaths = 0, pcADeaths = 0;
float AvgHabitat = 0;

int main(int argc, char *argv[]){
	Parameters Pm;
    LParameters LPm;
	
	load_parameters(&Pm, &LPm);

	init_run(&Pm,&LPm);
	Pm.Run = 1;
	init_rep(&Pm,&LPm);
	Pm.Rep = 0;
	
	while((TotalSize > Pm.QThreshold) && (Pm.Year < Pm.MaxYears))
	{
		step_rep(&Pm,&LPm);
	}

	close_rep(&Pm,&LPm);
	close_run(&Pm,&LPm);

}

int init_model(Parameters *Pm, LParameters *LPm){
  return 0;
}

int init_run(Parameters *Pm, LParameters *LPm){
  /*initialisation tasks that only need to be done once for a given parameter set */

  /* initialize deviation table for fuzzy direction finding */
  /* turn distribution table initialised for each individual in Birth */
  initfuzzytable(Pm);
  return 0;
}

int init_rep(Parameters *Pm, LParameters *LPm){
  /* initialisation tasks that need to be redone with every replicate */
  Pm->Year = 0;

  Dispersers = InitQueue();

  Population = InitList();

  InitPop(Population,Pm);         /* initialize the population */

  if (((Pm->Run == 1)&&(Pm->Rep == 0)) || (Pm->NewLand))
  {
      AdaBlock = InitLandScape(Pm,LPm);
  }

  AvgHabitat = LayerStats(Pm->Run, Pm->Rep,Pm->Qpivot,AdaBlock->H[0],Pm->Xsize,Pm->Ysize,NULL);
  MakeDTrnsfrmTable(Pm, AvgHabitat);
  MakeBTrnsfrmTable(Pm, AvgHabitat);

  /* put individuals on the landscape */
  if ((Pm->Xsize * Pm->Ysize) < Pm->AgeStruc[Pm->MaxAgeClass])
  {
      error("Too many individuals for this landscape");
  }

  for(CurInd = Population->NextInd; CurInd->Age >= 0; CurInd = CurInd->NextInd){
      CurInd->Where = GetCell(AdaBlock,1);
      CurInd->Quality = AdaBlock->H[0][CurInd->Where];
      AdaBlock->Map[CurInd->Where] = CurInd;
  }


  if (PrintSurvivalProb) {
      OutputSurvivalProb(Pm,AdaBlock);
  }

  /* get initial population size */
  TotalSize = census(Population, PopStruc, Competitors, Pm);

  /* reset cohort sampling */
  if (Pm->SampleCohort) sample_cohort(RESETCOHORT,NULL,Pm);

  /* obtain stratified sample of landscape always, whether used in habitat
      sampling or not */
  switch(Pm->Stratify){
      case STRATIFIED : StratifySites(AdaBlock,Pm); break;
      case REGULAR : RegularSites(AdaBlock,Pm); break;
      case RANDOMSAMPLE : RandomSites(AdaBlock,Pm); break;
      default : error("bad site sampling code in model()"); break;
  }

  return 0;
}

int step_rep(Parameters *Pm, LParameters *LPm){
  int i,j;    /* counter for map emptying loop */
  float  temppcDD;
  unsigned int temp;
  float err, Remember;    /* used to store error estimate from dfridr */


  Pm->Year++;


  /* population gives birth */
  NewBorns = natality(AdaBlock,Population,Dispersers,Pm);

#ifdef DEBUG_IT
  printf("dispersing...");
#endif

  /* disperse newborns */
  DispDeaths = dispersal(AdaBlock,Population,Dispersers,Pm);

#ifdef DEBUG_IT
  printf("finished\n");
#endif

  /* calculate per capita death rate during dispersal */
  if (NewBorns > 0)
  {
    temppcDD = (float) DispDeaths/NewBorns;
  }
  else
    temppcDD = 0.0;

  pcDDeaths = temppcDD;

#ifdef DEBUG_IT
  printf("aging...");
#endif
  /* age and kill population */
  AgeDeaths = mortality(AdaBlock,Population,Pm);
#ifdef DEBUG_IT
  printf("finished\n");
#endif

  /* calculate per capita death rate during aging */
  if (TotalSize > 0)
  {
    temppcDD = (float) AgeDeaths/(TotalSize + NewBorns - DispDeaths);
  }
  else
    temppcDD = 0.0;

  pcADeaths = temppcDD;

  /* census the population */
  TotalSize = census(Population, PopStruc, Competitors, Pm);

  /* and report results to report.txt */
  if (Pm->ReportEveryYear)
  {
    freport(Pm, AvgHabitat, pcDDeaths, pcADeaths, TotalSize, PopStruc, Competitors);
  }

  if (Pm->Year >= Pm->SampleHabitat) /* not extinct, and habitat sample desired */
  {
      /* if sample_habitat returns 0, set iExtinct flag */
      if (!sample_habitat(AdaBlock,Pm)) iExtinct = 1;
	  /* regardless, get the GIS expert to make a habitat map in ESRI grid format */
	  GISexpert(AdaBlock,Pm);	/* slightly inefficient if sampling done in more than one year */
  }

  if (Pm->SampleLocation){
    LocateAll(Pm->Run, Pm->Rep, Pm->Year, Pm->Xsize, Population, AdaBlock);
    TransectOccupancy(Pm->Run, Pm->Rep, Pm->Year, Pm->Xsize, AdaBlock);
  }

  /* the 1 is a flag that tells reportvariance to collect the sum */
  if (Pm->ReportVariance && (Pm->Year > 175))
      reportvariance(1, PopStruc, Pm);

  if (Pm->Year > 10075) SpatialCV(Pm->Run, Pm->Rep, Pm->Year, Pm->Xsize, AdaBlock);

  if (Pm->Year == 100) esttransitions(RESETTRANSITIONS,AdaBlock,Pm);
  if (Pm->Year > 100) esttransitions(SUMTRANSITIONS,AdaBlock,Pm);

  if (Pm->competition)
  {  /* break out if there is a winner */
    numtypes = 0;
    for (i=0;i<6;i++) if (Competitors[i] > 0) numtypes++;
  }

  if ((Pm->invasion) && (Pm->Year >= 200))
  {
    if (Pm->Year == 200)
    {   /* add in the new strategy */
      InitInvader(Population,AdaBlock,Pm,1);
      TotalSize = census(Population, PopStruc, Competitors, Pm);
    }
    numtypes = 0;
    for (i=0;i<6;i++) if (Competitors[i] > 0) numtypes++;
  }

  return 0;

}   /* end model() */

int close_rep(Parameters *Pm, LParameters *LPm)
{
  /* If Population went extinct, report its size and structure */
  if (TotalSize <= Pm->QThreshold)
  {
      iExtinct = 1;
  }

  if (Pm->CensusHabitat) summary_habitat(AdaBlock,Pm);
  
  if (Pm->MapAll) MapAll(Pm->Run, Pm->Rep, Pm->Year, Pm->Xsize, Population, AdaBlock, Pm);

  /* output transition probabilities by habitat category */
  esttransitions(OUTPUTTRANSITIONS,AdaBlock,Pm);

  /* we want to increment the reps also if extinction distribution desired */
  if ((!iExtinct) || (Pm->ExtinctTest))
  {
      /* increment # of reps 'cuz everything worked */
      Pm->Rep++;
  }

#ifdef DEBUG_IT
  printf("printing reports...");
#endif

  if (Pm->SampleCohort)
      printf("%-6d",sample_cohort(COHORTSUMMARY,NULL,Pm));

  /* report final state regardless of extinction etc. */
  freport(Pm, AvgHabitat, pcDDeaths, pcADeaths, TotalSize, PopStruc, Competitors);
  report(Pm, AvgHabitat, pcDDeaths, pcADeaths, TotalSize, PopStruc, Competitors);

  /* the 0 tells reportvariance to calculate and print out the results */
  if (Pm->ReportVariance && !iExtinct)
      reportvariance(0, PopStruc, Pm);



#ifdef DEBUG_IT
  printf("finished\n");
#endif

#ifdef DEBUG_IT
  printf("rep: %-d run: %-d SampleReps: %-d\n",Pm->Rep,Pm->Run,SampleReps[Pm->Run]);
#endif

  if (Pm->LocateAll) LocateAll(Pm->Run, Pm->Rep, Pm->Year, Pm->Xsize, Population, AdaBlock);

#ifdef DEBUG_IT
  printf("tidying up...");
#endif

 /* remove dispersal queue */
  DestroyQueue(&Dispersers);

  /* loop over map and set all pointers to zero */
  /* do this to reuse the same landscape */
  for(int i=0;i<=(Pm->Xsize*Pm->Ysize);i++){
    AdaBlock->Map[i] = NULL;
  }

  /* remove population list */
  DestroyList(&Population);

  /* destroy landscape on last run, or every run */
  if (Pm->NewLand)
    DestroyLandscape(&AdaBlock);

#ifdef DEBUG_IT
  printf("finished\n");
#endif

 /* in an extinction distribution run, always return 0 */
  if (Pm->ExtinctTest) iExtinct = 0;

  return iExtinct;

}

int close_run(Parameters *Pm, LParameters *LPm)
{
  /* fuzzytable is initialised in init_run - but no memory to be freed */
  return 0;
}

int close_model(Parameters *Pm, LParameters *LPm)
{
  return 0;
}

void load_parameters(Parameters *Pm, LParameters *LPm)
{
  FILE *infile;
  unsigned randseed;
  time_t t;
  FILE *RandomFile;

  /* initialize parameters */
  if ((infile = fopen("habitat.in","r")) == NULL){
    return;
    /*("Unable to open input file in Habitat model");*/
  }

  rdsvar(infile,"UseRandomSeed",&(Pm->UseRandomSeed),INTEGER);
  rdsvar(infile,"RandomSeed",&(Pm->RandomSeed),INTEGER);
  rdsvar(infile,"competition",&(Pm->competition),INTEGER);
  rdsvar(infile,"invasion",&(Pm->invasion),INTEGER);
  rdsvar(infile,"acceptfunction",&(Pm->acceptfunction),INTEGER);
  rdsvar(infile,"dispersaltype",&(Pm->dispersaltype),INTEGER);
  rdsvar(infile,"MinSearch",&(Pm->MinSearch),INTEGER);
  rdsvar(infile,"AcceptSensitivity",&(Pm->AcceptSensitivity),REAL);
  rdsvar(infile,"BaseAccept",&(Pm->BaseAccept),REAL);
  rdsvar(infile,"invaderdt",&(Pm->invaderdt),INTEGER);
  rdsvar(infile,"invaderaf",&(Pm->invaderaf),INTEGER);
  rdsvar(infile,"MaxYears",&(Pm->MaxYears),INTEGER);
  rdsvar(infile,"ExtinctionTest",&(Pm->ExtinctTest),INTEGER);
  rdsvar(infile,"MaxReps",&(Pm->MaxReps),INTEGER);
/*    rdsvar(infile,"PrintSurv",&PrintSurvivalProb,INTEGER);*/
  rdsvar(infile,"NewLand",&(Pm->NewLand),INTEGER);
  rdsvar(infile,"Xsize",&(Pm->Xsize),INTEGER);
  rdsvar(infile,"Ysize",&(Pm->Ysize),INTEGER);
  rdsvar(infile,"Sides",&(Pm->Sides),INTEGER);
  rdsvar(infile,"QThreshold",&(Pm->QThreshold),INTEGER);
  rdsvar(infile,"SampleSep",&(Pm->SampleSep),INTEGER);
  rdsvar(infile,"SampleInt",&(Pm->SampleInt),INTEGER);
  rdsvar(infile,"nSamples",&(Pm->nSamples),INTEGER);
  rdsvar(infile,"Stratify",&(Pm->Stratify),INTEGER);
  rdsvar(infile,"DeathWatch",&(Pm->DeathWatch),INTEGER);
  rdsvar(infile,"SampleLocation",&(Pm->SampleLocation),INTEGER);
  rdsvar(infile,"ReportEveryYear",&(Pm->ReportEveryYear),INTEGER);
  rdsvar(infile,"ReportVariance",&(Pm->ReportVariance),INTEGER);
  rdsvar(infile,"CensusHabitat",&(Pm->CensusHabitat),INTEGER);
//    rdsvar(infile,"HabitatErr",&(Pm->HabitatErr),REAL);
//    rdsvar(infile,"OccupancyErr",&(Pm->OccupancyErr),REAL);
  rdsvar(infile,"CHabitatFile",Pm->CHabitatFile,STRING);
  rdsvar(infile,"SampleHabitat",&(Pm->SampleHabitat),INTEGER);
  rdsvar(infile,"SampleCohort",&(Pm->SampleCohort),INTEGER);
  rdsvar(infile,"SampleRep",&(Pm->SampleRep),INTEGER);
  rdsvar(infile,"LocateAll",&(Pm->LocateAll),INTEGER);
  rdsvar(infile,"MapAll",&(Pm->MapAll),INTEGER);
/*    rdsvar(infile,"RandomSensitivity",&RandomSensitivity,INTEGER);
  rdsvar(infile,"RandomFile",randvar,STRING);*/
  rdsvar(infile,"HabitatFile",Pm->HabitatFile,STRING);
  rdsvar(infile,"LogisticFile",Pm->LogisticFile,STRING);
  rdsvar(infile,"StratifyFile",Pm->StratifyFile,STRING);
  rdsvar(infile,"TransitFile",Pm->TransitFile,STRING);
  rdsvar(infile,"ESRIFile",Pm->ESRIFile,STRING);
  rdsvar(infile,"NumLayers",&(Pm->NumLayers),INTEGER);
/*    rdsvar(infile,"PatchLayer",&(PatchLayer),INTEGER);*/
  rdsvar(infile,"Qpivot",&(Pm->Qpivot),INTEGER);
  rdsvar(infile,"SurvLayer",&(Pm->SurvLayer),INTEGER);
  rdsvar(infile,"SurvSlope",&(Pm->SurvSlope),REAL);
  rdsvar(infile,"BirthLayer",&(Pm->BirthLayer),INTEGER);
  rdsvar(infile,"BirthSlope",&(Pm->BirthSlope),REAL);
/*    rdsvar(infile,"FMultiplier",&FecundityMultiplier,REAL);*/
  rdsvar(infile,"DispMort",&(Pm->DispMort),REAL);
  rdsvar(infile,"randomalpha",&(Pm->randomalpha),INTEGER);
  rdsvar(infile,"alpha",(&Pm->alpha),REAL);
  rdsvar(infile,"minalpha",(&Pm->minalpha),REAL);
  rdsvar(infile,"maxalpha",(&Pm->maxalpha),REAL);
  rdsvar(infile,"MaxAgeClass",&(Pm->MaxAgeClass),INTEGER);
  rdvector(infile,"BirthRate",Pm->BirthRate,REAL);
  rdvector(infile,"DeathRate",Pm->DeathRate,REAL);
  rdvector(infile,"AgeStruc",Pm->AgeStruc,INTEGER);

  fclose(infile);

  if ((infile = fopen("makeland.in","r")) == NULL){
      return;/*error("Unable to open landscape parameter file in Habitat model");*/
  }
  rdsvar(infile,"hdim",&(LPm->Hdim),REAL);
  rdsvar(infile,"fracvar",&(LPm->fracvar),REAL);
  rdvector(infile,"meanval",LPm->meanval,INTEGER);
  rdvector(infile,"varval",LPm->varval,REAL);
  rdvector(infile,"varbump",LPm->varbump,REAL);
  rdvector(infile,"lowcut",LPm->lowcut,INTEGER);
  rdvector(infile,"highcut",LPm->highcut,INTEGER);
  rdvector(infile,"numdents",LPm->numdents,INTEGER);
  fclose(infile);

  if (Pm->UseRandomSeed)
    randseed = Pm->RandomSeed;
  else
    randseed = (unsigned) time(&t);

  /* write the seed into a file, appending if it is already present */
  RandomFile = fopen("seedfile.txt","a");
  fprintf(RandomFile,"%i %s",randseed,ctime((time_t *) &randseed));
  fclose(RandomFile);

  /* set the seed */
  sgenrand(randseed);
  return;
}
/*float tGrowthRate(float x){
/* changes index I in br or dr (all globals), then calls IGR with those
    vectors. Called by dfridr */

/*  Vector[Index] = x;
    return IGR(MSize,br,dr);
} */

int natality(LandScape *L, LinkList *Population, Queue *Dispersers, Parameters *Pm){
 /* population gives birth */
    LinkList *CurInd = NULL;
    LinkList *Baby = NULL;
    int birthcount = 0;
    float birthrate;
    HQuality HabValue;
    for(CurInd = Population->NextInd; CurInd->Age >= 0; CurInd = CurInd->NextInd){
        if (Pm->BirthSlope > 0.0001) HabValue = GetValue(L,Pm->BirthLayer,CurInd->Where);
        else HabValue = MAXHAB;
        birthrate = BTable[CurInd->AgeClass][HabValue];
        if (rand0to1() <= birthrate){
            if (CurInd->marked) sample_cohort(BIRTHRECORD,CurInd,Pm);
            birthcount++;
            CurInd->numbabies++;
            Baby = Birth(CurInd, Pm);
            Baby->Where = CurInd->Where;
            Baby->avghabquality = habitatvalue(GetValue(L,Pm->SurvLayer,
                                         Baby->Where),HabValue,Pm->MaxAgeClass);
            EnQueue(Dispersers,Baby);
            if (Pm->SampleCohort)
            {
                if ((Pm->Year > (Pm->MaxYears - MAXAGERECORDED*2)) &&
                      Pm->Year <=(Pm->MaxYears - MAXAGERECORDED))
                {
                    sample_cohort(MARK,Baby,Pm);
                }
            }
        }
    }
    return birthcount;
}

int dispersal(LandScape *L, LinkList *Population, Queue *Dispersers, Parameters *Pm){
/* put newborns into linked list, and place them on the map */
/* this function now implements synchronous dispersal, so that individuals
   suffer no bias from their original location in dispersal queue? */
    Individual *Returned = NULL;
    Individual *Displaced = NULL;
    int dispersecount=0, repeatcount=0;
    Location NewCell;
    int killcount = 0;

    while((Returned = DeQueue(Dispersers)) != NULL){
        dispersecount++;
      /* if this is first trip thru queue, remember position */
        if (Returned->queuepos == 0) Returned->queuepos = dispersecount;

/*        printf("%i %i %f\n",Returned->BornHere,GetValue(L,0,Returned->Where),Returned->avghabquality);*/

        while(Returned){
            repeatcount++;
         Returned->numsteps++;

                if (rand0to1() <= Pm->DispMort){
                /* didn't survive the journey */
                        killcount++;
                        Returned->died = Dispersal;
/*                  if (Pm->Year == Pm->MaxYears){
                        DeathWatch(Pm->Run,Pm->Rep,Pm->Year,Returned,L);
                     }*/
                        if (Returned->marked)
                        {
                    sample_cohort(DEATHRECORD,Returned,Pm);
                     if (Pm->DeathWatch) DeathWatch(Pm->Run,Pm->Rep,Pm->Year,Returned,L);
                        }
                        free(Returned);
                        Returned = NULL;
                }
                else{
                    NewCell = (disperse[Returned->dispersaltype])(L,Returned,Pm);
                    switch ((accept[Returned->selectortype])(L,NewCell,Returned,Pm))
                    {   /*accept is a POINTER to a function */
                  case -1 : /* displace previous owner */
                     Displaced = L->Map[NewCell];
                            Displaced->numbumped++;
                            Returned->numbumps++;
                            DeleteList(Population,Displaced);
                            EnQueue(Dispersers,Displaced);
							L->Map[NewCell] = NULL;
                            /* no break */
                        case 1 :
                    /* newcell is not occupied and is suitable */
                            Returned->Where = NewCell;
                            L->Map[NewCell] = Returned;
                            Returned->Quality = L->H[0][NewCell];   /* remember local habitat quality */
                            InsertList(Population,Returned);
                            Returned = NULL;    /* move on to the next individual */
                            break;
                        case 0 :
                            Returned->Where = NewCell;  /* move thru and keep looking */
                    break;
                    }   /* end switch */
                }   /* else survived */
        }   /* end while loop */
    }   /* end while dispersers left */
    return killcount;
}

int accept0(LandScape *L, Location NewCell, Individual *newbie, Parameters *Pm)
{
    /* this rule is now a direct implementation of accept1 with AcceptSensitivity=0
       and BaseAccept ~ infinity */
    float tempavg = newbie->avgstepstounoccupied;
    int tempsteps = newbie->numsteps;
   int temptries;
    HQuality currentsq,currentfq;
   float tempquality = newbie->avghabquality, currentvalue;

    if (L->Map[NewCell])    /* if cell occupied... */
    {
        if ((L->Map[NewCell])->Age > newbie->Age)   /* ...by someone older */
            return 0;   /*reject occupied cell*/
    }

    /* only do the following if cell unoccupied or occupied by a newborn */
    temptries = ++newbie->numtries; /* only have 1 sample for step counter */
   /* numtries is always going to be the # of unoccupied habitat cells encountered */
    newbie->avgstepstounoccupied = tempavg * (float)(temptries - 1)/temptries;
    newbie->avgstepstounoccupied += (float) tempsteps /temptries;
   newbie->totalsteps += newbie->numsteps;
   newbie->numsteps = 0;             /* reset step counter */
    /* there is an extra sample for the habitat quality average; the natal cell */
   newbie->avghabquality = tempquality * (float) temptries / (temptries+1);
   currentsq = GetValue(L,Pm->SurvLayer,NewCell);
   currentfq = GetValue(L,Pm->BirthLayer,NewCell);
   currentvalue = habitatvalue(currentsq,currentfq,Pm->MaxAgeClass);
    newbie->avghabquality += currentvalue / (temptries + 1);

    if (L->Map[NewCell])    /* cell occupied by another newborn */
   {
        if (rand0to1() <= 0.5) return -1;   /* displace owner */
      else return 0; /* lost out, keep moving */
   }
   else /* empty, suitable cell */
   {
        return 1;   /* accept landscape cell */
    }

}   /* end function accept0 */

int accept1(LandScape *L, Location NewCell, Individual *newbie, Parameters *Pm)
{
/*  This function assumes gliders reject occupied territories, and accept
    territories where the expected LRS is higher than average, including the
    territory where they were born, and discounted by the average # of steps to
    reach an unoccupied territory. Only unoccupied territory qualities are
    incorporated into their average.

    new twist - parameter MinSearch specifies the minimum number of unoccupied territories
    that must be searched in order to begin accepting - akin to a "best-of-n" rule without
    the returning to the best part.

    new twist - feed the ln(currentvalue/expectedvalue) through a logistic transform to
    come up with a probability of acceptance. Parameter AcceptSensitivity controls the
    steepness of the transform 0 => accept all with .5, >> 100 is near perfect cutoff

    BaseAccept controls the acceptance rate when there is no sensitivity

*/
   float tempavg = newbie->avgstepstounoccupied;
   int tempsteps = newbie->numsteps;
   int temptries;
   float tempquality = newbie->avghabquality;
   HQuality currentsq, currentfq;
    float currentvalue = 0;
    float expectedvalue = 0;
    float paccept = 0.0;

    if (L->Map[NewCell])
    {
        if ((L->Map[NewCell])->Age > newbie->Age)
/*        printf("%i\n",NewCell);*/
        return 0;   /*reject occupied cell*/
   }

    /* does the following only for unoccupied cells or cell with newborns in them */
    temptries = ++newbie->numtries; /* only have 1 sample for step counter */
    /* numtries is always going to be the # of unoccupied habitat cells encountered */
    newbie->avgstepstounoccupied = tempavg * (float)(temptries - 1)/temptries;
    newbie->avgstepstounoccupied += (float) tempsteps /temptries;
    newbie->totalsteps += newbie->numsteps;  /* keep track of total # of steps */
    newbie->numsteps = 0;             /* reset step counter */
    /* there is an extra sample for the habitat quality average; the natal cell */
    newbie->avghabquality = tempquality * (float) temptries / (temptries+1);
    currentsq = GetValue(L,Pm->SurvLayer,NewCell);
    currentfq = GetValue(L,Pm->BirthLayer,NewCell);
    currentvalue = habitatvalue(currentsq,currentfq,Pm->MaxAgeClass);
    newbie->avghabquality += currentvalue / (temptries + 1);
    expectedvalue = newbie->avghabquality;
    if (newbie->avgstepstounoccupied <= 0) error("Avgstepstounoccupied less than 0");
    expectedvalue *= pow((1-Pm->DispMort),floor(newbie->avgstepstounoccupied));

    /* BaseAccept = 0 -=> paccept = 0.5 at origin,
        AcceptSensitivity = 0 -=> paccept = BaseAccept in all cells,
        AcceptSensitivity >> 100 -=> sharp cutoff at currentvalue = expected value */
    paccept = Pm->AcceptSensitivity * log(currentvalue/expectedvalue) + Pm->BaseAccept;
    paccept = exp(paccept) / (1 + exp(paccept));

/*    printf("%i %i %f %f %f %f %f\n",NewCell,currentsq,newbie->avghabquality,
        newbie->avgstepstounoccupied,currentvalue, expectedvalue, paccept);*/

    /* always sample at least Pm->MinSearch unoccupied cells before settling */
    if ((newbie->numtries > Pm->MinSearch) && (rand0to1() < paccept))
   {
        if (L->Map[NewCell])    /* cell occupied by another newborn */
        /* because cells occupied by older individuals rejected above */
      {
        if (rand0to1() <= 0.5) return -1;   /* displace owner */
            else return 0; /* lost out, keep moving */
        }
      else  /* empty, suitable cell */
      {
          return 1; /* accept landscape cell */
        }
    }
   else return 0;   /* reject unsuitable cell */

}

Location disperse0(LandScape *L, Individual *newbie, Parameters *Pm)
{   /* localized random walk dispersal */
    return GetNextCell(L,newbie->Where);
}

Location disperse1(LandScape *L, Individual *newbie, Parameters *Pm)
{  /* global dispersal */
    return GetCell(L,0);
}

Location disperse2(LandScape *L, Individual *newbie, Parameters *Pm)
{   /* localized hill climbing dispersal, with occupation of destination reducing
        habitat quality according to how many steps have been taken to reach
        that spot */
    HQuality neighbours[FULLSQR+1] = {0};
    int ties[FULLSQR+1] = {-1};
    int temp;
   int numties = 0;
   int sides;
   int i;
   int bestdirection = -1;
    int bestquality = L->H[Pm->SurvLayer][newbie->Where];   /* local quality */

   /* fill neighbours with the qualities of surrounding cells
    decrement the quality of occupied cells by the number of
        steps taken + 1, so they start out slightly pessimistic */
   sides = GetNeighbours(L,newbie->Where,neighbours,newbie->numsteps+1);

   /* loop for maximum quality, remember tied directions in ties */
   for (i = 0; i < sides; i++)
    {
        if (neighbours[i] > bestquality)    /* if better than any previous */
      {
        bestdirection = i;
            bestquality = neighbours[i];
         numties = 0;
         ties[0] = bestdirection;   /* reset tie array */
      }
        else if (neighbours[i] == bestquality)  /* if tied with a previous */
        {
            ties[++numties] = i;    /* add direction to ties array */
      }
    }

   if (numties && (bestdirection >= 0)) /* there is a tie but not with current cell */
   {
      temp = numties + 1;
        temp = random(temp);
        bestdirection = ties[temp];
        if ((bestdirection > sides) || (bestdirection < 0))
        error("Bad tie-breaker in disperse2()?\n");
   }
    else if (numties && (bestdirection < 0)) /* there is a tie with the current cell */
   {
        if (numties > 1)
      {
         temp = random(numties) + 1;
            bestdirection = ties[temp];
        }
        else
        bestdirection = ties[1];
    }

    if (bestdirection < 0)  /* no neighbour better than current */
    return GetNextCell(L,newbie->Where);  /* pick a direction at random */
   else                     /* move to best unoccupied cell */
        return GetCellinDir(L,newbie->Where,bestdirection);
}

Location disperse3(LandScape *L, Individual *newbie, Parameters *Pm)
{   /* localized hill climbing dispersal, with random deviations from beta dist */
    HQuality neighbours[FULLSQR+1] = {0};
   int ties[FULLSQR+1] = {-1};
   int temp;
    int numties = 0;
    int sides;
    int i;
    int bestdirection = -1, newdirection = -1;
    int bestquality = L->H[Pm->SurvLayer][newbie->Where];   /* local quality */

    /* fill neighbours with the qualities of surrounding cells
        decrement the quality of occupied cells by the number of
        steps taken + 1, so they start out slightly pessimistic */
    sides = GetNeighbours(L,newbie->Where,neighbours,newbie->numsteps+1);

    /* loop for maximum quality, remember tied directions in ties */
    for (i = 0; i < sides; i++)
    {
        if (neighbours[i] > bestquality)    /* if better than any previous */
        {
            bestdirection = i;
            bestquality = neighbours[i];
            numties = 0;
            ties[0] = bestdirection;    /* reset tie array */
        }
        else if (neighbours[i] == bestquality)  /* if tied with a previous */
        {
            ties[++numties] = i;    /* add direction to ties array */
        }
    }

    if (numties && (bestdirection >= 0)) /* there is a tie but not with current cell */
    {
        temp = numties + 1;
        temp = random(temp);
        bestdirection = ties[temp];
        if ((bestdirection > sides) || (bestdirection < 0))
            error("Bad tie-breaker in disperse2()?\n");
    }
    else if (numties && (bestdirection < 0)) /* there is a tie with the current cell */
    {
        if (numties > 1)
        {
            temp = random(numties) + 1;
            bestdirection = ties[temp];
        }
        else
            bestdirection = ties[1];
    }


    if (bestdirection < 0)  /* no neighbour better than current */
        return GetNextCell(L,newbie->Where);  /* pick a direction at random */
    else                        /* otherwise move fuzzily in best direction */
    {   /* but if edge encountered, move randomly */
        newdirection =
            GetCellinDir(L,newbie->Where,fuzzydirection(bestdirection,newbie));
        if (newdirection < 0) newdirection = GetNextCell(L,newbie->Where);
        return newdirection;
    }
}

Location disperse4(LandScape *L, Individual *newbie, Parameters *Pm)
{   /* localized hill climbing dispersal with probability weighted according to habitat quality
       relative to other neighbours
       this takes no regard of the current quality, and effectively breaks ties randomly
       if Pm->alpha = 0, all directions equally likely
       if Pm->alpha = 1, directions directly proportional to habitat quality proportion
       if Pm->alpha > 1, good directions get more weight */
    HQuality neighbours[FULLSQR+1] = {0};
    double p_direction[FULLSQR+1] = {0.0};
   int temp;
    int sides;
    int i;
    int bestdirection = -1, newdirection = -1;
    int sumquality = 0;
    float sum_p = 0;

    /* fill neighbours with the qualities of surrounding cells */
    sides = GetNeighbours(L,newbie->Where,neighbours,0);

    /* loop for sum */
    for (i = 0; i < sides; i++)
    {
        /* only valid directions will have qualities >= 0 */
        if (neighbours[i] >= 0) sumquality += ++neighbours[i];
        else neighbours[i] = 1; /* screen out the invalid directions */
        /* because negative neighbours will wreak havoc with pow() below
           and 0 quality directions SHOULD be ignored with 0 probability */
    }

    /* check the result for sense! */
    /* there is some kind of overflow problem here ... not completely fixed yet */
    if (sumquality <= 1)
    {
       for(i=0; i< sides; i++) printf(" %i",neighbours[i]);
       printf(" %i %i\n",sumquality,newbie->Where);
       error("bad sumquality in disperse4");
    }

    /* loop to get non-normalised probabilities */
    for (i = 0; i < sides; i++)
    {
        sum_p += p_direction[i] = pow((float) neighbours[i] / sumquality,Pm->alpha);
    }

    /* loop to normalise probabilities */
    for (i = 0; i < sides; i++)
    {
        p_direction[i] /= sum_p;
    }
    /* put in a buffer to catch overruns */
    p_direction[sides] = 1.0;

    /* choose a direction randomly with roulette wheel selection *
    /* but if edge encountered, move randomly */
    newdirection = GetCellinDir(L,newbie->Where,randompick(p_direction,sides));
    if (newdirection < 0) newdirection = GetNextCell(L,newbie->Where);
    return newdirection;

}   /* end disperse4 */

int initfuzzytable(Parameters *Pm)
{
    /* initializes the deviation table for SAW's */
    int initial, index;
   int nd;

    for (initial = 0; initial < Pm->Sides; initial++)
    {
        for (index = 0; index <= Pm->Sides; index++)
        {
            nd = index - Pm->Sides / 2; /* hopefully integer division */
            nd += initial;
            if (nd < 0) nd += Pm->Sides;    /* deviation is negative in this case */
            else if (nd >= Pm->Sides) nd -= Pm->Sides;  /* deviation positive */
            deviations[initial][index] = nd;
        }
    }
    return 0;
}   /* end initfuzzytable */

int fuzzydirection(int bd, Individual *p)
{
    /* applies a beta distribution to expected direction using lookup table */
    float key;
    int index = -1;
    int nd, deviation;

    key = rand0to1();       /* get a random number */
    while (p->tdist[++index] < key);    /* search table for next greater value */

    /* index must now be translated into a deviation from the bestdirection */

    return deviations[bd][index];

}   /* end fuzzydirection */

int mortality(LandScape *L, LinkList *Population, Parameters *Pm){
    /* apply mortality and ageing */
    LinkList *CurInd = NULL;
    int killcount = 0;
    int popcount = 0;
    float survrate;

    /* this loop terminating condition works because there is a sentinal at the end of the
      population list with age < 0. new individuals inserted at head of list */
    for (CurInd = Population->NextInd; CurInd->Age >= 0; CurInd = CurInd->NextInd){
        popcount++;
        survrate = DTable[CurInd->AgeClass][GetValue(L,Pm->SurvLayer,CurInd->Where)];

      /* survrate is actually the probability of dying, just to be confusing */
        if (rand0to1() < survrate){
            CurInd->died = Aging;
/*          if ((Pm->Year == Pm->MaxYears)&&(Pm->DeathWatch)){
            DeathWatch(Pm->Run,Pm->Rep,Pm->Year,CurInd,L);
         }*/
            if (CurInd->marked)
            {
            sample_cohort(DEATHRECORD,CurInd,Pm);
            if (Pm->DeathWatch) DeathWatch(Pm->Run,Pm->Rep,Pm->Year,CurInd,L);
            }
            Death(&CurInd,survrate,Population,L);   /* kill the indivdual off */
            killcount++;
        }   /* end if */
        else {
            CurInd->Age++;
            if (CurInd->AgeClass < Pm->MaxAgeClass){
                CurInd->AgeClass++; /* once adulthood reached this doesn't increment anymore */
            }   /* end if */
        }
    }   /* end for */
    return killcount;
}

void MakeDTrnsfrmTable(Parameters *Pm, float HabitatPivot){
    float curdeathrate;
    float temp;
    int i,j;

    if (Pm->SurvSlope > 0.0001)
   {
        for(i=0;i<=Pm->MaxAgeClass;i++){
            curdeathrate = (float) log((double) Pm->DeathRate[i]/(1-Pm->DeathRate[i]));
            for(j=0;j<=MAXHAB;j++){
                temp = -Pm->SurvSlope*(((float)j/HabitatPivot) - 1);    /* habitat values less than AVGHAB are positive */
                temp += curdeathrate;           /* add in the transformed survival rate */
                temp = (float) exp((double) temp);                  /* begin backtransformation */
                DTable[i][j] = temp / (temp + 1);   /* complete backtransformation and store */
            }
        }
   }
   else
    {
        for(i=0;i<=Pm->MaxAgeClass;i++){
            for(j=0;j<=MAXHAB;j++){
                DTable[i][j] = Pm->DeathRate[i];
            }
        }
    }

    return;
}

void MakeBTrnsfrmTable(Parameters *Pm, float HabitatPivot){
    float curbirthrate;
    float temp;
    int i,j;

    if (Pm->BirthSlope > 0.0001) {
        for(i=0;i<=Pm->MaxAgeClass;i++){
            temp = 1.99 * Pm->BirthRate[i];
            if (temp > 0.0){
                curbirthrate = (float) log((double) temp/(1-temp));
                for(j=0;j<=MAXHAB;j++){
                    temp = Pm->BirthSlope*(((float)j/HabitatPivot) - 1);    /* habitat values less than AVGHAB are positive */
                    temp += curbirthrate;           /* add in the transformed survival rate */
                    temp = (float) exp((double) temp);                  /* begin backtransformation */
                    temp = temp / (temp + 1);   /* complete backtransformation and store */
                    BTable[i][j] = temp / 2.0;
                }
            }
            else{
                for(j=0;j<=MAXHAB;j++){
                    BTable[i][j] = 0.0;
                }
            }
        }
    }
    else {
        for(i=0;i<=Pm->MaxAgeClass;i++){
            for(j=0;j<=MAXHAB;j++){
                BTable[i][j] = Pm->BirthRate[i];
            }
        }
    }
    return;
}

void OutputSurvivalProb(Parameters *Pm, LandScape *L){
    int x;
    HQuality HabVal;

    for(x = 0; x<(Pm->Xsize*Pm->Ysize); x++){
            HabVal = GetValue(L,Pm->SurvLayer,x);
        printf("%-5d%-5d%-5d%-7.4f%-7.4f%-7.4f\n",Pm->Run,Pm->Rep,HabVal,
                DTable[0][HabVal],
                DTable[1][HabVal],
                DTable[2][HabVal]);
    }
    return;
}

float fchop(float a, float min, float max){
    if (a > max){
    return max;
  }
  else{
    if (a < min){
        return min;
     }
  }
  return a;
}

float habitatvalue(HQuality sq, HQuality fq, int AgeClass)
{
    float temp;
   float S;
    float F;

   S = 1 - DTable[AgeClass][sq];
   F = BTable[AgeClass][fq];

    temp = (S*F) / (1-S);

    return temp;
}

