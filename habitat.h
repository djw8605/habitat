//---------------------------------------------------------------------------

#ifndef habitatH
#define habitatH

/**********************************************************************
    This file contains all the data structure definitions for the entire
  habitat model.

  In version 2, Locations on the map are represented with a single unsigned
  int. The maximum grid size is then 256 x 256.

Last Modified: 9/7/97
Written by: Drew Tyre, Dept. of Environmental Science and Management,
                University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
**********************************************************************/

#include <math.h>
#include <omp.h>
/* Number of neighbors for different neighborhoods */
#define ORTHO 4
#define HEX 6
#define FULLSQR 8
#define MAXSIDES FULLSQR

/* Maximum number of age categories that can be represented */
#define MAXAGECLASS 5

/* Habitat Quality maximum and average */
/***************************************************************
   NB!! if you change these values you must ensure that they are also
   changed in MAKELAND.IN and HABITAT.IN!!
****************************************************************/
#define MAXHAB 100
#define AVGHAB 50
#define MAXLAYERS 5

/* parameter structure */
struct Parm{
    int UseRandomSeed;
    unsigned RandomSeed;
   int competition;
   int invasion;
    int dispersaltype;
   int acceptfunction;
   int MinSearch;
   float AcceptSensitivity;
   float BaseAccept;
   int invaderdt;
   int invaderaf;
    int MaxYears;
   int ExtinctTest;
  int MaxReps;
  int Run;
  int Rep;
  int Year;
  int NewLand;
  int Xsize;
  int Ysize;
  int Sides;
  int QThreshold;
  int NumLayers;
  int Qpivot;
  int SurvLayer;
  float SurvSlope;
  int BirthLayer;
  float BirthSlope;
  int SampleSep;
  int SampleInt;
  int nSamples;     /* for stratified sampling */
  int Stratify;
  int DeathWatch;
  int SampleLocation;
  int ReportVariance;
  int ReportEveryYear;
  int CensusHabitat;
  int SampleHabitat;
  int SampleCohort;
  int SampleRep;
  int LocateAll;
  int MapAll;
  float HabitatErr;
  float OccupancyErr;
  char CHabitatFile[20];
  char HabitatFile[20];
  char LogisticFile[20];
  char StratifyFile[20];
  char TransitFile[20];
  char ESRIFile[20];
  float DispMort;
  int randomalpha;
  float alpha; /*TurnDist[MAXSIDES+1];*/
  float minalpha;
  float maxalpha;
  int MaxAgeClass;
  float DeathRate[MAXAGECLASS+1];
  float BirthRate[MAXAGECLASS+1];
  int AgeStruc[MAXAGECLASS+1];
};
typedef struct Parm Parameters;

typedef int HQuality;

/* Landscape parameter structure */
struct LParm {
    HQuality meanval[MAXLAYERS];
   float fracvar;
    float Hdim;
    float varval[MAXLAYERS];
    float varbump[MAXLAYERS];
    HQuality lowcut[MAXLAYERS];
    HQuality highcut[MAXLAYERS];
    int numdents[MAXLAYERS];
};
typedef struct LParm LParameters;

enum SexType {Male, Female};

/* The following is for Greater Gliders, I think */
enum AgeType {NewBorn=0, Juvenile=1, Adult=2};

enum DeathType {Alive = 0, Dispersal=1, Aging=2};
typedef unsigned int Location;

struct Ind{
    struct Ind *PrevInd;
    struct Ind *NextInd;
    int Age;    
    int AgeClass;   /* this is the age class */
    enum SexType Sex;
    Location Where;
    Location BornHere;
    int Quality;
    float avghabquality;
    float avgstepstounoccupied;
    int numtries;
    int numsteps;
    int totalsteps;
    int numbabies;
    enum DeathType died;
    int dispersaltype;
   int selectortype;
   int marked;
   int numbumps;
   int numbumped;
    int queuepos;
    float tdist[MAXSIDES + 1];  /* this is the turning distribution */

    // OpenMP lock
    omp_lock_t mp_indiv_lock;
};

typedef struct Ind Individual;

typedef Individual * LandCell;

struct LS{
    HQuality *H[MAXLAYERS]; /* array of pointer to int */
    LandCell *Map;      /* uses tntutil style array allocation */
   float r[MAXHAB+1];   /* array of intrinsic growth rates for each habitat quality */
    int unused;
    int Sides;  /* can take values of ORTHO, HEX, FULLSQR */
    int NumLayers;
    int Xsize;
    int Ysize;

    omp_lock_t *CellLocks;
};
typedef struct LS LandScape;

/* the functions that run the model! */
int init_model(Parameters *Pm, LParameters *LPm);
int init_run(Parameters *Pm, LParameters *LPm);
int init_rep(Parameters *Pm, LParameters *LPm);
int step_rep(Parameters *Pm, LParameters *LPm);
void load_parameters(Parameters *Pm, LParameters *LPm);
int close_rep(Parameters *Pm, LParameters *LPm);
int close_run(Parameters *Pm, LParameters *LPm);

/* STATE VARIABLES "published" by the habitat module */
extern int TotalSize;

#include "mt19937.h"
//---------------------------------------------------------------------------
#endif
