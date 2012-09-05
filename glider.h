//---------------------------------------------------------------------------

#ifndef gliderH
#define gliderH
/*************************************************************************
This file contains declarations for greater gliders
Last Modified: 09/05/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/

#include "habitat.h"
#include "linked.h"

LinkList *InitPop(LinkList *L, Parameters *Pm);
LinkList *InitInvader(LinkList *L, LandScape *Land, Parameters *Pm,
	int adults);
Individual *Birth(Individual *Parent, Parameters *Pm);
int Death(Individual **I, float rate, LinkList *Population, LandScape *L);
float IGR(int N, float *br, float *dr);


//---------------------------------------------------------------------------
#endif
