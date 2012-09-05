//---------------------------------------------------------------------------

#ifndef linkedH
#define linkedH
//---------------------------------------------------------------------------
/*************************************************************************
This file contains declarations for the Linked List of individuals
Last Modified: 30/04/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdlib.h>
#include "habitat.h"

typedef struct Ind LinkList;

LinkList *InitList(void);
Individual *InsertList(LinkList *, Individual *);
Individual *DeleteList(LinkList *, Individual *);
void DestroyList(LinkList **);

#endif
