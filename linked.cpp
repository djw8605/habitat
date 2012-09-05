//---------------------------------------------------------------------------
#include "linked.h"

//---------------------------------------------------------------------------

/*************************************************************************
This file implements a linked list of individuals.

Last Modified: 10/05/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "util98.h"
#include "habitat.h"
#include "linked.h"

LinkList *InitList(void){
/* create a doubly linked list of individuals with a sentinel */
	LinkList *NewList = NULL;

  if ((NewList = (LinkList *) malloc(sizeof(LinkList))) != NULL){
  	NewList->NextInd = NewList;
    NewList->PrevInd = NewList;
    NewList->Age = -1;					/* Age < 0 identifies the sentinel */
  }
  else{
  	error("Insufficient Memory in InitList");
  }
  return NewList;
}

Individual *InsertList(LinkList* L, Individual *I){
/* insert an individual at the head of the list */
	I->NextInd = L->NextInd;
  L->NextInd->PrevInd = I;
  L->NextInd = I;
  I->PrevInd = L;
  return I;
}

Individual *DeleteList(LinkList *L, Individual *I){
/* delete a specific individual from the list */
/* assumes that I is part of a linked list */
	if (I == NULL) error("bad individual pointer in DeleteList");
	I->PrevInd->NextInd = I->NextInd;
  I->NextInd->PrevInd = I->PrevInd;
  I->PrevInd = NULL;
  I->NextInd = NULL;
  return I;
}

void DestroyList(LinkList **L){
/* destroy the list and set the sentinel pointer to NULL */
	/* delete all individuals left in the list */
	while((*L)->NextInd != *L) free((void *) DeleteList(*L,(*L)->NextInd));
  /* free memory for sentinal */
	free((void *) *L);
  *L = NULL;
  return;
}
