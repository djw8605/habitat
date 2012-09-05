//---------------------------------------------------------------------------

#include "queue.h"

//---------------------------------------------------------------------------
/*************************************************************************
This file contains the code necessary to implement a queue of individuals
in the habitat model.
Last Modified: 10/05/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "util98.h"	/* This may cause problems on other computers */
#include "habitat.h"
#include "queue.h"

Queue *InitQueue(void){
/* Initializes a new Queue by allocating memory and setting pointers to Null */
	Queue *NewQueue;


  if ((NewQueue = (Queue *) malloc(sizeof(Queue))) != NULL)
	{
	  NewQueue->Head = NULL;
  	NewQueue->Tail = NULL;
  }
  else {			/* insufficient memory */
  	error("Insufficient Memory for New Queue!");
  }
  return NewQueue;
}	/* end InitQueue */

Individual *EnQueue(Queue *Q, Individual *I){
/* places I into Q. This queue cannot overflow. */

	if (Q->Tail){	 /* queue has at least one member */
  	Q->Tail->NextInd = I;
		I->PrevInd = Q->Tail;
	  Q->Tail = I;
  }
  else{					/* queue is empty */
  	Q->Head = I;	/* I assume that NextInd is NULL */
    Q->Tail = I;
  }
  return I;
}	/* end EnQueue */

Individual *DeQueue(Queue *Q){
/* If an individual is in the queue, returns it, otherwise returns NULL */
/* assumes that the last individual in the queue has NextInd = NULL */
	Individual *I = NULL;

	if (Q->Head){
  	I = Q->Head;
    Q->Head = Q->Head->NextInd;
  }

  if (!Q->Head) Q->Tail = NULL; /* Last Individual DeQueued */

  return I;
}	/* end Dequeue */

void DestroyQueue(Queue **Q){
/* destroys the queue, and sets the queue pointer to NULL */
  Individual *I = NULL;
  /* DeQueue all remaining individuals, and free them */
  while((I = DeQueue(*Q))!= NULL) free(I);
	free((void *) *Q);
  *Q = NULL;
	return;
}
