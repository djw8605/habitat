//---------------------------------------------------------------------------

#ifndef queueH
#define queueH
//---------------------------------------------------------------------------
/*************************************************************************
This file contains the definitions necessary to implement a queue of
individuals in the habitat model.
Last Modified: 30/04/96
Written by: Drew Tyre, Dept. of Environmental Science and Management,
						University of Adelaide, Roseworthy Campus, Roseworthy 5371 SA
            Ph: (08) 303 7931 Fax: (08) 303 7956
            email: dtyre@roseworthy.adelaide.edu.au
*************************************************************************/

#include "habitat.h"

struct Q{
  Individual *Head;
  Individual *Tail;
};

typedef struct Q Queue;

Queue *InitQueue(void);
Individual *EnQueue(Queue *,Individual *);
Individual *DeQueue(Queue *);
void DestroyQueue(Queue **);

#endif
