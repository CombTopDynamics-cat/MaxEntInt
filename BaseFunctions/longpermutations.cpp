/*********************************************************
   Class LongBasePermutation code of functions
   Used by mapentint.cpp and maxentallbounds.cpp
   This is another instance of class BasePermutations
   that allows lengths up to LONGMAXPERM (=1000).

   For comments we refer to the class BasePermutations
**********************************************************/#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/longpermutations.h"
#define DIR(x,y) ((x) < (y) ?  +1 : ((x) > (y) ? -1 : 0))

void LongBasePermutation::AllocateLongBasePermutation(unsigned int order)
{
  register unsigned int i;
  if (order < 2 || order > LONGMAXPERM){fprintf (stderr, "Error: Permutation length not allowed!!\n");  exit (1);}
  if(perm) free(perm);
  size = order; sm1 = size - 1; sm2 = size - 2;
  perm = (unsigned int *) malloc (sizeof (unsigned int) * size);
  if (!perm) { fprintf (stderr, "Memory problems initializing the permutations vector.\n"); exit (1); }
  for(i=0; i<size; i++) perm[i] = 0;
  return;
};

int LongBasePermutation::DefineManual(int order, int * map){
   register unsigned int i;

   LongBasePermutation::AllocateLongBasePermutation((unsigned int) order);
   for(i=0; i<size; i++) { if( map[i] < 1 || map [i] > order) { return 0; } perm[i] = (unsigned int) (map[i]-1); map[i] = 0; }
   for(i=0; i<size; i++) {if(map[perm[i]]) return 0; else map[perm[i]] = 1; }
   return 1;
}

int LongBasePermutation::IsCycle (void)
{ unsigned int i=0, r=1;
  while ((i=perm[i]) != 0 && r <= size ) {r++;}
  if( r != size ) return 0; else return 1;
}

int LongBasePermutation::IsMaximodal (void)
{
  register unsigned int i;
  int signo = DIR(perm[0],perm[1]), signn;

  for (i = 1; i < sm1; i++)
  {
    signn=DIR(perm[i],perm[i+1]);
    if (signo == signn) return 0;
    signo = signn;
  }
  return 1;
 }

int LongBasePermutation::ComputeNumberOfExceptionsIfMaximodal(void)
{ register unsigned int i; unsigned int startmin = 0, separation = (sm2)/2, count = 0;

  if(DIR(perm[0],perm[1]) < 0) startmin = 1;

  for(i=startmin; i < size; i+=2) if(perm[i] <= separation) count++;
  return separation-count+1;
}
