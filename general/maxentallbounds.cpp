/*
For every period between 4 and a maximum period entered as a parameter
three numbers are computed:
* The Misiurewicz-Nitecki bound
* The maximum entropy of permutations of period n
* The maximum entropy of n-cycles (including conjectured families)

These numbers are printed in table form with four columns.
The first column corresponds to the period.

It works for period smaller than 1000.

Usage: maximum-length(>4)
*/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/longpermutations.h"
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/maxent.h"

#define TOL 1.e-8

double MiNibound(unsigned int);
void maxperm(unsigned int, int *);
void maxcycle(unsigned int, int *);

int main (int argc, char *argv[])
{
  int map[LONGMAXPERM];
  unsigned int maxper=300, p;
  double entperm, entcyc;
  LongBasePermutation perm;
  IntervalMarkovMatrix MarMat (4);

  if( argc != 2 || atoi (argv[1]) < 4 || (maxper = atoi (argv[1]))  > LONGMAXPERM)
         Usage(argv[0], "maximum-length(>4)");

  for(p=4; p <= maxper ; p++){ MarMat.ResizeMarkovMatrix(p);
    maxperm(p,map);
    if(!perm.DefineManual(p, map)){
      fprintf (stderr, "\nIn permutation of lenght %u", p); ErrorAborting("unexpected problem");
    }
/*    if(!perm.IsMaximodal() || perm.ComputeNumberOfExceptionsIfMaximodal()){
       fprintf (stderr, "\nALARM: problems in permutation of lenght %u\n\n",p);
    } */
    MarMat.MapCompute(perm.perm);
    entperm = MarMat.SpectralRadius ();

    maxcycle(p,map);
    if(!perm.DefineManual(p, map)) {
      fprintf (stderr, "\nIn cycle of lenght %u", p); ErrorAborting("unexpected problem");
    }
/*    if(!perm.IsMaximodal() || perm.ComputeNumberOfExceptionsIfMaximodal() || !perm.IsCycle()){
       fprintf (stderr, "\nALARM: problems in cycle of lenght %u\n\n", p);
    } */
    MarMat.MapCompute(perm.perm);
    entcyc = MarMat.SpectralRadius ();

    printf("%u %.16lf %.16lf %.16lf\n", p, MiNibound(p), entperm, entcyc);
  }
  return 0;
}

/*
  See reference [13] (Misiurewicz, Nitecki) from the paper.
*/
#define PI 3.141592653589793238462643383279502884197
double MiNibound(unsigned int length){return ((2*length)/PI); }

/*
Functions to Generate the maximal permutations, cycles and conjectures

The formula for the "maxperm" of odd period can be found in
the reference [5] (Geller, Tolosa) from the paper. The case even come
from references [7] (Geller, Zhang) and [9] (King). In some cases we
take the reverse instead of the original permutation. They have the
same entropy.

*/

void maxperm(unsigned int length, int *map){
  unsigned int i;
  if (length<4) ErrorAborting("permutation length must be at least 4!!");
  if(length%2){ unsigned int l = (length-1)/4 ; l= 2*l;
    for(i=1; i <= length; i++){ int ii = i-1;
      if(i%2){ if(i <= length-l-2) map[ii] = length - l - i; else map[ii] = i - length + l + 1; }
      else   { if(i <= l) map[ii] = length - l + i - 1;      else map[ii] = length + l - i + 2; }
    }
  } else { unsigned int k = length/2;
    for(i=1; i <= length; i++){ int ii = i-1;
      if(i%2){ if(i <= k) map[ii] = k - i + 1; else map[ii] = i - k;       }
      else   { if(i <= k) map[ii] = k + i;     else map[ii] = 3*k - i + 1; }
    }
  }
}

/*
The cycles of maximum entropy for odd periods are the maximum
entropy permutations of the same period. This is proved in the
reference [5] (Geller, Tolosa) from the paper.
*/
void maxoodcycle(unsigned int period, int *map){
  if(period < 5 || !(period%2)) ErrorAborting("cycle length must be odd larger than 3!!");
  maxperm(period, map);
}

/*
The formula for the maximum entropy cycle of period 4k can be found in
the reference [11] (King, Strantzen) from the paper.
*/
void maxcycle4k(unsigned int period, int *map){
  unsigned int i, k=period/4;
  if(period == 4) { maxperm(period, map); return;}
  if(period < 8 || period%4) ErrorAborting("cycle length must be multiple of 4!!");

  for(i=1; i <= period; i++){ int ii = i-1;
    if(i%2){
      if(i <= k+1) map[ii] = 2*k - i + 1;
      else if(i <= 2*k+1) map[ii] = 2*k - i + 2;
      else if(i <= 3*k) map[ii] = i - 2*k - 1;
      else map[ii] = i - 2*k;
    } else {
      if(i <= k+1) map[ii] = 2*k + i;
      else if(i <= 2*k) map[ii] = 2*k + i - 1;
      else if(i <= 3*k) map[ii] = 6*k - i + 2;
      else map[ii] = 6*k - i + 1;
    }
  }
}

/*
These are the conjectured maximum entropy cycle of period 4k+2. See
Definition 5.1 from the paper.
*/
#define PAR(x) ((x)%2 ?  -1 : +1)
void conjecturemaxcycle4kp2(unsigned int period, int *map){
  unsigned int i, k=(period-2)/4;
  if(period == 6) {map[0]=3; map[1]=6; map[2]=2; map[3]=5; map[4]=1; map[5]=4; return ; }
  if(period < 10 || period%2 || (period-2)%4) ErrorAborting("cycle length must be of the form 4*k+2.!!");

  if(k%2){
    for(i=1; i <= period; i++){ int ii = i-1;
      if(i%2){
        if(i <= k-2) map[ii] = 2*k + 2 - i;
        else if(i <= k+2) map[ii] = 2*k + 1 - i;
        else if(i <= 2*k+3) map[ii] = 2*k + 4 - i;
        else if(i <= 3*k) map[ii] = i - 2*k - 3;
        else if(i == 3*k+2) map[ii] = k + 2;
        else map[ii] = i - 2*k - 1;
      } else {
        if(i <= k-1) map[ii] = 2*k + 1 + i;
        else if(i <= k+3) map[ii] = 2*k + i;
        else if(i <= 2*k+2) map[ii] = 2*k - 1 + i;
        else if(i <= 3*k+1) map[ii] = 6*k + 6 - i;
        else if(i == 3*k+3) map[ii] = 3*k + 2;
        else map[ii] = 6*k + 4 - i;
      }
    }
  } else { unsigned int kk=k/2;
    for(i=1; i <= period; i++){ int ii = i-1;
      if(i%2){
        if(i <= kk+1) map[ii] = 2*k + 1 + i;
        else if(i <= 3*kk+1) map[ii] = 2*k + 1 + i + PAR(kk+1);
        else if(i <= 2*k+3) map[ii] = 2*k - 1 + i;
        else if(i <= period - 3*kk + 1) map[ii] = 6*k + 6 - i;
        else if(i <= period - kk - 1) map[ii] = 6*k + 4 - i + PAR(kk);
        else map[ii] = 6*k + 4 - i;
      } else {
        if(i <= kk+1) map[ii] = 2*k + 2 - i;
        else if(i <= 3*kk+1) map[ii] = 2*k + 2 - i + PAR(kk+1);
        else if(i <= 2*k+3) map[ii] = 2*k + 4 - i;
        else if(i <= period - 3*kk + 1) map[ii] = i - 2*k - 3;
        else if(i <= period - kk - 1) map[ii] = i - 2*k - 1 + PAR(kk);
        else map[ii] = i - 2*k - 1;
      }
    }
  }
}
/*
This function summarizes the three previous ones and gives the
maximum entropy of a cycle of any period including the conjectured
ones.
*/
void maxcycle(unsigned int period, int *map){
  if (period<4)ErrorAborting("permutation length must be at least 4!!");
  if(period%2) { maxoodcycle(period, map); return; }
  if(period%4 == 0) { maxcycle4k(period, map); return; }
  conjecturemaxcycle4kp2(period, map); return;
}
