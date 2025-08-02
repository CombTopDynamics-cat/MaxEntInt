/*
Computes the maximum entropy of the cross-product permutations
a <cross> b of the form descibed below WHICH ARE CYCLES (Section 3):

On input a file that must be in a valid A-list format must be
specified. Optionally, a segment of the list can be specified
(the default segment is the whole list) and a "report-each-%-of-work"
can also be specified (default is 1%).

Then the program generates all products a <cross> b with:
a varying in the chosen segment
b \ge a in the order of the list
AND KEEPS ONLY THOSE WHICH ARE MMAAXXIIMMOODDAALL CYCLES.

For these cycles the program computes the associated Markov Matrix
and then its spectral radius. During this process the absolute
maximumis selected. The output gives the entropy and the product cycle
giving it. The cycle is specified by giving the position occupied by
a and b in the alpha list. To know the actual elements a and b, the
cycle in map form and the Markov matrix with all row and colum sums
together with the total sum the program singleprodinfo.exe must be used.

When a full product cicle of the program is done (that is multiplying an
element a by all b's such that b \ge a) and a percentage of the task larger
or equal than "report-each-%-of-work" since the last time that a report has
been generated, then the program writes to the file the current maximum
spectral radius with all cycles having it. Of course this is a security
measure for long lasting programs to be aesily continued in case of non
expected interruption with a minimal loss.

On input the parameter entropy-filter-level is 1.0 by default. We only
look for cycles with entropy larger than this specified entropy-filter-level.
This is a device to increase the speed of the program in critical cases when
we know in advance that the entropy must be larger than an a priori value.

This program is essentially maxentcrossprod.cpp except for the following
CHANGES:
* Different (and more powerful and versatile) specification of parameters
* Report of work each certain percentage of work done. This is a security
  device for long lasting runs
* Checking that the products are maximodal. This is a special requirement
  for the Y and Z list (see Lemmas 4.5 and 4.9 and Remarks 4.6 and 4.10)

Usage: [options] list-inputfile outputfile [options]

  Options:
     -f=entropy-filter-level
     -i=ini-interval
     -e=end-interval
     -r=report-each-%-of-work
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // For the function access
#include <limits.h> // For ULONG_MAX
#include <math.h> // For the function floor
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/maxent.h"
#define TOL 1.e-8
char OptParam(const char *);

int main (int argc, char *argv[])
{
  int i; char ParamOK=1;
  int lif=0, of=0;
  unsigned long int ini=0L, end=0L, endlist, repnum, count=0L;
  register unsigned long int a, b;
  double rep=100.0, maxentini=1.0, maxent=1.0, topent, totaltask;
  FILE *fp;
  char Usage_msg[] = "[options] list-inputfile outputfile [options]\n\n\
  Options:\n\
     -f=entropy-filter-level\n\
     -i=ini-interval\n\
     -e=end-interval\n\
     -r=report-each-%-of-work";

  for(i=1; i< argc; i++){ long int aux;
    switch(OptParam(argv[i])){
      case 'f': maxent = maxentini = atof(argv[i]+3); break;
      case 'i': aux = atol(argv[i]+3); if(aux < 1) { ParamOK=0; break; } ini = aux - 1; break;
      case 'e': aux = atol(argv[i]+3); if(aux < 1) { ParamOK=0; break; } end = aux; break;
      case 'r': rep = atof(argv[i]+3); break;
      default : if(!lif) lif = i; else if(!of) of = i ; else ParamOK=0; break;
    }
  }
  if(!ParamOK || !lif || !of || maxent < 1.0 || rep <= 0.0 || rep > 100.0 ) Usage(argv[0], Usage_msg);
  if (!access (argv[of], F_OK)) ErrorAborting("output file already exists");

  ProductPermutation alphabeta(argv[lif]); endlist = alphabeta.SizeOfAlphaList();
  if(!end) end = endlist;
  if( end > endlist ) ErrorAborting("end-interval too large");
  if(end <= ini) Usage(argv[0], Usage_msg);

/* Given an interval [i,j) = {l: i <= l < j } of elements of the list N(i,j)
   denotes the number of operations associated to this interval.
   That is: N(i,j) = \sum_{l=i}^{j-1} (M - (l-1)) = (j-i)(2M+3-i-j)/2
   Observation: according to the internal notation: i=ini+1 and j=end+1.
   Thus, j-i = end - ini and 3-i-j = 1 - ini - end.
*/
  topent = rep/100.0; rep = 2.0*((double) endlist) + 1.0 - ((double) ini);
  totaltask = floor(((rep - end)*(end-ini))/2.0 + 0.5);
  topent = floor(totaltask*topent); repnum = (topent < ULONG_MAX) ? ((unsigned long int) topent) : ULONG_MAX;
  if(repnum < endlist) repnum = endlist;

  IntervalMarkovMatrix MarMat(alphabeta.Period());

  for(a=ini; a < end ; a++) { alphabeta.Readalpha(a);
    for (b=a; b < endlist ; b++){ alphabeta.Readbeta(b);
      if(!alphabeta.Cross() || !alphabeta.IsMaximodal()) continue;
      MarMat.Compute((BasePermutation) alphabeta);
      topent = MarMat.ComputeSpectralRadiusIfNotSmallerThan (maxent);
      if (topent < maxent - TOL) MarMat.StopIfErrorInPowerMethod (maxent,(BasePermutation) alphabeta);
      else if (topent > maxent + TOL) { maxent = topent; ResetInfoResults(a, b, 'x', topent); }
      else AddInfoToResults(a, b, 'x', topent);
    } b = endlist - a;
    if(count < repnum - b) count += b;
    else { count = 0L;
      fp = sfopen(argv[of], "w", "when opening output file");
      PrintInfoResults(fp, argv[lif], alphabeta.Period(), ini, end, maxentini);
      fprintf(fp, "PARTIAL RESULT: Alpha interval REALLY checked: %ld to %ld (%.2f\%)\n\n",
                  ini+1, a+1, floor(((rep - a)*(a-ini))/2.0 + 0.5)*100.0/totaltask );
      fclose(fp);
    }
  }

  fp = sfopen(argv[of], "w", "when opening output file");
  PrintInfoResults(fp, argv[lif], alphabeta.Period(), ini, end, maxentini);
  fclose(fp);
  return 0;
}

#include <ctype.h>
char OptParam(const char *parm){
 if(parm[0] != '-' || !isalpha(parm[1]) || parm[2] != '=') return 0;
 return parm[1];
}
