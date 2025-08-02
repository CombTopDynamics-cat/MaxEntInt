/*
This program is necessary because maxentcrossprod.exe and
maxentdotprod.exe do not display enough information on the
cycles of maximum entropy. They only display the maximum
spectral radius and the position of a and b in the list.
They do not display the actual permutation of those elements
the permutation of the resulting cycle and the Markov matrix.

On input a file that must be in a valid A-list format (Section 3)
must be specified, two valid indices a and b in the file, and one
type of product:
x for cross
. for dot

The program tests that both products of these two elements are cycles
and maximodal. In that case it performs the chosen product and displays
the actual permutations corresponding to a and b, the product cycle in
map form, the topological entropy and the Markov matrix with all row
and colum sums, together with the total sum.

Usage: filename alphaindex producttype(x,.) betaindex
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/maxent.h"
#define TOL 1.e-8

int main (int argc, char *argv[])
{
  unsigned long int a, b, endlist;
  double topent;
  int test;
  char message[255];

  // Testing parameters and open output file
  if( argc != 5 || atol(argv[2]) < 1 || atol(argv[4]) < 1 \
      || strlen(argv[3]) > 1 || ( argv[3][0] != 'x' && argv[3][0] != '.') )
           Usage(argv[0], "filename alphaindex producttype(x,.) betaindex");

  ProductPermutation alphabeta(argv[1]); endlist = alphabeta.SizeOfAlphaList();
  a = (unsigned long int) atol(argv[2])-1 ; b = (unsigned long int) atol(argv[4])-1 ;

  sprintf(message, "alphaindex > %lu (the number of elements in the list)", endlist);
  if(a >= endlist) ErrorAborting(message);
  if(b >= endlist){ strncpy(message, " betaindex ", 10); ErrorAborting(message); }

  alphabeta.Readalpha(a); alphabeta.Readbeta(b);
  if (argv[3][0] == 'x')  test = alphabeta.Cross(); else  test = alphabeta.Dot();
  if(!test) ErrorAborting("the product of these two permutations is not a cycle");
  if(!alphabeta.IsMaximodal()) ErrorAborting("the product of these two permutations is not maximodal");


  IntervalMarkovMatrix MarMat(alphabeta.Period());
  MarMat.Compute((BasePermutation) alphabeta);
  topent = MarMat.SpectralRadius ();

  alphabeta.PrintAll (stdout);
  fprintf(stdout, "Entropy = log %.16lf\n", topent);
  MarMat.Print (stdout);

  if (topent < -10) ErrorAborting("no convergence in Power Method");
  else return 0;
}
