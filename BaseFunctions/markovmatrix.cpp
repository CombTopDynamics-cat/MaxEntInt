/*********************************************************
  Class MarkovMatrix (general) and
  class IntervalMarkovMatrix code of functions

  Used by maxent-1ex-crossprod.cpp, mapentint.cpp,
  maxentint.cpp, maxentallbounds.cpp, singleprodinfo.cpp,
  maxentcrossprod.cpp, maxentdotprod.cpp
**********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/maxent.h"
#define MIN(x,y) ((x) < (y) ? (x) : (y))

/*******************************************/
/* IntervalMarkovMatrix   Functions        */
/*******************************************/
/*
   Computes the Markov Matrix from a permutation. Uses the base function
   FillImageRowofIntervalMarkovMatrix
*/
void IntervalMarkovMatrix::Compute(BasePermutation perm){ register unsigned int r;
  for (r = 0; r < size; r++)
    FillImageRowofIntervalMarkovMatrix (r, (unsigned int) perm.perm[r], (unsigned int) perm.perm[r + 1]);
}
/*
   Computes the Markov Matrix from a permutation written in map from. Uses the base function
   FillImageRowofIntervalMarkovMatrix
*/
void IntervalMarkovMatrix::MapCompute (unsigned int *map){ register unsigned int r;
       for (r = 0; r < size; r++)  FillImageRowofIntervalMarkovMatrix (r, map[r], map[r + 1]);
}

/* Given the images of the two endpoints of a basic (Markov) interval and assuming that the map is an
   interval map it computes the corresponding row of the Markov Matrix of the connect-the-dots with
   the specified images.
*/

void IntervalMarkovMatrix::FillImageRowofIntervalMarkovMatrix (unsigned int row, unsigned int fa, unsigned int fb)
{ register unsigned int aux; register unsigned char *msrow;

  if (fa > fb) { aux = fa; fa = fb; fb = aux; }
  for (msrow = mm + row * size, aux = 0; aux < size; msrow++, aux++) *msrow = 0;
  for (msrow = mm + row * size + fa, aux = fa; aux < fb; msrow++, aux++) *msrow = 1;
}

/*******************************************/
/* MarkovMatrix   Functions                */
/*******************************************/

void MarkovMatrix::InitializeMemory (unsigned int order)
{
  size = order - 1;
  size2 = size * size;

  mm = (unsigned char *) malloc (sizeof (unsigned char) * size2);
  v = (double *) malloc (sizeof (double) * size);
  av = (double *) malloc (sizeof (double) * size);
  if (!mm || !v || !av) ErrorAborting("memory problems initializing the Markov Matrix");
  return;
}

void MarkovMatrix::StopIfErrorInPowerMethod (double rho , BasePermutation perm)
{  if (rho >= -10) return;
   fprintf (stderr, "\n"); ErrorMessage("no convergence in Power Method:\n\n");
   perm.Print (stderr); Print(stderr);
   fprintf (stderr, "Aborting...\n\n"); exit(1);
}

void MarkovMatrix::PrintAll(FILE *fp, double rho , BasePermutation perm)
{ perm.Print(fp); fprintf (fp, "Entropy = log %.16lf\n", rho); Print(fp); }

void MarkovMatrix::InitializePowerMethod ()
{ unsigned register int i;
  for (i = 0; i < size; i++)  av[i] = 1.0;
}

unsigned int MarkovMatrix::PowerMethodIterate ()
{
  unsigned char *matrix = mm;
  register unsigned int i, j;
  double ep, vmax = -1.0;
  unsigned int imax = 0;

  for (i = 0; i < size; i++)
    {
      ep = 0.0;
      for (j = 0; j < size; j++)
	{
	  ep += (*matrix) * av[j];
	  matrix++;
	}
      v[i] = ep;
      if (ep > vmax)
	{
	  vmax = ep;
	  imax = i;
	}
    }

  for (i = 0; i < size; i++)
    av[i] = v[i] / vmax;
  return imax;
}

#define MaxIterRound 10000000
#define Transit 10
double MarkovMatrix::PowerMethod (double prec)
{
  unsigned int imax, imaxv;
  long unsigned int iter;
  double ent, entv = -1.0;

  for (iter = 0; iter < Transit; iter++) imaxv = PowerMethodIterate ();

  for (iter = 0; iter < MaxIterRound; iter++)
    {
      imax = PowerMethodIterate ();
      ent = v[imaxv];

      if (fabs (ent - entv) < prec)
	return ent;

      imaxv = imax;
      entv = ent;
    }
  ErrorMessage("convergence problems in Power Method. Abandoning...\n");
  return -11.0;
}

/*
   The next two functions implement a displacement trick to assure convergence
   when the -\rho is also an eigenvalue.
*/
void MarkovMatrix::DisplaceMatrix (int disp)
{
  unsigned char *matrix = mm;
  unsigned int i;

  for (i = 0; i < size; i++, matrix += size + 1)
    *matrix += disp;
}

#define PrecRho 1.e-12
double MarkovMatrix::SpectralRadius ()
{
  double rho;
  DisplaceMatrix (1);
  InitializePowerMethod ();
  rho = PowerMethod (PrecRho) - 1;
  DisplaceMatrix (-1);
  return rho;
}

/*
The next function implements the efficient computation of the spectral radius
of matrices in a maximization process. As a parameter it is entered the current
maximum level of spectral radius achieved.

* First it is computed the upper bound for the spectral radius based in the 1
and the infinity-norms. If this bound is found smaller than the current value
then the matrix can be discarded (and a fail value of -1 is returned to the
calling program).

* If this test fails it is performed a (quick) Power Method with lower
precision If this rough estimated value of the spectral radius is smaller
than "current - SecurityFactor" then we also can safely discard this matrix
(and a fail value of -1 is returned to the calling program).

* otherwise we can continue the Power Method from the data obtained in the
first rough step to perform a second step to compute the spectral radius with
high precision.Then this value is delivered to the calling program.

Of course this procedure is much slower than just performing the standard
Power Method. However,  when doing massive computations that there are a
lot of cases to examine after the "current" value has grown a lot or when we
know an "a priori" lower bound or the spectral radius and we enter it to the
program this strategy has revealed to be extremely efficient compared with the
straightforward use of the Power Method and choosing the largest spectral radius
at each analyzed permutation.
*/

#define ToleranceEstimate 1.e-04
// define SecurityFactor 2500*ToleranceEstimate
#define SecurityFactor 0.25
double MarkovMatrix::ComputeSpectralRadiusIfNotSmallerThan (double current)
{
  double rho;

  if (UpperBoundForSpectralRadius () < current) return -1.0;

  DisplaceMatrix (1);
  InitializePowerMethod ();
  rho = PowerMethod (ToleranceEstimate) - 1;
  if (rho < current - SecurityFactor) {
      if (rho > -10) return -1.0;
      ErrorAborting("unexpected error of convergence in Power Method");
  }
  rho = PowerMethod (PrecRho) - 1;
  DisplaceMatrix (-1);
  return rho;
}

unsigned int RowSums[1024];
void MarkovMatrix::Print (FILE * fp)
{
  unsigned register int i, j, s=0;
  unsigned int len, maxlen = 0;
  char aux[1025];

  for (i = 0; i < size; i++)
    RowSums[i] = 0;

  fprintf (fp, "\nAssociated Markov Matrix:\n");

  for (i = 0, s = 0; i < size2; i += size)
    { len = 0;
      for (j = 0; j < size; j++)
	{
	  unsigned char m;
	  m = mm[i + j];
	  s += m;
	  RowSums[j] += m;
	  len += fprintf (fp, "%d ", m);
	}
      len += fprintf (fp, "= %d\n", s);
      s = 0;
      if(len > maxlen) maxlen = len;
    }

  maxlen = MIN (maxlen, 1024);
  aux[maxlen] = 0;
  memset(aux, '-', maxlen); fprintf (fp, "%s\n", aux);

  s = 0;
  for (j = 0; j < size; j++)
    {
      s += RowSums[j];
      fprintf (fp, "%d ", RowSums[j]);
    }
  memset (aux, '=', maxlen); fprintf (fp, "= %d\n%s\n\n", s, aux);
}


/*
It computes the minimum of the norms of the Markov Matrix with the 1-norm
and with the infinity-norm. This is an upper bound for the spectral radius
of the Markov Matrix.
*/

unsigned int MarkovMatrix::UpperBoundForSpectralRadius ()
{
  unsigned int ub1 = 0, ub2 = 0;
  unsigned register int i, j, s;

  for (i = 0; i < size; i++)
    RowSums[i] = 0;

  for (i = 0, s = 0; i < size2; i += size)
    {
      for (j = 0; j < size; j++)
	{
	  unsigned char m;
	  m = mm[i + j];
	  s += m;
	  RowSums[j] += m;
	}
      if (ub1 < s)
	ub1 = s;
    }

  for (j = 0; j < size; j++)
    if (ub2 < RowSums[j])
      ub2 = RowSums[j];

  return MIN (ub1, ub2);
}
