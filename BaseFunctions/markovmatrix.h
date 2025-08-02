/*********************************************************
  Class MarkovMatrix (general) and
  class IntervalMarkovMatrix declarations

  Used by maxent-1ex-crossprod.cpp, mapentint.cpp,
  maxentint.cpp, maxentallbounds.cpp, singleprodinfo.cpp,
  maxentcrossprod.cpp, maxentdotprod.cpp
**********************************************************/

class BasePermutation;

class MarkovMatrix{
 unsigned int size2;
 double *v, *av;
 void InitializePowerMethod();
 unsigned int PowerMethodIterate();
 double PowerMethod(double);
 void DisplaceMatrix(int);
 unsigned int UpperBoundForSpectralRadius();
 void InitializeMemory(unsigned int);
protected:
  unsigned char *mm;
  unsigned int size;

public:
  MarkovMatrix (unsigned int order){ InitializeMemory(order); };
  void Print (FILE *);
  void PrintAll(FILE *, double, BasePermutation);
  void ResizeMarkovMatrix(unsigned int order){free(mm); free(v); free(av); InitializeMemory(order); };

  double SpectralRadius ();
  double ComputeSpectralRadiusIfNotSmallerThan (double);
  void StopIfErrorInPowerMethod (double, BasePermutation);
};

class IntervalMarkovMatrix:public MarkovMatrix
{
  void FillImageRowofIntervalMarkovMatrix (unsigned int, unsigned int, unsigned int);
public:
  IntervalMarkovMatrix (unsigned int order):MarkovMatrix::MarkovMatrix (order) {};
  void MapCompute (unsigned int *);
  void Compute(BasePermutation);
};
