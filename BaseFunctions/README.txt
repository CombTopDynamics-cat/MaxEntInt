********* maxent.h --- maxent.cpp **************
---------------------------------------------------------------------
General auxiliary functions declarations
Used in maxent.cpp, permutations.cpp and markovmatrix.cpp

Contains:
=========
* A function to decode integer parameters in a predefined interval
* Error message functions
Routines to store and print results for maxentprods
---------------------------------------------------------------------

********* permutations.h --- permutations.cpp **************
---------------------------------------------------------------------
These are the files that implement the basic structures for
generating permutations, lists, dealing with products.... They are
the kernel of the whole set of programs. A big effort has been made
in documenting these files (mainly permutations.cpp)

    class BasePermutation (base class)
    class CiclicPermutation
    class MaximodalCiclicIntervalPermutation
    class AlphaPermutation
    class RestrictedAlphaPermutation
    class ProductPermutation
    class AlphaLists1Exc

    All of them derived from class BasePermutation (base class)

Used by GenerateYZLists1Except.cpp, maxent-1ex-crossprod.cpp,
  mapentint.cpp, maxentint.cpp, maxentallbounds.cpp,
  CountJungreisCycles.cpp, singleprodinfo.cpp, maxentcrossprod.cpp,
  maxentdotprod.cpp, ReadAList.cpp, GenerateFreeAList.cpp,
  GenerateRestrictedAlphaList.cpp
---------------------------------------------------------------------

********* markovmatrix.h --- markovmatrix.cpp **************
---------------------------------------------------------------------
Class MarkovMatrix (general) and
class IntervalMarkovMatrix declarations

Used by maxent-1ex-crossprod.cpp, mapentint.cpp,
maxentint.cpp, maxentallbounds.cpp, singleprodinfo.cpp,
maxentcrossprod.cpp, maxentdotprod.cpp

Contains:
=========
* IntervalMarkovMatrixvFunctions
* MarkovMatrix Functions
---------------------------------------------------------------------

********* longpermutations.h --- longpermutations.cpp **************
---------------------------------------------------------------------
Class LongBasePermutation declarations
Used by mapentint.cpp and maxentallbounds.cpp
This is another instance of class BasePermutations
that allows lengths up to LONGMAXPERM (=1000).

For comments we refer to the class BasePermutations
---------------------------------------------------------------------
