********* maxentint.cpp **************
---------------------------------------------------------------------
Computes the maximum entropy of Jungreis cycles by the
method of generating all maximodal cyclic permutations of a given
length and then checking whether they have zero or one exceptions.
Then, computing the associated Markov Matrix and then its spectral
radius. During this process the absolute maximum is selected.

The output gives the entropy, the cycle in map form and the Markov
matrix with all row and colum sums together with the total sum.

It works for period smaller than 14 (17?). Then the number of cases is to
big to finish in reasonable time.

Usage: (3 <) period (<= 200) output-filename [entropy-filter-level]
---------------------------------------------------------------------

********* mapentint.cpp **************
---------------------------------------------------------------------
Computes the entropy of a cycle specified by a permutation in
map form. The output gives the entropy, the cycle
in map form and the Markov matrix with all row and colum sums
together with the total sum.

It works for period smaller than 1000.

Usage: <map vector i_1 i_2 i_3 ... i_n (2 < length)>
---------------------------------------------------------------------

********* maxentallbounds.cpp **************
---------------------------------------------------------------------
For every period between 4 and a maximum period entered as a parameter
three numbers are computed:
* The Misiurewicz-Nitecki bound
* The maximum entropy of permutations of period n
* The maximum entropy of n-cycles (including conjectured families)

These numbers are printed in table form with four columns.
The first column crresponds to the period.

It works for period smaller than 1000.

Usage: maximum-length(>4)
---------------------------------------------------------------------

********* CountJungreisCycles.cpp **************
---------------------------------------------------------------------
Computes all non-conjugate cycles of a given period.
For each maximal cycle it is specified this property
and the number of exceptions

It works for period smaller than 200.

Usage: (3 <) period (<= 200)
---------------------------------------------------------------------
