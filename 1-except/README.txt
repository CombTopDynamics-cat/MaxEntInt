********* maxentint1except.cpp **************
---------------------------------------------------------------------
Computes the maximum entropy of cyclic maximodal permutations with one exception
by the method of generating all maximodal cycles of a given length and checking
whether they have one exception. Then, computing the associated Markov Matrix
and spectral radius. During this process the absolute maximum is selected.
The output gives the entropy, the permutation in map form and the Markov matrix
with all row and colum sums together with the total sum.

It works for period smaller than 200 (in fact <= 17; then the number of cases is
to big to finish in reasonable time.

Usage: (3 <) period (<= 200) output-filename [entropy-filter-level]
---------------------------------------------------------------------

********* GenerateYZLists1Except.cpp **************
---------------------------------------------------------------------
Generates the Y and Z lists two obtain the maximum entropy cycle
in C^1_n (Section 4).

The output is the number of elements in Y and Z lists (stdout)
respectively, and the files containing the actual lists.
Their names are "freealphalist-1ex-Y.<odd-length>" and
"freealphalist-1ex-Z.<odd-length>"

The files are direct acces and are written in an "A-list format".
This format consists in initializing the file with the following information
at the beginning:
odd-length: unsigned char
the number <n> of elements of the list: unsigned long int
<n> permutations in the list: <n> x (odd-length x unsigned char)

Usage: odd-length
---------------------------------------------------------------------


********* maxent-1ex-crossprod.cpp **************
---------------------------------------------------------------------
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
---------------------------------------------------------------------
