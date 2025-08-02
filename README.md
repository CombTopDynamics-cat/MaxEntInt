# MaxEntInt
This repository contains the programs used in the paper

"A lower bound for the maximum topological entropy of 4k+2-cycles"
                             by
               Lluis Alseda, David Juher and Deborah King


The directories:
     1-except,
     general,
     general-products and
     restricted-products
contain all programs classified by subject.

The folder 'BaseFunctions' contains functions common to all programs.

Each folder contains file named 'README.txt' with a brief description
of the contents of the folder and its use.

Each of the first 4 folders (namely: 1-except, general,
general-products and restricted-products) contains a standard
'Makefile' that allows the easy compilation and linking of the
programs.

To obtain all executables in a folder go to the corresponding
folder and at the shell prompt execute the command:
    'make all'

As usual the command 'make clean' will clean the created object files.

The executables have been created without any problem in different
Linux (Debian based --- as different versions of Debian and Ubuntu)
operating systems. We do not know about other Linux Systems and
Windows.
