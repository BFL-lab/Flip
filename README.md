
This is flip, version 2.1.2.

flip - translates and reformats DNA sequences

This package is based on old work that was done at the [OGMP](http://megasun.bch.umontreal.ca/ogmp/),
we decide to put this work open source in order to put [MFannot](http://megasun.bch.umontreal.ca/RNAweasel/) open source.
For sure lot of other bioinformatics tools does the same thing.


PACKAGE CONTENTS
----------------

The package should contain the following files

  - README.md
  - flip.1 (file describing how to use flip in a Unix-like manpage format)
  - src/Makefile
  - src/code.h 
  - src/flip.h
  - src/flip.c


REQUIREMENTS
------------

- The Unix C compiler "gcc" or "cc".


HOW TO BUILD THE EXECUTABLE
---------------------------

  gcc -o flip flip.c   or
   cc -o flip flip.c


LIMITATIONS
-----------

- Flip can handle sequences of at most 400000 nucleotides. If this is
  too much of a restriction, try increasing MAX_SEQ_LG in flip.h.

- Flip will truncate lines that are longer than 5000 characters in the
  input file supplied on the command line. If it's a limitation for you, try
  increasing MAX_LINE_LG in flip.h.


NOTES
-----

It is possible that the flip.1 file is not up to date. 
