
This is flip, version 2.1.2.

flip - translates and reformats DNA sequences

This package is based on activities of the [OGMP](http://megasun.bch.umontreal.ca/ogmp/) (i.e, priori to 2002), and
becomes open source as part of [MFannot](http://megasun.bch.umontreal.ca/RNAweasel/).
Note: Original coder was B.F Lang and Nicolas Brossard

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

- Flip can handle sequences of at most 400,000 nucleotides. To accept longer sequences, try
   increasing MAX_SEQ_LG in flip.h..

- Flip will truncate lines that are longer than 5000 characters in the input file. For more
   characters per line, try increasing MAX_LINE_LG in flip.h.

NOTES
-----

It is possible that the flip.1 file is not up to date. 
