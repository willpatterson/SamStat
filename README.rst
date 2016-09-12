*******
SamStat
*******

SamStat is a custom toolkit used in a bioinformatics pipeline.

QNAME Statistics:
-----------------

``qstat`` preforms the following operations on QNAMES from a SAM file:

1. Number of lines per QNAME
2. Number of unquie RNAMEs per QNAME
3. Range of RNAMEs per QNAME
4. Find feature classification (exonic, intronic, intergenetic, combination) of alignments using a GFF3 file

True Direction (``truedir``):
-----------------------------

``truedir`` deduces the true direction of each sequence per-RNAME for every QNAME in a SAM file using a GFF3 file

Install:
--------

**requires python3**
Developled using Python 3.4.1

::

  git clone https://github.com/willpatterson/SamStat.git
  cd SamStat
  python setup.py install


Usage:
------

::

  #Help
  samstat --help 
  samstat -h

  #qstat usage
  samstat qstat <SAM_filepath> <GFF3_filepath> <Out_filepath>

  #truedir usage
  samstat truedir <SAM_filepath> <GFF3_filepath> <Out_filepath>

