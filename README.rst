*******
SamStat
*******

SamStat is a custom script used in a bioinformatics pipeline.

SamStat preforms the following operations on corrispoding .sam file and .gff 
files:

1. Number of lines per QNAME
2. Number of forward/reverse per QNAME
3. Number of unquie RNAMEs per QNAME
4. Range of RNAMEs per QNAME
5. Link attribute type of sequences using Gnomon GFF files

Install:
--------

**requires python3**

::

  git clone https://github.com/willpatterson/SamStat.git
  cd SamStat
  python setup.py install

Usage:
------

::

  samstat <SAM_filepath> <GFF_filepath> <Out_filepath>
  samstat --help
