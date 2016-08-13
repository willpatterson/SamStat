*******
SamStat
*******

SamStat is a custom toolkit used in a bioinformatics pipeline.


SamStat preforms the following operations on qnames from a SAM file and corrispoding GFF3 file  
files:

1. Number of lines per QNAME
2. Number of unquie RNAMEs per QNAME
3. Range of RNAMEs per QNAME
4. Link attribute type of sequences using Gnomon GFF3 files

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

  samstat <SAM_filepath> <GFF3_filepath> <Out_filepath>
  samstat --help 
  samstat -h
