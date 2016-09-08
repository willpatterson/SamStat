import unittest
import os
import warnings

from samstat.maps import Region
from samstat.maps import RegionMap

IN_GFF = '/disk/bioscratch/Will/Drop_Box/GCF_001266775.1_Austrofundulus_limnaeus-1.0_genomic_andMITO.gff'
class TestRegionMap(unittest.TestCase):
    """Test Class for RegionMap"""
    def setUpClass(cls):
        """Reads in GFF file and creates Region Map object"""
        cls.rm = RegionMap(IN_GFF)
        cls.rm_list = list(rm.rmap.items())


if __name__ == '__main__':
    test_classes = (TestRegionMap)
    test_suite = unittest.TestSuite()
    for test_class in test_classes:
        test_suite.addTest(test_class())

    runner = unittest.TextTestRunner()
    runner.run(test_suite)
