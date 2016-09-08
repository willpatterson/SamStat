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

COORD = (5, 15)
FULL_MATCH = (5, 15)
LOWER_MATCH = (10, 20)
UPPER_MATCH = (1, 10)


class TestRegion(unittest.TestCase):
    """Tests for Region Class"""
    #def setUpClass(cls):
    def test_coordinate_relations_full(self):
        lower, upper = Region.coordinate_relations(COORD, FULL_MATCH)
        self.assertTrue(lower)
        self.assertTrue(upper)


if __name__ == '__main__':
    test_classes = (TestRegionMap, TestRegion)
    test_suite = unittest.TestSuite()
    for test_class in test_classes:
        test_suite.addTest(test_class())

    runner = unittest.TextTestRunner()
    runner.run(test_suite)
