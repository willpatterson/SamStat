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


Gene = Region.Gene
COORDS = {'asedf':Gene((1,10),1,1),
          'ade':Gene((20, 30),1,1),
          '4d':Gene((40, 50),1,1),
          '34d':Gene((60, 70),1,1),
          'w3':Gene((80, 100),1,1)}

FULL_MATCH_BINARY = (20, 30)

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

    def test_coordinate_relations_lower(self):
        lower, upper = Region.coordinate_relations(COORD, LOWER_MATCH)
        self.assertTrue(lower)
        self.assertFalse(upper)

    def test_coordinate_relations_upper(self):
        lower, upper = Region.coordinate_relations(COORD, UPPER_MATCH)
        self.assertFalse(lower)
        self.assertTrue(upper)

    def test_binary_coordinate_match_full(self):
        sorted_coords = sorted(COORDS.values())
        match = Region.binary_coordinate_match(sorted_coords, FULL_MATCH_BINARY)
        self.assertTrue(match.lower)
        self.assertTrue(match.upper)
        self.assertEqual(FULL_MATCH_BINARY, match.value.location)
        self.assertEqual(1, match.index)




if __name__ == '__main__':
    test_classes = (TestRegionMap, TestRegion)
    test_suite = unittest.TestSuite()
    for test_class in test_classes:
        test_suite.addTest(test_class())

    runner = unittest.TextTestRunner()
    runner.run(test_suite)
