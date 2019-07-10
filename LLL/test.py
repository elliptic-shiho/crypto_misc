from __future__ import print_function, division, absolute_import
import unittest
import test_vector
import test_lattice

def main():
  suite = unittest.TestSuite()
  suite.addTests(unittest.makeSuite(test_vector.TestVector))
  suite.addTests(unittest.makeSuite(test_lattice.TestLattice))
  unittest.TextTestRunner().run(suite)

main()
