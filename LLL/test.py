from __future__ import print_function, division, absolute_import
import unittest
import test_vector

def main():
  suite = unittest.TestSuite()
  suite.addTests(unittest.makeSuite(test_vector.TestVector))
  unittest.TextTestRunner().run(suite)

main()
