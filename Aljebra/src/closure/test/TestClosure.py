'''
Created on Jul 15, 2013

@author: williamdemeo
'''
import unittest
from closure.closure import Closure

class Test(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testAdd(self):
        cl = Closure()
        result = cl.add(opa=2,opb=3)
        self.assertEquals(result, 5, "add method broken")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()