'''
Created on Jul 15, 2013

@author: williamdemeo
'''
import unittest
from org.uacalc.alg.conlat import BasicPartition
from closure.closure import Closure

class Test(unittest.TestCase):

    def setUp(self):
        
        # Test Case 0 (first, very basic test):
        aa0 = BasicPartition("|0,1|2,3|4,5|")
        aa1 = BasicPartition("|0,2|1,4|3,5|")
        self.partsAA = aa0, aa1
        
        # Test Case 1 (another basic test):
        a0 = BasicPartition("|0|1,2,3,4|5,6|7,8|")
        a1 = BasicPartition("|0,1,2|3,4,5|6,7,8|")
        a2 = BasicPartition("|0,1,8|2,3|4,5|6,7|")
        a3 = BasicPartition("|0,1,2,3,4|5,6,7,8|")
        self.partsA = a0, a1, a2, a3
            
        # Test Case 2 (PJ11):
        H = BasicPartition("|0,1,6|2,7,35|3,22,36|4,8,17|5,9,18|19,26,56|10,20,57|11,21,58|12,45,59|13,46,60|16,25,39|51,66,96|32,70,84|14,23,37|15,24,38|34,48,62|28,43,80|40,49,77|41,50,78|33,47,61|65,74,95|27,42,79|54,90,99|29,44,81|86,91,103|30,68,82|31,69,83|87,92,104|64,73,94|75,88,105|52,67,97|100,102,107|55,71,85|76,101,106|53,89,98|63,72,93|")
        A = BasicPartition("|0,1,6,15,24,38,55,71,85|2,7,11,21,29,35,44,58,81|3,22,36,53,54,89,90,98,99|4,8,16,17,25,33,39,47,61|5,9,14,18,23,34,37,48,62|19,26,40,49,56,63,72,77,93|10,20,28,43,52,57,67,80,97|12,13,45,46,59,60,76,101,106|27,42,51,66,75,79,88,96,105|30,31,32,68,69,70,82,83,84|41,50,64,73,78,86,91,94,103|65,74,87,92,95,100,102,104,107|")
        B = BasicPartition("|0,1,4,5,6,8,9,14,15,16,17,18,23,24,25,33,34,37,38,39,47,48,55,61,62,71,85|2,7,10,11,20,21,27,28,29,35,42,43,44,51,52,57,58,66,67,75,79,80,81,88,96,97,105|3,12,13,22,30,31,32,36,45,46,53,54,59,60,68,69,70,76,82,83,84,89,90,98,99,101,106|19,26,40,41,49,50,56,63,64,65,72,73,74,77,78,86,87,91,92,93,94,95,100,102,103,104,107|")
        C = BasicPartition("|0,1,2,3,6,7,19,22,26,35,36,56|4,8,10,12,17,20,40,45,49,57,59,77|5,9,11,13,18,21,41,46,50,58,60,78|16,25,29,32,39,44,65,70,74,81,84,95|33,47,51,53,61,66,86,89,91,96,98,103|14,23,27,30,37,42,63,68,72,79,82,93|15,24,28,31,38,43,64,69,73,80,83,94|34,48,52,54,62,67,87,90,92,97,99,104|55,71,75,76,85,88,100,101,102,105,106,107|")
        K = BasicPartition("|0,4,14|1,9,25|2,51,52|3,13,32|5,15,33|6,61,62|7,43,88|22,45,68|8,24,48|26,73,102|10,11,75|19,86,87|12,31,54|27,28,29|16,34,55|35,57,79|36,83,106|17,18,85|20,44,66|21,42,67|56,78,95|46,69,89|37,38,39|70,90,101|23,47,71|49,74,91|50,72,92|40,41,100|63,64,65|30,53,76|58,80,96|59,84,98|60,82,99|77,94,104|81,97,105|93,103,107|")
        self.partsB = H,A,B,C,K

        # Test Case 3 (parallel sum of M3's):
        ad = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,18,19|8,9,20,21|10,11,22,23|30,31|32,33|34,35|")
        bd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,30,31|8,9,32,33|10,11,34,35|18,19|20,21|22,23|")
        cd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7|8,9|10,11|18,19,30,31|20,21,32,33|22,23,34,35|")
        da = BasicPartition("|0,2,4,6,8,10|1,3,7,9|5,11|12,14,16,18,20,22|13,15,19,21|17,23|24,26,28,30,32,34|25,27,31,33|29,35|")
        db = BasicPartition("|0,2,4,6,8,10|1,5,7,11|3,9|12,14,16,18,20,22|13,17,19,23|15,21|24,26,28,30,32,34|25,29,31,35|27,33|")
        dc = BasicPartition("|0,2,4,6,8,10|1,7|3,5,9,11|12,14,16,18,20,22|13,19|15,17,21,23|24,26,28,30,32,34|25,31|27,29,33,35|")
        self.partsC = ad,bd,cd,da,db,dc

    def tearDown(self):
        pass

    def test_sd_embedding(self):
        cl = Closure()
        
        fun_name = "sd_embedding()"
        print "\n===== Testing", fun_name, "====="
        # ----Test 0----  
        current_test_number = 0
        correct_ans = [[0,0], [0,1], [1,0], [1,2], [2,1], [2,2]]
        print "\n--- Test", current_test_number, "---"
        cl.set_partitions(self.partsAA)
        print "    cl.partitions = ", cl.partitions
        ans = cl.sd_embedding()
        print "    The subdirect embedding is:", ans
        self.assertEquals(ans, correct_ans, "Test "+str(current_test_number)+": " + fun_name + "seems broken")

        # ----Test 1----
        current_test_number+=1
        correct_ans = [[0,0,0,0], [1,0,0,0], [1,0,1,0], [1,1,1,0], [1,1,2,0], [2,1,2,1], [2,2,3,1], [3,2,3,1], [3,2,0,1]]
        print "\n--- Test", current_test_number, "---"
        cl.set_partitions(self.partsA)
        print "    cl.partitions = ", cl.partitions
        ans = cl.sd_embedding()
        print "    The subdirect embedding is:", ans
        self.assertEquals(ans, correct_ans, "Test "+str(current_test_number)+": sd_embedding seems broken")

    def test_optimal_sdf_subset(self):
        cl = Closure()
        fun_name = "optimal_sdf_subset()"
        print "\n===== Testing", fun_name, "====="

        # ----Test 0----  
        current_test_number = 0
        print "\n--- Test", current_test_number, "---"
        cl.set_partitions(self.partsAA)
        correct_ans = self.partsAA
        print "    cl.partitions = ", cl.partitions
        ans = cl.optimal_sdf_subset()
        print "    The optimal sdf subset is:", ans
        self.assertEquals(ans, correct_ans, "Test "+str(current_test_number)+": " + fun_name + "seems broken")

        # ----Test 1----  
        current_test_number+=1
        print "\n--- Test", current_test_number, "---"
        cl.set_partitions(self.partsA)
        #correct_ans = self.partsA
        print "    cl.partitions = ", cl.partitions
        ans = cl.optimal_sdf_subset()
        print "    The optimal sdf subset is:", ans
        #self.assertEquals(ans, correct_ans, "Test "+str(current_test_number)+": " + fun_name + "seems broken")


        # ----Test 2----  
        current_test_number+=1
        print "\n--- Test", current_test_number, "---"
        cl.set_partitions(self.partsB)
        #correct_ans = self.partsA
        print "    cl.partitions = ", cl.partitions
        ans = cl.optimal_sdf_subset()
        print "    The optimal sdf subset is:", ans
        #self.assertEquals(ans, correct_ans, "Test "+str(current_test_number)+": " + fun_name + "seems broken")


        # ----Test 3----  
        current_test_number+=1
        print "\n--- Test", current_test_number, "---"
        cl.set_partitions(self.partsC)
        #correct_ans = self.partsA
        print "    cl.partitions = ", cl.partitions
        ans = cl.optimal_sdf_subset()
        print "    The optimal sdf subset is:", ans
        #self.assertEquals(ans, correct_ans, "Test "+str(current_test_number)+": " + fun_name + "seems broken")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()