'''
Created on Jul 15, 2013

@author: williamdemeo
'''
import unittest
import inspect  # for getting name of current function
from org.uacalc.alg.conlat import BasicPartition
from closure.closure import Closure
# from json.tests.test_encode_basestring_ascii import CASES
from compiler.transformer import asList

class Test(unittest.TestCase):

    @staticmethod
    def compare_orderable(s, t):
        return sorted(s)==sorted(t)
    
    @staticmethod
    def compare(s, t):
        t = list(t)   # make a mutable copy
        try:
            for elem in s:
                t.remove(elem)
        except ValueError:
            return False
        return not t
    
    def setUp(self):
        
        self.parts = ()
        
        # Test Case 0 (first, very basic test):
        aa0 = BasicPartition("|0,1|2,3|4,5|")
        aa1 = BasicPartition("|0,2|1,4|3,5|")
        self.parts = self.parts+ ([aa0, aa1],)
        
        # Test Case 1 (another basic test):
        a0 = BasicPartition("|0|1,2,3|4,5,6|7,8|")
        a1 = BasicPartition("|0,1,2|3,4,5|6,7,8|")
        a2 = BasicPartition("|0,1,8|2,3|4,5|6,7|")
        a3 = BasicPartition("|0,1,2,3,4|5,6,7,8|")
        self.parts = self.parts + ([a0, a1, a2, a3],)
            
        # Test Case 2 (PJ11):
        H = BasicPartition("|0,1,6|2,7,35|3,22,36|4,8,17|5,9,18|19,26,56|10,20,57|11,21,58|12,45,59|13,46,60|16,25,39|51,66,96|32,70,84|14,23,37|15,24,38|34,48,62|28,43,80|40,49,77|41,50,78|33,47,61|65,74,95|27,42,79|54,90,99|29,44,81|86,91,103|30,68,82|31,69,83|87,92,104|64,73,94|75,88,105|52,67,97|100,102,107|55,71,85|76,101,106|53,89,98|63,72,93|")
        A = BasicPartition("|0,1,6,15,24,38,55,71,85|2,7,11,21,29,35,44,58,81|3,22,36,53,54,89,90,98,99|4,8,16,17,25,33,39,47,61|5,9,14,18,23,34,37,48,62|19,26,40,49,56,63,72,77,93|10,20,28,43,52,57,67,80,97|12,13,45,46,59,60,76,101,106|27,42,51,66,75,79,88,96,105|30,31,32,68,69,70,82,83,84|41,50,64,73,78,86,91,94,103|65,74,87,92,95,100,102,104,107|")
        B = BasicPartition("|0,1,4,5,6,8,9,14,15,16,17,18,23,24,25,33,34,37,38,39,47,48,55,61,62,71,85|2,7,10,11,20,21,27,28,29,35,42,43,44,51,52,57,58,66,67,75,79,80,81,88,96,97,105|3,12,13,22,30,31,32,36,45,46,53,54,59,60,68,69,70,76,82,83,84,89,90,98,99,101,106|19,26,40,41,49,50,56,63,64,65,72,73,74,77,78,86,87,91,92,93,94,95,100,102,103,104,107|")
        C = BasicPartition("|0,1,2,3,6,7,19,22,26,35,36,56|4,8,10,12,17,20,40,45,49,57,59,77|5,9,11,13,18,21,41,46,50,58,60,78|16,25,29,32,39,44,65,70,74,81,84,95|33,47,51,53,61,66,86,89,91,96,98,103|14,23,27,30,37,42,63,68,72,79,82,93|15,24,28,31,38,43,64,69,73,80,83,94|34,48,52,54,62,67,87,90,92,97,99,104|55,71,75,76,85,88,100,101,102,105,106,107|")
        K = BasicPartition("|0,4,14|1,9,25|2,51,52|3,13,32|5,15,33|6,61,62|7,43,88|22,45,68|8,24,48|26,73,102|10,11,75|19,86,87|12,31,54|27,28,29|16,34,55|35,57,79|36,83,106|17,18,85|20,44,66|21,42,67|56,78,95|46,69,89|37,38,39|70,90,101|23,47,71|49,74,91|50,72,92|40,41,100|63,64,65|30,53,76|58,80,96|59,84,98|60,82,99|77,94,104|81,97,105|93,103,107|")
        self.parts = self.parts + ([H,A,B,C,K],)

        # Test Case 3 (parallel sum of M3's):
        ad = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,18,19|8,9,20,21|10,11,22,23|30,31|32,33|34,35|")
        bd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,30,31|8,9,32,33|10,11,34,35|18,19|20,21|22,23|")
        cd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7|8,9|10,11|18,19,30,31|20,21,32,33|22,23,34,35|")
        da = BasicPartition("|0,2,4,6,8,10|1,3,7,9|5,11|12,14,16,18,20,22|13,15,19,21|17,23|24,26,28,30,32,34|25,27,31,33|29,35|")
        db = BasicPartition("|0,2,4,6,8,10|1,5,7,11|3,9|12,14,16,18,20,22|13,17,19,23|15,21|24,26,28,30,32,34|25,29,31,35|27,33|")
        dc = BasicPartition("|0,2,4,6,8,10|1,7|3,5,9,11|12,14,16,18,20,22|13,19|15,17,21,23|24,26,28,30,32,34|25,31|27,29,33,35|")
        self.parts = self.parts + ([ad,bd,cd,da,db,dc],)

        # Test Case 4 (another basic test):
        x0 = BasicPartition("|0,1,2|3,4,5|")
        x1 = BasicPartition("|0,3|1,4|2,5|")
        x2 = BasicPartition("|0,4|1,5|2,3|")
        self.parts = self.parts + ([x0, x1, x2],)

        # Test Case 5 (a closed M_4 in Eq(6)):
        x0 = BasicPartition("|0,3,5|1,2,4|")
        x1 = BasicPartition("|0,4|1,5|2,3|")
        x2 = BasicPartition("|0,1|2,5|3,4|")
        x3 = BasicPartition("|0,2|1,3|4,5|")
        self.parts = self.parts + ([x0, x1, x2, x3],)
        
        # Test Case 6 (a possibly closed L_7 in Eq(30)):
        K = BasicPartition("|0,1,2,3,4|5,6,7,8,9|10,11,12,13,14|15,16,17,18,19|20,21,22,23,24|25,26,27,28,29|")
        M1 = BasicPartition("|0,5,12,19,24,27|1,6,13,17,22,28|2,7,14,18,23,29|3,8,10,15,21,26|4,9,11,16,20,25|")
        M2 = BasicPartition("|0,5,10,15,20,25|1,6,11,16,21,26|2,7,12,17,22,27|3,8,13,18,23,28|4,9,14,19,24,29|")
        J1 = BasicPartition("|0,10,20|1,16,26|2,12,17|3,13,18|4,14,24|5,15,25|6,11,21|7,22,27|8,23,28|9,19,29|")
        self.parts = self.parts + ([K, M1, M2],)
            
        # Test Case 7 (M4 in Eq(16)):
        b0 = BasicPartition("|0,1,2,3|4,5,6,7|8,9,10,11|12,13,14,15|")
        b1 = BasicPartition("|0,4,8,12|1,5,9,13|2,6,10,14|3,7,11,15|")
        b2 = BasicPartition("|0,5,10,15|1,4,11,14|2,7,8,13|3,6,9,12|")
        b3 = BasicPartition("|0,7,9,14|1,6,8,15|2,5,11,12|2,4,10,13|")
        self.parts = self.parts + ([b0, b1, b2, b3],)
            
        
   
    def tearDown(self):
        pass


    def test_decomp_size(self):
        fun_name = "decomp_size()"

        # Which examples to test:
        #test_cases = [0,1,3]
        test_cases = []  # (run no tests)

        correct_answer = {0: (3**3) * (3**3), 1: (4**4) * (3**3) *(4**4) * (2**2), 3: (9**9)**6 }
        # NB: the correct_answers are the values we expect to be returned by decomp_size 
        # __for the given partitions__, not necessarily for the optimal set of partitions 
        # that will be found later by the optimal_sdf_subset() function.

        if len(test_cases)>0:
            print "\n===== Testing", fun_name, "====="
            for case_number in test_cases:
                print "\n--- Test", case_number, "---"
                print "    partitions = ", self.parts[case_number]
                ans = Closure.decomp_size(self.parts[case_number])
                print "    decomp size for this set of partitions:", ans
                self.assertEquals(ans, correct_answer[case_number], "Test "+str(case_number)+": " + fun_name + "seems broken")

        
    def test_compute_sd_embedding(self):
        fun_name = "compute_sd_embedding()"

        # Which examples to test:
        test_cases = []  # (run no tests)
        
        correct_answer = {0: [[0,0], [0,1], [1,0], [1,2], [2,1], [2,2]], 1: [[0,0,0,0], [1,0,0,0], [1,0,1,0], [1,1,1,0], [2,1,2,0], [2,1,2,1], [2,2,3,1], [3,2,3,1], [3,2,0,1]]}
        
        if len(test_cases)>0:
            print "\n===== Testing", fun_name, "====="
            for case_number in test_cases:
                print "\n--- Test", case_number, "---"
                cl = Closure(self.parts[case_number])
                print "    cl.partitions = ", cl.partitions
                ans = cl.compute_sd_embedding()
                print "    The subdirect embedding is:", ans
                print "    The correct answer is:", correct_answer[case_number]
                self.assertEquals(ans, correct_answer[case_number], "Test "+str(case_number)+": " + fun_name + "seems broken")


    def test_compute_optimal_sdf_subset(self):
        # removing this for now, since it might be causing memory issues (though unlikely)
        test_fun_name = inspect.stack()[0][3]  # name of current function
        fun_name = test_fun_name[5:]  # name of function to be tested
        #fun_name = "compute_optimal_sdf_subset()"

        # Which examples to test:
        test_cases = []  
        #test_cases = [1,3]

        p30 = self.parts[3][0].join(self.parts[3][1])
        p31 = self.parts[3][3].join(self.parts[3][4])
        correct_answer = {0: self.parts[0],
                          1: self.parts[1], 
                          3: (p30, p31), 
                          4: (self.parts[4][0], self.parts[4][1]),
                          5: (self.parts[5][0], self.parts[5][1])}

        if len(test_cases)>0:
            print "\n===== Testing", fun_name, "====="

            for case_number in test_cases:
                print "\n--- Test", case_number, "---"
                cl = Closure(self.parts[case_number])
                print "    cl.partitions = ", cl.partitions
                ans = cl.compute_optimal_sdf_subset()
                print "    The optimal sdf subset is:", ans
                print "        ==> decomposition size:", Closure.decomp_size(ans)
                print "    The correct answer is:", correct_answer[case_number]
                print "        ==> decomposition size:", Closure.decomp_size(correct_answer[case_number])
                self.assertEquals(sorted(ans), sorted(correct_answer[case_number]), "Test "+str(case_number)+": " + fun_name + "seems broken")


    def test_compute_sd_Fix(self):
        test_fun_name = inspect.stack()[0][3]  # name of current function
        fun_name = test_fun_name[5:]  # name of function to be tested
        #fun_name = "compute_sd_Fix"

        # Which examples to test:
        #test_cases = [3]  # 3 is the parallel sum of M_3. universe size: 36 (it will take a while to finish)
        test_cases = [] # 
        
        if len(test_cases)>0:
            print "\n===== Testing", fun_name, "====="

            for case_number in test_cases:
                print "\n--- Test", case_number, "---"
                partitions = self.parts[case_number]
                #correct_ans = BasicPartition.unaryPolymorphisms(partitions, None)
                correct_ans = []
                if not case_number==3:  # case 3 takes too long to compute the old way
                    A = BasicPartition.unaryPolymorphismsAlgebra(partitions)
                    ops = A.operations()
                    for op in ops:
                        correct_ans.append(asList(op.getTable()))

                cl = Closure(partitions)
                print "    cl.partitions = ", cl.partitions
                FF = cl.compute_sd_Fix([], [])
                print "    There are ", len(FF), " unary polymorphisms",
                if not case_number == 3:
                    print FF
                    print "CORRECT ANS:"
                    print "    There are ", len(correct_ans), "unary polymorphisms:", correct_ans
                    self.assertEquals(sorted(FF), sorted(correct_ans), "Test "+str(case_number)+": " + fun_name + "seems broken")
                else:
                    print "\n"

    def test_algebra_from_unary_polymorphisms_filebased(self):
        test_fun_name = inspect.stack()[0][3]  # name of current function
        fun_name = test_fun_name[5:]  # name of function to be tested
        #fun_name = "basic_algebra_from_unary_polymorphisms()"

        # Which examples to test:
        #test_cases = [0,1,4,5]  # these all passed  
        #test_cases = [3]  # 3 is the parallel sum of M_3. universe size: 36 (it will take a while to finish)
        test_cases = [7]
        # Case 3 now works.  computed the closure on this 36 element algebra in just 5549 seconds.  
        #test_cases = [5]  # a closed M4
        #test_cases = [] # (run no tests)
        
        if len(test_cases)>0:
            print "\n===== Testing", fun_name, "====="

            for case_number in test_cases:
                print "\n--- Test", case_number, "---"

                partitions = self.parts[case_number]
                cl = Closure(partitions)
                print "    cl.partitions = ", cl.partitions

                FF = cl.compute_sd_Fix([], [])
                #A = cl.algebra_from_unary_polymorphisms_filebased(FF, "M3ParallelSum.ua")
                #A = cl.algebra_from_unary_polymorphisms_filebased(FF, "Elusive7.ua")
                A = cl.algebra_from_unary_polymorphisms_filebased(FF, "M4over.ua")
                print "Created algebra A with universe: ", A.universe()
                # compute congruence lattices, check they are equal to each other and to original set of partitions
                print "|ConA| = ", len(A.con().universe())
                #print "ConA = ", A.con().universe()


                if not case_number==7:  # case 6 takes too long to compute the old way
                    correct_ans = BasicPartition.unaryPolymorphismsAlgebra(partitions)
                    print "Correct algebra B has universe: ", correct_ans.universe()
                    print "|ConB| = ", len(correct_ans.con().universe())
                    print "ConB = ", correct_ans.con().universe()

                    self.assertEquals(A.universe(), correct_ans.universe(), "Test "+str(case_number)+": " + fun_name + "seems broken")
                    self.assertEquals(A.con().universe(), correct_ans.con().universe(), "Test "+str(case_number)+": " + fun_name + "seems broken")


    @staticmethod
    def unique_items(L):
        ans = []
        for p in L:
            if not sorted(p) in ans:
                ans.append(sorted(p))
        return ans

    def test_findMn(self):
        test_fun_name = inspect.stack()[0][3]  # name of current function
        fun_name = test_fun_name[5:]  # name of function to be tested

        # Which examples to test:
        test_cases = []
#         test_cases = {0:[3,2], # M_2 in Eq(3)
#                       1:[3,3], # M_3 in Eq(3)
#                       2:[4,3]} # M_3 in Eq(4)
        correct_ans = {0:3, # there are 3 M_2's in Eq(3)
                       1:1, # there is 1 M_3 in Eq(3)
                       2:7} # there are 7 M_3's in Eq(4)
        
        if len(test_cases)>0:
            print "\n===== Testing", fun_name, "====="

            for case in range(len(test_cases)):
                print "\n--- Test finding M_"+str(test_cases[case][1]), "in Eq("+str(test_cases[case][0])+") ---"
                Mns = Closure.findMn(test_cases[case][0],test_cases[case][1])
                Mns = Test.unique_items(Mns)
                print "\nThere are", len(Mns), "M_"+str(test_cases[case][1]), "in Eq("+str(test_cases[case][0])+")   Expected:", correct_ans[case]
                print "\nThey are:", Mns
                self.assertEquals(len(Mns), correct_ans[case], fun_name+" seems broken")
#                 
#         
        for N in range(3,11):
            for n in range(2,N+2):
                Mns = Closure.findMn(N,n)
                if len(Mns)==0:
                    break
                print N, n,
                #print 'nonunique:', len(Mns), 'unique:',
                UniqueMns = Test.unique_items(Mns)
                print len(UniqueMns)

#         for p in UniqueMns:
#             print p

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()