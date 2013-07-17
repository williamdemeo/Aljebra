''' Utilities for finding congruence lattice representations of small lattices.
Created on Jun 12, 2013
@author: williamdemeo at gmail
@see: Itertools documentation at http://docs.python.org/2/library/itertools.html
'''
from org.uacalc.alg.conlat import BasicPartition
from org.uacalc.lat import BasicLattice
from org.uacalc.alg.conlat import CongruenceLattice
from org.uacalc.io import AlgebraIO
from org.uacalc.alg import BasicAlgebra
from OperationFactory import Operation
import itertools
from itertools import chain, combinations

# For a set of partitions p, compute the closure, that is,
# the set of all partitions respected by all unary ops that
# respect all partitions in p.
def SlowClosure(p):
    A = BasicPartition.unaryPolymorphismsAlgebra(p)
    print "|ConA| = ", len(A.con().universe())
    AlgebraIO.writeAlgebraFile(A, "/tmp/A.ua")
    return A.con()

def maximal_element(ind, BL):
    '''Return a maximal element of the given basic lattice among those with indices in ind.'''
    pars = BL.getUniverseList()
    max_par = pars[ind[0]]
    for i in ind:
        if BL.leq(max_par, pars[i]):
            max_par = pars[i]
    return max_par


def all_unary_maps(n):
    '''Return a list of all n^n unary maps from the set {0,1,...,n-1} to itself.
    Examples:
        all_unary_maps(2) --> [(0, 0), (0, 1), (1, 0), (1, 1)]
        all_unary_maps(3) --> [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 1), (0, 1, 2),...
    Recipes:
        All pairs of maps from {0,1,2} to {0,1,2}. 
            F = all_unary_maps(3)
            pairs = list(itertools.combinations(F,2))
    @see http://docs.python.org/2/library/itertools.html
    '''
    return list(itertools.product(range(n),repeat=n))

def powerset(iterable):
    '''Return a of all subsets of a set.
    Example: list(powerset([1,2,3])) --> [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    @see http://docs.python.org/2/library/itertools.html
    '''
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def subsets(s):
    return map(set, powerset(s))

def nary_meet(pars):
    '''Return the meet of all elements in the list pars'''
    ans = pars[0]
    for p in pars:
        ans = ans.meet(p)
    return ans

def decomp_size(pars):
    '''Return a^a * b^b * c^c *... where a, b, c,... are the numbers of blocks in 
    pars[0], pars[1], pars[2], ... respectively.'''
    answer = 1
    for p in pars:
        n = p.numberOfBlocks()
        answer = answer * (n**n)
    return answer

def optimal_sdf_subset(pars):
    '''Return a subset S of partitions pars that is optimal for computing 
    the closure of pars using a subdirect factorization.
    For S to result in a subdirect decomposition, the partitions in S must meet to 0.
    For the subdirect decomposition to be optimal for computing the closure,
    decomp_size = a^a * b^b * c^c * ... should be as small as possible, 
    where a, b, c,... denote the numbers of blocks in partitions S[0], S[1], S[2], ...
    '''
    pars = BasicPartition.joinClosure(pars)
    N = len(pars)
    n = pars[0].universeSize()
    d_size = n**n  # the number to minimize
    #ZeRo = BasicPartition.zero(n)
    for k in range(2,N+1):
        parsubs = list(combinations(pars, k))
        for psub in parsubs:
            m = nary_meet(psub)
            if (m.isZero()):
                temp = decomp_size(psub)
                if (temp < d_size):
                    answer = psub
                    d_size = temp
    return answer

def sd_embedding(pars):
    '''Return an array with (i,j) entry equal to the index of the block of pars[j] containing i.'''
    answer = []
    n = pars[0].universeSize()
    for i in range(n):
        temp = []
        for j in range(len(pars)):
            temp.append(pars[j].blockIndex(i))
        answer.append(temp)
    return answer


class PartialFunction:
    def __init__(self, dom, ran):
        self.dom = dom
        self.ran = ran
        self.map = []

    def isFull(self):
        if (self.deficiency()==0):
            return 0
        return 1

    def deficiency(self):
        return self.dom - len(self.map)
    
    def getMap(self):
        return self.map
        
        
class DecomposedUnaryMap:
    
    def __init__(self, partitions):
        self.partitions = partitions
        self.maps = []
        for p in partitions:
            self.maps.append(PartialFunction(p.numberOfBlocks(),p.numberOfBlocks()))
        
    def numberOfParts(self):
        return len(self.partitions)
        
    def isFull(self):
        for f in self.maps:
            if not f.isFull():
                return 0
        return 1

    def indexOfMaxDeficiency(self):
        maxdef = 0
        answer = 0
        for i in range(self.numberOfParts()):
            d =self.maps[i].deficiency() 
            if d>maxdef:
                maxdef=d
                answer = i
        return answer
    
    def getMap(self,i):
        return self.maps[i].getMap()
        

def build_unary_polys(H, F, pars):
    '''Example:
    oppars = 
    H = DecomposedUnaryMap(pars)
    F = build_unary_polys(H,[],pars)
    '''
    if H.isFull():
        # the partial functions in H are all actual functions
        F.append(H)
        return F

    # Get index of partial function in H that is farthest from being a full function
    i = H.indexOfMaxDeficiency()
    for k in range(H.maps[i].range):
        # make a temp copy of map i ---make sure we are really getting a copy!!!
        trial_map = H.getMap(i)
        trial_map.append(k)
        if respects(trial_map, F, pars):
            H.maps[i].map.append(k)
            return build_unary_polys(H,F,pars)
        

def respects(trial_f, F, pars):
    '''check if the maps in trial_f and F respect all partitions in pars'''
        
        
'''UNIT TESTS'''
def setup_partitions(example):
    if (example==100):
        # A first, very basic test:
        a0 = BasicPartition("|0|1,2,3,4|5,6|7,8|")
        a1 = BasicPartition("|0,1,2|3,4,5|6,7,8|")
        a2 = BasicPartition("|0,1,8|2,3|4,5|6,7|")
        a3 = BasicPartition("|0,1,2,3,4|5,6,7,8|")
        pars = a0, a1, a2, a3
        return pars
    # Example 1 (PJ11):
    if (example=='optimal_sdf_subset_test0'):
        H = BasicPartition("|0,1,6|2,7,35|3,22,36|4,8,17|5,9,18|19,26,56|10,20,57|11,21,58|12,45,59|13,46,60|16,25,39|51,66,96|32,70,84|14,23,37|15,24,38|34,48,62|28,43,80|40,49,77|41,50,78|33,47,61|65,74,95|27,42,79|54,90,99|29,44,81|86,91,103|30,68,82|31,69,83|87,92,104|64,73,94|75,88,105|52,67,97|100,102,107|55,71,85|76,101,106|53,89,98|63,72,93|")
        A = BasicPartition("|0,1,6,15,24,38,55,71,85|2,7,11,21,29,35,44,58,81|3,22,36,53,54,89,90,98,99|4,8,16,17,25,33,39,47,61|5,9,14,18,23,34,37,48,62|19,26,40,49,56,63,72,77,93|10,20,28,43,52,57,67,80,97|12,13,45,46,59,60,76,101,106|27,42,51,66,75,79,88,96,105|30,31,32,68,69,70,82,83,84|41,50,64,73,78,86,91,94,103|65,74,87,92,95,100,102,104,107|")
        B = BasicPartition("|0,1,4,5,6,8,9,14,15,16,17,18,23,24,25,33,34,37,38,39,47,48,55,61,62,71,85|2,7,10,11,20,21,27,28,29,35,42,43,44,51,52,57,58,66,67,75,79,80,81,88,96,97,105|3,12,13,22,30,31,32,36,45,46,53,54,59,60,68,69,70,76,82,83,84,89,90,98,99,101,106|19,26,40,41,49,50,56,63,64,65,72,73,74,77,78,86,87,91,92,93,94,95,100,102,103,104,107|")
        C = BasicPartition("|0,1,2,3,6,7,19,22,26,35,36,56|4,8,10,12,17,20,40,45,49,57,59,77|5,9,11,13,18,21,41,46,50,58,60,78|16,25,29,32,39,44,65,70,74,81,84,95|33,47,51,53,61,66,86,89,91,96,98,103|14,23,27,30,37,42,63,68,72,79,82,93|15,24,28,31,38,43,64,69,73,80,83,94|34,48,52,54,62,67,87,90,92,97,99,104|55,71,75,76,85,88,100,101,102,105,106,107|")
        K = BasicPartition("|0,4,14|1,9,25|2,51,52|3,13,32|5,15,33|6,61,62|7,43,88|22,45,68|8,24,48|26,73,102|10,11,75|19,86,87|12,31,54|27,28,29|16,34,55|35,57,79|36,83,106|17,18,85|20,44,66|21,42,67|56,78,95|46,69,89|37,38,39|70,90,101|23,47,71|49,74,91|50,72,92|40,41,100|63,64,65|30,53,76|58,80,96|59,84,98|60,82,99|77,94,104|81,97,105|93,103,107|")
        pars = H,A,B,C,K
        return pars
    # Example 2 (parallel sum of M3's):
    if (example=='optimal_sdf_subset_test1' or example=='sd_embedding_test1'):
        ad = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,18,19|8,9,20,21|10,11,22,23|30,31|32,33|34,35|")
        bd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,30,31|8,9,32,33|10,11,34,35|18,19|20,21|22,23|")
        cd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7|8,9|10,11|18,19,30,31|20,21,32,33|22,23,34,35|")
        da = BasicPartition("|0,2,4,6,8,10|1,3,7,9|5,11|12,14,16,18,20,22|13,15,19,21|17,23|24,26,28,30,32,34|25,27,31,33|29,35|")
        db = BasicPartition("|0,2,4,6,8,10|1,5,7,11|3,9|12,14,16,18,20,22|13,17,19,23|15,21|24,26,28,30,32,34|25,29,31,35|27,33|")
        dc = BasicPartition("|0,2,4,6,8,10|1,7|3,5,9,11|12,14,16,18,20,22|13,19|15,17,21,23|24,26,28,30,32,34|25,31|27,29,33,35|")
        pars = ad,bd,cd,da,db,dc
        return pars
    if (example=='sd_embedding_test0'):
        # A first, very basic test:
        pi0 = BasicPartition("|0,1|2,3|4,5|")
        pi1 = BasicPartition("|0,2|1,4|3,5|")
        pars = pi0, pi1
        return pars
    print "Error: unrecognized example identifier."
    return None


def test_sd_embedding(pars):
    answer = sd_embedding(pars)
    '''For first test, should be: [[0,0], [0,1], [1,0], [1,2], [2,1], [2,2]]'''
    print "\n    The subdirect embedding is:", answer
    
def test_optimal_sdf_subset(pars):
    N = pars[0].universeSize()
    print "pars: "
    for p in pars:
        n = p.numberOfBlocks()
        print "    ", p
        print "     num blocks:", n, "   nb^nb:", n**n
    oppars = optimal_sdf_subset(pars)
    red_num = 1
    print "\noptimal subset:"
    for p in oppars:
        n = p.numberOfBlocks()
        red_num = red_num * (n**n)
        print "    ", p
        print "     num blocks:", n, "   nb^nb:", n**n
        #print "    ", p, "  num blocks:", n, "   nb^nb:", n**n
    print "\nTotal number of unary maps:  |X|^|X| = ", str(N) + "^" + str(N), " = %.3g" % N**N  
    print "Number of decomposed unary maps: a^a * b^b * c^c *... = %.3g" % red_num


'''Here is where we actually run the tests.'''
if 1:
    count=0

    test_name='optimal_sdf_subset_test0'
    print "\n\n--------------- TEST", count, test_name
    test_optimal_sdf_subset(setup_partitions(test_name))
    count=count+1
    
    test_name = 'optimal_sdf_subset_test1'
    print "\n\n--------------- TEST", count, test_name
    test_optimal_sdf_subset(setup_partitions(test_name))
    count=count+1
    
    test_name = 'sd_embedding_test0'
    print "\n\n--------------- TEST", count, test_name,
    test_sd_embedding(setup_partitions(test_name))

    test_name = 'sd_embedding_test1'
    print "\n\n--------------- TEST", count, test_name
    pars = optimal_sdf_subset(setup_partitions(test_name))
    print "     Using optimal partition set with numbers of blocks:", 
    for p in pars:
        print p.numberOfBlocks(),
    test_sd_embedding(pars)

    
if 0:
    a0 = BasicPartition("|0,1|2|3|4|")
    a1 = BasicPartition("|0|1|2|3,4|")
    a2 = BasicPartition("|0,1|2|3,4|")
    a3 = BasicPartition("|0,1,2|3,4|")
    a4 = BasicPartition("|0,1|2,3,4|")
    pars = a1, a3, a4
    print "nary_meet(pars)= |0|1|2|3,4| =?= ", nary_meet(pars)
    pars = a0, a1, a2, a3, a4
    
    parpower = powerset(pars)
    parpowerlist = list(parpower)
    print parpowerlist
    print "|P| = ", len(parpowerlist)
    print "2^5 = ", 2**5
    for s in parpower:
        print "s = ", s
        print "len(s) = ", len(s)
        
    print list(powerset([1,2,3]))
    #L = Closure(pars)
    #A = BasicPartition.unaryPolymorphismsAlgebra(pars)
    #print "Con(A) (", A.con().cardinality(), "):", A.con().universe()
    #L = BasicLattice("ConA", A.con(), 0)
    # AlgebraIO.writeAlgebraFile(A, "/tmp/A.ua")
    # AlgebraIO.writeAlgebraFile(L, "/tmp/ConA.ua")



'''Test conjecture about homomorphic images of congruence lattices'''
def test_conjecture(L,debug):
    Luni = L.getUniverseList()
    # IMPORTANT: you must use BL.getUniverseList() here,
    # instead of BL.universe(), if you want the elements 
    # to be in the order corresponding to the numbers in 
    # the congruence lattice.

    ConL = CongruenceLattice(L)  # Con(Con(A))

    if debug:
        print "\nL: "
        print "    ", L.getUniverseList()
        print "Con(L) = "
        print "    ", ConL.universe()
        print "Con(L).getAlgebra() = "
        print "    ", ConL.getAlgebra().getUniverseList()

    for theta in ConL.universe():
        if debug:
            print "\n\n------ theta: ", theta, " -----------"
        blks= theta.getBlocks()
        kept_pars = []
        #print "theta.getBlocks():", blks
        for b in blks:
            if len(b)>1:
                m = maximal_element(b, L)
                part = BasicPartition(m.toString())
                kept_pars.append(part)
            else:
                part = BasicPartition(Luni[b[0]].toString())
                kept_pars.append(part)

        if debug:
            print "kept_pars (", len(kept_pars), "): ", kept_pars
        if len(kept_pars) > 2:
            A = BasicPartition.unaryPolymorphismsAlgebra(kept_pars)
            CL = BasicLattice("ConA", A.con(), 0)
            if debug:
                print "closure   (", CL.cardinality(), "): ", CL.universe()
            if CL.cardinality()!=len(kept_pars):
                print "\nCOUNTEREXAMPLE FOUND!!!\n\n "
    

'''Example 1'''
if 0:
    alphat = BasicPartition("|0,1,2,6|3,5|4|7,8|")
    alpstr = BasicPartition("|0,1,2,6|3|5|4|7,8|")
    bethat = BasicPartition("|0,3,4,8|1|2,7|5,6|")
    betstr = BasicPartition("|0,3,4,8|1|2|7|5,6|")
    gamhat = BasicPartition("|1,4,5,7|2,3|0|6,8|")
    gamstr = BasicPartition("|1,4,5,7|2,3|0|6|8|")

    pars = alphat, alpstr, bethat, betstr, gamhat, gamstr
    A = BasicPartition.unaryPolymorphismsAlgebra(pars)
    #print "Con(A): ", A.con().universe()
    L = BasicLattice("ConA", A.con(), 0)
    # AlgebraIO.writeAlgebraFile(A, "/tmp/A.ua")
    # AlgebraIO.writeAlgebraFile(L, "/tmp/ConA.ua")

    L = SlowClosure(pars)  # Con(A)
    print "L = ", L.universe()
    BL = BasicLattice("ConA",L,0)

    test_conjecture(BL,1)
            
            
            
        
'''Example 2: PJ19'''
if 0:
    #7,6,6,7,3,2,2,3
    #0,1,1,0,4,5,5,4
    #0,2,3,1,0,2,3,1

    print "\n---- Example 2 ----"
    # the algebra will have universe {0, 1, ..., 7}, and the following unary operations: 
    f0 = 7,6,6,7,3,2,2,3
    f1 = 0,1,1,0,4,5,5,4
    f2 = 0,2,3,1,0,2,3,1

    # turn these into functions (not very elegant, but it works)
    def fun0(args):
        return f0[args[0]]
    def fun1(args):
        return f1[args[0]]
    def fun2(args):
        return f2[args[0]]

    # use the functions above to construct UACalc operations
    op0 = Operation(fun0, "f0", 1, 8)
    op1 = Operation(fun1, "f1", 1, 8)
    op2 = Operation(fun2, "f2", 1, 8)

    # check that the operations give what we expect
    print "Operations:"
    print "  f0:", 
    for i in range(8):
        print op0.intValueAt([i]),
    print "\n  f1:", 
    for i in range(8):
        print op1.intValueAt([i]),
    print "\n  f2:", 
    for i in range(8):
        print op2.intValueAt([i]),

    # make a list of the operations we want in the algebra
    ops = op0, op1, op2

    # construct the algebra
    A = BasicAlgebra("MyUnaryAlgebra", 8, ops)

    # quick check that we constructed an algebra
    print "\nA.getName() = ", A.getName()        
    print "A.universe() = ", A.universe()


    BL = BasicLattice("ConA", A.con(), 0)
    # AlgebraIO.writeAlgebraFile(A, "/tmp/A.ua")
    # AlgebraIO.writeAlgebraFile(L, "/tmp/ConA.ua")

    test_conjecture(BL, 1)
    
    
    
    
'''Example 3: PJ31'''
if 0:
    print "\n---- Example 3 ----"

    # the algebra will have universe {0, 1, 2, 3, 4}, and the following unary operations: 
    f0 = (0,1,1,0,0) 
    f1 = (1,1,2,2,2) 
    f2 = (3,2,2,4,4)
    fns = f0, f1, f2
    
    ops = []
    # use the functions above to construct UACalc operations
    for i in range(len(fns)):
        ops.append(Operation(fns[i], "f"+str(i), 1, 5))

    # check that the operations give what we expect
    print "Operations:"
    for i in range(len(ops)):
        print "   " + ops[i].symbol().name() + ":", 
        for j in range(5):
            print ops[i].intValueAt([j]),
        print " "

    # construct the algebra
    A = BasicAlgebra("MyUnaryAlgebra", 5, ops)

    # quick check that we constructed an algebra
    print "\nA.getName() = ", A.getName()        
    print "A.universe() = ", A.universe()


    BL = BasicLattice("ConA", A.con(), 0)
    AlgebraIO.writeAlgebraFile(A, "/tmp/A.ua")
    # AlgebraIO.writeAlgebraFile(L, "/tmp/ConA.ua")
    
    test_conjecture(BL, 1)