'''Closure class for computing the closure of a set of partitions.
@author: williamdemeo at gmail'''

from org.uacalc.alg.conlat import BasicPartition
from itertools import combinations
from org.uacalc.io import AlgebraIO
from org.uacalc.alg.op import Operation
from org.uacalc.alg import BasicAlgebra

class UnaryPolynomial:
    '''A simple class for representing a unary function'''
    def __init__(self, n=0):
        self.domain = n
        self.table = []
    
    # Return True if this operation's table is completely filled in.    
    def is_full(self):
        if len(self.table)==self.domain:
            return True
        return False
    
    # What fraction of this operation's table has been filled in?
    def progress(self):
        return len(self.table)/self.domain
    

class Closure(object):

    def __init__(self, parts):
        self.partitions = parts

        # sd_embedding will store the map from X to X/part_1 x X/part_2 x ... x X/part_r
        # where part_1, ..., part_r are the partitions in optimal_sdf_subset 
        # (not all the partitions in self.partitions) 
        self.sd_embedding = None  # will be an n x r array where X = {0, 1, ..., n-1}.
        self.has_sd_embedding = False

        # optimal_sdf_subset will store the partitions that should be used in the decomposition
        self.optimal_sdf_subset = None
        self.has_optimal_sdf_subset = False

        # Inv will store the partitions in the closure
        self.Inv = None
        self.has_Inv = False

        # Fix will store the unary functions that respect all partitions in self.partitions 
        self.Fix = None
        self.has_Fix = False
        

    # check if it's okay to have added k to position j of F[i] 
    # (F[i][j] = k means the ith function will map j to k.
    def in_range(self,F,i,j,k):
        if self.has_sd_embedding==False:
            self.compute_sd_embedding()
            
        n = len(self.sd_embedding)
        r = len(F)
        
        for x in range(n):
            y=[-1]*r
            # find each row of embedding where ith column is j 
            if (self.sd_embedding[x][i]==j):
            # (these must be mapped to a row of sd_embedding with k in the ith column)
                for p in range(r):
                    if len(F[p].table) < self.sd_embedding[x][p]:
                        # then we can't determine if it's okay to put k in pos j of F[i]
                        return -1
                    y[p] = F[p].table[self.sd_embedding[x][p]]
            # now check that the y vector is in the sd_embedding range of values
            if not (y in self.sd_embedding):
                return False
        
        return True


    def synthesize_function(self, F):
        '''Given F, the decompositions of a function, return that function'''
        if self.has_sd_embedding==False:
            self.compute_sd_embedding()
        n = len(self.sd_embedding)
        r = len(F)
        f = [-1]*n
        y = [-1]*r
        for k in range(n):
            sde = self.sd_embedding[k]
            for m in range(r):
                y[m] = F[k][sde[m]]
            for i in range(n):
                if (self.sd_embedding[i]==y):
                    f[k] = i
        if -1 in f:
            return -1
        return f
                    
    
    def Fix(self,F,FF):
        '''Recursively compute the set FF of all unary functions that respect the partitions in optimal_sdf_subset
        @param F: the decomposed version of the unary function we are currently building
        @param FF: the collection of unary functions already built that respect all partitions in optimal_sdf_embedding.'''
        i = Closure.get_short_index(F)
        
        # if the tables of all functions in F are completely filled in, add this F (or rather, the function f 
        # with decomposition F), to FF.
        if (i==-1):
            f = self.synthesize_function(F)
            if f==-1:
                print "Error: the unary function f we expected could not be constructed"
                return -1
            FF.append(f)
            return FF
        
        # otherwise, continue building F
        for k in range(F[i].domain):
            j = len(F[i].table)
            F[i].table.append(k)  # add k to position j of F[i]
            # check if this is okay (i.e. whether F[i].table[j] = k allows preserving partitions) 
            if self.in_range(F,i,j,k):
                # if so, leave it there and continue by recursion
                return Closure.Fix(F,FF)
            # otherwise, drop this k and continue (try k+1 next)
                F[i].table.pop()
    

    @staticmethod
    def get_short_index(F):
        '''Given F, a list of UnaryPolynomials, get the index of the element of F that is farthest from complete.'''
        answer = 0
        min_prog = 1
        for k in range(len(F)):
            if not F[k].is_full():
                if(F[k].progress()<min_prog):
                    min_prog = F[k].progress()
                    answer = k
        if min_prog==1:
            return -1  # this means all functions are full
        return answer
            
            
    @staticmethod
    def slow_closure(p):
        '''For a set of partitions p, compute the closure, that is, the set of 
        all partitions respected by all unary ops that respect all partitions in p.  
        This is the old, slower way of computing the closure.'''
        A = BasicPartition.unaryPolymorphismsAlgebra(p)
        print "|ConA| = ", len(A.con().universe())
        AlgebraIO.writeAlgebraFile(A, "/tmp/A.ua")
        return A.con()


    @staticmethod
    def basic_algebra_from_unary_polymorphisms(fns):
        '''Given a list fns, of unary functions on a set X = {0, 1,..., n-1}, 
        return the BasicAlgebra object representing <X, fns>'''
        f0 = fns[0]
        n = len(f0)
        ops = []
        # use the given functions in fns to construct UACalc operations
        for i in range(len(fns)):
            ops.append(Operation(fns[i], "f"+str(i), 1, n))

        #===========================================================================
        # # DEBUGGING: check that the operations give what we expect
        # print "Operations:"
        # for i in range(len(ops)):
        #     print "   " + ops[i].symbol().name() + ":", 
        #     for j in range(n):
        #         print ops[i].intValueAt([j]),
        #     print " "
        #===========================================================================

        # construct and return the algebra
        return BasicAlgebra("", n, ops)


    @staticmethod
    def maximal_element(ind, BL):
        '''Return a maximal element of the given basic lattice among those with indices in ind.'''
        pars = BL.getUniverseList()
        max_par = pars[ind[0]]
        for i in ind:
            if BL.leq(max_par, pars[i]):
                max_par = pars[i]
        return max_par

    @staticmethod
    def nary_meet(pars):
        '''Return the meet of all elements in the list pars'''
        ans = pars[0]
        for p in pars:
            ans = ans.meet(p)
        return ans

    @staticmethod
    def decomp_size(pars):
        '''Return a^a * b^b * c^c *... where a, b, c,... are the numbers of blocks in
        pars[0], pars[1], pars[2], ... respectively.'''
        answer = 1
        for p in pars:
            n = p.numberOfBlocks()
            answer = answer * (n**n)
        return answer

    def compute_sd_embedding(self):
        '''Return an array with (i,j) entry equal to the index of the block of pars[j] containing i.
            So, for example, if
                pars = |0,1|2,3|4,5|, |0,2|1,4|3,5|
            then return
                [[0,0], [0,1], [1,0], [1,2], [2,1], [2,2]]
            See Also
        '''
        answer = []
        n = self.partitions[0].universeSize()
        for i in range(n):
            temp = []
            for j in range(len(self.partitions)):
                temp.append(self.partitions[j].blockIndex(i))
            answer.append(temp)
        self.sd_embedding = answer
        self.has_sd_embedding = True
        return answer



    def compute_optimal_sdf_subset(self):
        '''Return a subset S of partitions pars that is optimal for computing
        the closure of pars using a subdirect factorization.
        For S to result in a subdirect decomposition, the partitions in S must meet to 0.
        For the subdirect decomposition to be optimal for computing the closure,
        decomp_size = a^a * b^b * c^c * ... should be as small as possible,
        where a, b, c,... denote the numbers of blocks in partitions S[0], S[1], S[2], ...
        :param pars:
        '''
        n = self.partitions[0].universeSize()
        pars = BasicPartition.joinClosure(self.partitions)
        N = len(pars)
        d_size = n**n  # the number to minimize
        #ZeRo = BasicPartition.zero(n)
        for k in range(2,N+1):
            parsubs = list(combinations(pars, k))
            for psub in parsubs:
                m = Closure.nary_meet(psub)
                if (m.isZero()):
                    temp = Closure.decomp_size(psub)
                    if (temp < d_size):
                        answer = psub
                        d_size = temp
        self.optimal_sdf_subset=answer
        self.has_optimal_sdf_subset=True
        return answer
