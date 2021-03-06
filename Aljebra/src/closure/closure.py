'''Closure class for computing the closure of a set of partitions.
@author: williamdemeo at gmail'''

from org.uacalc.alg.conlat import BasicPartition
from itertools import combinations
from org.uacalc.io import AlgebraIO
from OperationFactory import Operation
#from org.uacalc.alg.op import Operation
from org.uacalc.alg import BasicAlgebra


class ClosureIndexError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
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
        return float(len(self.table))/float(self.domain)

    # Allow easy printing a list L of UnaryPolynomial objects (just print the table)
    def __str__(self):
        return str(self.table)
    def __repr__(self):
        return str(self)
    def __getitem__(self,i):
        return self.table[i]
        

def debug_print(args):
    if Closure.GLOBAL_DEBUG:
        for m in args:
            print m,
        print "\n"



class Closure(object):

    GLOBAL_DEBUG = False

    def __init__(self, parts):

        self.partitions = parts

        self.universe_size = parts[0].universeSize()

        # sd_embedding: the map from X to X/part_1 x X/part_2 x ... x X/part_r where part_1,...
        # ..., part_r are the partitions in optimal_sdf_subset (not all partitions in self.partitions)
        # Example: for partitions |0,1|2,3|4,5| and |0,2|1,4|3,5|, 
        # sd_embedding is [[0,0], [0,1], [1,0], [1,2], [2,1], [2,2]] 
        self.sd_embedding = None  # will be an n x r array where X = {0, 1, ..., n-1}.
        self.has_sd_embedding = False

        # optimal_sdf_subset: the partitions that should be used in the decomposition
        self.optimal_sdf_subset = None
        self.has_optimal_sdf_subset = False

        # nonoptimal_idemdecs: partitions in self.partitions that are not in optimal_sdf_subset.
        # N.B. this field stores paritions in their idemdec representation.
        self.nonoptimal_idemdecs = None
        self.has_nonoptimal_idemdecs = False
        
        # Inv: the partitions in the closure
        self.Inv = None
        self.has_Inv = False

        # Fix will store the unary functions that respect all partitions in self.partitions 
        self.Fix = None
        self.has_Fix = False



    def subarray(self, i,j):
        ''' return the rows of the sd_embedding that have a j in the ith column'''
        if self.has_sd_embedding==False:
            self.compute_sd_embedding()
        indx = [p for p, row in enumerate(self.sd_embedding) if row[i]==j]
        return [self.sd_embedding[i] for i in indx]

    # check if it's okay to have added k to position j of F[i] 
    # (F[i][j] = k means the ith function will map j to k.
    def in_range(self,F,i,j,k):
        if self.has_sd_embedding==False:
            self.compute_sd_embedding()
            
        r = len(F)

        Dom = self.subarray(i,j)  # get those rows with j in column i
        Ran = self.subarray(i,k)  # get those rows with k in column i
        
        # Putting the value k at position F[i][j] means the ith function in F maps j to k.
        # The set of rows of sd_embedding with j in column i must get mapped to 
        # the set of rows of sd_embedding with k in column i.
        ColsD = zip(*Dom)
        ColsR = zip(*Ran)
        
        fun_indx = range(r)
        del fun_indx[i]     # Don't want to check the ith function!
        for p in fun_indx:
            for x in ColsD[p]:
                if x < len(F[p].table) and not (F[p].table[x] in ColsR[p]):
                    return 0
        debug_print([">>>>> Returning 1 <<<<<"])
        return 1


    def synthesize_function(self, F):
        '''Return the function f that has subdirectly decomposed representation F'''
        # Don't need the next two lines.  (Couldn't get here without having an sd_embedding, I think.)
        #if self.has_sd_embedding==False:
        #   self.compute_sd_embedding()
        # n = len(self.sd_embedding)
        n = self.universe_size
        r = len(F)
        f = [-1]*n
        y = [-1]*r
        for x in range(n):
            sde = self.sd_embedding[x]
            for m in range(r):
                y[m] = F[m].table[sde[m]]
            for i in range(n):
                if (self.sd_embedding[i]==y):
                    f[x] = i
        if -1 in f:
            return -1
        return f
                    
    
    def isRespector(self, f):
        '''Check whether the map F respects the remaining partitions in self.partitions.'''
        # n = len(f)
        n = self.universe_size
        for p in self.partitions:
            if p not in self.optimal_sdf_subset:
                for x in range(n-1):
                    u = f.table[x]
                    rx = p.representative(x);
                    for y in range(x+1,n):
                        if (rx == p.representative(y)) and not p.isRelated(u, f.table[y]):
                            return False
        return True
    

    def get_vector_value(self, F, s):
        '''Return [ F[0][s[0]], F[1][s[1]], ..., F[r][s[r]] ]'''
        ans = []
        assert len(F)==len(s)
        for i in range(len(F)):
            if len(F[i].table) <= s[i]:
                raise ClosureIndexError("F["+str(i)+"][x] not defined for x="+str(s[i]))
            else:
                ans.append(F[i].table[s[i]])
        return ans
    
    
    def respects_nonoptimal_idemdecs(self, F):
        
        if self.has_nonoptimal_idemdecs:
            for t in self.nonoptimal_idemdecs:
                for x in range(self.universe_size):
                    s = self.sd_embedding[x]
                    y = self.sd_embedding[t[x]]
                    try:
                        # F is [F[0], F[1], ..., F[r-1]], a decomposed version of a function f on X.
                        u = self.get_vector_value(F,s)  # the sde of f(x)
                        v = self.get_vector_value(F,y)  # the sde of f(t(x)) 
                        inv_u = self.sd_embedding.index(u)  # the inverse of the sde of f(x)
                        inv_v = self.sd_embedding.index(v)  # the inverse of the sde of f(t(x))

                        if not t[inv_u]==t[inv_v]:
                            # t(f(t(x))) is not equal to t(f(x)) then f doesn't respect t:
                            return False
                    
                    except ClosureIndexError as e:
                        debug_print(e)

        return True
                    
    
    
    def compute_sd_Fix(self,F,FF):
        '''Recursively compute the set FF of all unary functions that respect the partitions in optimal_sdf_subset
        @param F: the decomposed version of the unary function we are currently building
        @param FF: the collection of unary functions already built that respect all partitions in optimal_sdf_embedding.'''
        if self.has_sd_embedding==False:
            self.compute_sd_embedding()
        if len(F)==0:
            # insert r new empty UnaryPolynomials into F
            for k in range(len(self.optimal_sdf_subset)):
                # Fhe k-th unary polynomial is a map from the set {0, 1, ..., b-1} to itself,
                # where b is the number of blocks of the k-th partition in optimal_sdf_subset:
                F.append(UnaryPolynomial(self.optimal_sdf_subset[k].numberOfBlocks()))
            i=0  # we will start by building up the 0-th member of F
            debug_print(["F: ",F])
        else: 
            # Otherwise, we build up the shortest member of F
            # Find out which one that is: 
            i = Closure.get_short_index(F)
        
        if (i==-1): # this means the tables of all functions in F are completely filled in.
            # In this case, we add this F (or rather, the function f with decomposition F), to FF, and return.
            f = self.synthesize_function(F)
            if f==-1:
                print "ERROR: the unary function f we expected could not be constructed"
                return -1

            # no longer need this, since our in_range function now does all the work (tests all partitions respected)
            # if self.isRespector(f):
                # FF.append(f)

            FF.append(f)
            return FF
        
        # otherwise, work on the shortest function in F (the one at index i)
        for k in range(F[i].domain):
            j = len(F[i].table)
            F[i].table.append(k)  # add k to position j of F[i]
            debug_print(["F: ",F])

            # check if this is okay (i.e. whether F[i].table[j] = k allows preserving partitions
            # in optimal_sdf_subset and respects other paritions as well)
            if self.in_range(F,i,j,k) and self.respects_nonoptimal_idemdecs(F):
                # if so, leave it there and continue by recursion
                FF = self.compute_sd_Fix(F,FF)
            # otherwise, drop this k and continue (try k+1 next)
            F[i].table.pop()

        return FF

    @staticmethod
    def get_short_index(F):
        '''Given F, a list of UnaryPolynomials, get the index of the element of F that is farthest from complete.'''
        answer = 0
        min_prog = 1
        for k in range(len(F)):
            debug_print(["F["+str(k)+"].progress() = "+str(F[k].progress())])
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
    def partition2idemdec(p):
        n = p.universeSize()
        ans = [-1]*n
        for x in range(n):
            ans[x] = p.representative(x)
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
        '''Return an array with (i,j) entry equal to the index of the block of pars[j] containing i,
        where pars is the set of partitions stored in optimal_sdf_subset.
        So, for example, if pars = |0,1|2,3|4,5|, |0,2|1,4|3,5|
        then this function should return [[0,0], [0,1], [1,0], [1,2], [2,1], [2,2]].
        '''
        if self.has_optimal_sdf_subset == False:
            self.compute_optimal_sdf_subset()
        
        answer = []
        n = self.optimal_sdf_subset[0].universeSize()
        for i in range(n):
            temp = []
            for j in range(len(self.optimal_sdf_subset)):
                blkindx = self.optimal_sdf_subset[j].blockIndex(i)
                temp.append(blkindx)
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

        # The self.partitions that are not in optimal_sdf_subset:
        nonopts = [p for p in self.partitions if p not in answer]
        if len(nonopts) > 0:
            self.has_nonoptimal_idemdecs=True
            # store the idemdec representations:
            self.nonoptimal_idemdecs = []
            for p in nonopts:
                self.nonoptimal_idemdecs.append(Closure.partition2idemdec(p))

        return answer
            
    
    
