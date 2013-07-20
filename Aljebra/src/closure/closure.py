from org.uacalc.alg.conlat import BasicPartition
from itertools import chain, combinations



class Closure(object):


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

    def set_partitions(self, partitions):
        self.partitions = partitions


    def sd_embedding(self):
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
        return answer



    def optimal_sdf_subset(self,pars):
        '''Return a subset S of partitions pars that is optimal for computing
        the closure of pars using a subdirect factorization.
        For S to result in a subdirect decomposition, the partitions in S must meet to 0.
        For the subdirect decomposition to be optimal for computing the closure,
        decomp_size = a^a * b^b * c^c * ... should be as small as possible,
        where a, b, c,... denote the numbers of blocks in partitions S[0], S[1], S[2], ...
        :param pars:
        '''
        pars = BasicPartition.joinClosure(pars)
        N = len(pars)
        n = pars[0].universeSize()
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
        return answer
