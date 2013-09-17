'''Utilities for constructing UACalc algebra objects in Jython.

Created on Jun 18, 2013
@see: AlgebraConstructionExample.py
@author: williamdemeo at gmail
'''
import sys
sys.path.append("Jars/uacalc.jar")
sys.path.append("Jars/LatDraw.jar")


from org.uacalc.alg.op import AbstractOperation

def operation_class_factory(fn, op_class_name, arity):
    '''Return a subclass of AbstractOperation with intValueAt() defined using the function fn.
    
    :param fn: the Python function used to define the intValueAt() method of the class.
    :param op_class_name: a string name for the class.
    :param arity: the arity of the operation (must be 1 if fn is a list or tuple)
    '''
    class OperationClass(AbstractOperation):
        # check if fn is a (callable) function
        if hasattr(fn, '__call__'):
            def intValueAt(self,args):
                return fn(args)
        # check if fn is a list or tuple, in which case the op is unary
        elif type(fn)==list or type(fn)==tuple:
            if arity!=1:
                raise Exception("If function is a vector (list or tuple), arity of operation must be 1.")
            def intValueAt(self,args):
                return fn[args[0]]
        else:
            raise Exception("Unsupported function type.")
            
    # Give the class a unique name:
    OperationClass.__name__= op_class_name
    return OperationClass


def Operation(fn, op_name, arity, univ_size):
    '''Return a UACalc Operation object defined using fn.
    
    :param fn: the Python function that will be used to define the UACalc operation.
    :param op_name: a string name of the resulting UACalc operation.
    :param arity: an integer arity of the operation.
    :param univ_size: an integer cardinality of the universe of the algebra.
    '''
    opClass = operation_class_factory(fn, op_name, arity)
    op = opClass(op_name, arity, univ_size)
    return op

'''
Example Usage: Constructing an algebra
'''
from org.uacalc.alg import BasicAlgebra
from org.uacalc.io import AlgebraIO

if __name__ == '__main__':

    # define a function
    def plus_mod5(args):
        result = 0
        for x in args:
            result = result + x
        return result % 5
     
    # use it to construct some UACalc operations
    op0 = Operation(plus_mod5, "binaryPlusMod5", 2, 5)
    op1 = Operation(plus_mod5, "ternaryPlusMod5", 3, 5)

    # quick check that the operations give what we expect
    print "4 + 10 mod 5 = ", op0.intValueAt([4,10])
    print "4 + 10 + 1 mod 5 = ", op1.intValueAt([4,10,1])

    # make a list of the operations we want in the algebra
    ops = op0, op1

    # construct the algebra
    alg = BasicAlgebra("MyAlgebra", 5, ops)

    # quick check that we actually constructed an algebra
    print "alg.getName() = ", alg.getName()        
    print "alg.universe() = ", alg.universe()

    import os.path
    if os.path.exists("../Algebras/"):
    
        # write the algebra to a UACalc file:
        AlgebraIO.writeAlgebraFile(alg, "../Algebras/Example_ConstructAlgebra.ua")

