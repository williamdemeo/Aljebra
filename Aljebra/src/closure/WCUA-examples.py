'''Examples demonstrating the construction of UACalc algebras in Jython
for the Workshop on Computational Universal Algebra

Created on Oct 3, 2013

@see: OperationFactory.py
@author: williamdemeo@gmail.com and ralph@math.hawaii.edu
'''


# Example 1: Importing UACalc Classes

from org.uacalc import alg

ad = conlat.BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,18,19|8,9,20,21|10,11,22,23|30,31|32,33|34,35|")

from org.uacalc.alg import conlat


from org.uacalc.alg.conlat import BasicPartition

bd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7,30,31|8,9,32,33|10,11,34,35|18,19|20,21|22,23|")
cd = BasicPartition("|0,1,12,13,24,25|2,3,14,15,26,27|4,5,16,17,28,29|6,7|8,9|10,11|18,19,30,31|20,21,32,33|22,23,34,35|")
da = BasicPartition("|0,2,4,6,8,10|1,3,7,9|5,11|12,14,16,18,20,22|13,15,19,21|17,23|24,26,28,30,32,34|25,27,31,33|29,35|")
db = BasicPartition("|0,2,4,6,8,10|1,5,7,11|3,9|12,14,16,18,20,22|13,17,19,23|15,21|24,26,28,30,32,34|25,29,31,35|27,33|")
dc = BasicPartition("|0,2,4,6,8,10|1,7|3,5,9,11|12,14,16,18,20,22|13,19|15,17,21,23|24,26,28,30,32,34|25,31|27,29,33,35|")


print "\n---- Example 1 ----"
print "Constructing an algebra with operations defined using a Python function."

# define a function
def plus_mod5(args):
    result = 0
    for x in args:
        result = result + x
    return result % 5

     
# use the function above to construct UACalc operations
op0 = Operation(plus_mod5, "binaryPlusMod5", 2, 5)

from OperationFactory import Operation

op1 = Operation(plus_mod5, "ternaryPlusMod5", 3, 5)

print "Quick sanity check of the operations:"
print "   binaryPlusMod5:  4 + 10 mod 5 = ", op0.intValueAt([4,10])
print "   ternaryPlusMod5:  4 + 10 + 1 mod 5 = ", op1.intValueAt([4,10,1])

# make a list of the operations we want in the algebra
ops = op0, op1

# construct the algebra
alg = BasicAlgebra("MyAlgebra", 5, ops)

print "Quick check that we constructed an algebra:"
print "   alg.getName() = ", alg.getName()        
print "   alg.universe() = ", alg.universe()

import os.path
if os.path.exists("../Algebras"):

    # write the algebra to a UACalc file
    AlgebraIO.writeAlgebraFile(alg, "../Algebras/Example1_ConstructAlgebra.ua")
    print "UACalc algebra file created: ../Algebras/Example1_ConstructAlgebra.ua"




print "\n\n---- Example 2 ----"
print "Constructing an algebra with unary operations defined 'by hand' as vectors."

# The algebra will have universe {0, 1, ..., 7}, and the following unary operations: 
f0 = 7,6,6,7,3,2,2,3
f1 = 0,1,1,0,4,5,5,4
f2 = 0,2,3,1,0,2,3,1

# Use the functions above to construct UACalc operations:
op0 = Operation(f0, "f0", 1, 8)
op1 = Operation(f1, "f1", 1, 8)
op2 = Operation(f2, "f2", 1, 8)

# Make a list of the operations we want in the algebra:
ops = op0, op1, op2

# The above is nice and easy, but the following is more general
# and useful if you want to define many operations:
# Suppose all the operations are stored as a list of lists, as in:
fns = f0, f1, f2
# Build the list of operations by looping through fns:
ops = []
for i in range(len(fns)):
    ops.append(Operation(fns[i], "f"+str(i), 1, 8))

    
# Next we check that the operations are as expected:
print "The unary operations are:"
for i in range(len(ops)):
    print "   " + ops[i].symbol().name() + ":", 
    for j in range(8):
        print ops[i].intValueAt([j]),
    print " "


# Finally, construct the algebra:
alg = BasicAlgebra("MyUnaryAlgebra", 8, ops)

print "\nQuick check that we constructed an algebra:"
print "   alg.getName() = ", alg.getName()        
print "   alg.universe() = ", alg.universe()

if os.path.exists("../../Algebras"):
    fqname2 = "../../Algebras/Example2_MutliunaryAlgebra.ua"
    # Optionally, write the algebra to a UACalc file that can be loaded into the gui.
    AlgebraIO.writeAlgebraFile(alg, fqname2)
    print "UACalc algebra file created:", fqname2
    


print "\n\n---- Example 3 ----"
print "Constructing an algebra with unary operations defined 'by hand' as vectors."
print "Quick and dirty method: write to a .ua file, then read with AlgebraIO"

def write_unaryalgebra(ops, fqname):
    '''Write a .ua file representing a multiunary algebra with operations given
    by the vectors in the argument ops.'''
    outfile = open(fqname, 'w')  # open file for appending (might change this to 'w')
    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<algebra>\n<basicAlgebra>\n<algName>"+name+"</algName>\n")
    outfile.write("<cardinality>"+ str(len(ops[0]))+"</cardinality>\n")
    outfile.write("<operations>\n")
    opcount=0
    
    for op in ops:
        outfile.write("<op>\n<opSymbol>\n<opName>op"+str(opcount)+"</opName>\n")
        outfile.write("<arity>1</arity>\n</opSymbol>\n<opTable>\n<intArray>\n<row>")
        for x in range(len(op)-1):
            outfile.write(str(op[x])+",")
        outfile.write(str(op[-1]))
        outfile.write("</row>\n</intArray>\n</opTable>\n</op>\n")
        opcount=opcount+1;
    outfile.write("</operations>\n</basicAlgebra>\n</algebra>\n")
    outfile.close()

# The algebra will have universe {0, 1, ..., 7}, and the following unary operations: 
f0 = 7,6,6,7,3,2,2,3
f1 = 0,1,1,0,4,5,5,4
f2 = 0,2,3,1,0,2,3,1
ff = f0, f1, f2

name = "Example3"    
if os.path.exists("../../Algebras"):
    fqname3 = '../../Algebras/'+name+'.ua'
else:
    fqname3 = name+'.ua'

print "Writing algebra to file ", fqname3
write_unaryalgebra(ff, fqname3)


alg = AlgebraIO.readAlgebraFile(fqname3)

print "\nQuick check that we constructed an algebra:"
print "   alg.getName() = ", alg.getName()        
print "   alg.universe() = ", alg.universe()
    




print "\n\n---- Example 4 ----"
print "The congruence lattice of a congruence lattice."
print "\nThe congruence lattice Con(A) of the algebra A is itself an algebra (specifically, a lattice)."
print "We represent it in UACalc as an object of the class BasicLattice.  UACalc represents the" 
print "universe of an algebra of cardinality n with integers {0, 1, ..., n-1}, and in some cases"
print "it is important to know to what elements these integers correspond."

print "\nFor example, suppose we read in Polin's algebra from the file polin.ua, and name this algebra P,"
P = AlgebraIO.readAlgebraFile("../../Algebras/polin.ua")

print "and suppose we then construct an algebra that is the congruence lattice of P (using the convenient"
print "and fast UACalc con() method), and call this algebra conP."
conP = BasicLattice("conP", P.con(), 0)

print "\nWe can print the universe of conP with"
print "\n    conP.universe(): ", conP.universe()

print "\nNow suppose we want the congruence lattice of the algebra conP.  "
print "\n    conP.con().universe() gives ", conP.con().universe()

print "\nTo make use of this result, we must know to which elements of conP.universe() the"
print "integers appearing in the blocks of conP.con().universe() correspond."

n = conP.cardinality()
print "\nThe elements in the universe of conP happen to be labeled in UACalc as follows:"
print "\n    ( k, conP.getElement(k) ):\n   ",
for k in range(n):
    print "( "+str(k)+",", conP.getElement(k), "), ",

print "\n\nPrinting the universe of conP..."
print "   ...with conP.universe():", conP.universe(), "   (lists elements in arbitrary order)"
print "   ...with conP.getUniverseList():", conP.getUniverseList(), "   (lists elements in correct order)"

print "\nSo, when you want to know how UACalc is labelling the elements,"
print "use either getUniverseList() or getElement(k), instead of universe()."

    
