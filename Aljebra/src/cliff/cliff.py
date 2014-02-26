'''
Created on Feb 25, 2014
@author: williamdemeo@gmail.com
'''
from org.uacalc.alg.Malcev import isCongruenceDistIdempotent
from org.uacalc.alg.Malcev import cpIdempotent
from org.uacalc.io import AlgebraIO

homedir = "/home/williamdemeo/git/UACalc/"

outfile1 = open(homedir+"Algebras/CIB5/isCD.txt", 'w')
outfile2 = open(homedir+"Algebras/CIB5/isCP.txt", 'w')

for k in range(1,220):
    algname = "CIB5-"+str(k)
    algfile = homedir+"Algebras/CIB5/"+algname+".ua"
    A = AlgebraIO.readAlgebraFile(algfile)
    outfile1.write(algname+"   "+str(isCongruenceDistIdempotent(A, None))+"\n")
    outfile2.write(algname+"   "+str((cpIdempotent(A, None)==None))+"\n")
        
outfile1.close()
outfile2.close()


        