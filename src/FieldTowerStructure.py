'''
Created on 08.10.2018

@author: Leonard
'''
from __future__ import division
"""
Helper methods to test whether a given logarithm/exponential is transcendendal or algebraic over a given field tower.
"""

import Matrix
import FieldExtension as FE
import Polynomial as Pol
import RationalFunction as Rat
import random

from Number import Rational


def logIsTranscendental(u, fieldTower):
    """
    returns: log(u) is transcendental over fieldTower (True/False)
    let log(u_j) the logarithmic extensions of fieldTower.
    Then: log(u) is algebraic over fieldTower iff there are k,k_j el Z, k!=0 where the product u^k * product(u_jk_j) is constant
    """
    # 0 = k*(u'/u) + sum(k_j*(u_j'/u))
    
    logs = [x.argFunction.reduceToLowestPossibleFieldTower() for x in fieldTower.getLogarithms()]
    #if sum(x.fieldTower.towerHeight for x in logs)>0:
     #   raise NotImplementedError()
    
    N = len(logs)
    if N==0:
        return not u.isConstant()
    
    """log_diffs = [u.logDifferentiate()]+[log.logDifferentiate() for log in logs]
    
    e = True
    while e:
        testPoints = [Rational(random.randint(0,1000),1000) for _ in range(N+1)]
        coefficients = Matrix.QuadMatrix(N+1)
        try:
            for i in range(N+1):
                for j in range(N+1):
                    coefficients.setElement(i,j,log_diffs[j].completeEvaluate(testPoints[i]))
        except ZeroDivisionError:
            continue
        e = False
        print(coefficients)
        (solution, unique) = Matrix.solveLinearSystem(coefficients, [0]*(N+1))
        print(solution,unique)"""
            
    
    vs = [u.differentiate()*reduce(lambda x,y:x*y, [logs[j] for j in range(N)])]+  [logs[i].differentiate()*u*reduce(lambda x,y:x*y, [1]+[logs[j] for j in range(N) if j!=i]) for i in range(N)]
    #v = u.differentiate()*reduce(lambda x,y:x*y, [logs[j] for j in range(N)])
    terms = []
    one = Rational(1,1)
    #terms.append([v]+[one]+[0]*N)
    i = 0
    todoTerms = zip(vs,range(N+1))
    while i<len(todoTerms):
        (term,ind) = todoTerms[i]
        if type(term)==Rat.RationalFunction:
            t = term.asPolynomial()
            if id(t)!=id(None):
                term=t
        # do it with multivariate polynomials -> TODO
        
    """for term in vs:
        x  =logHelperMultipleInList(term, terms)
        if x==False:
            terms.append([term] + [0]*(i+1) + [one]+ [0]*(N-i-1))
        else:
            (t,c) = x # c= term/t[0]
            t[2+i] += c
            
        i += 1"""
        
    """print(terms)
    M = Matrix.QuadMatrix(N+1)
    for i in range(N+1):
         for j in range(N+1):
            if i>=len(terms):
                M.setElement(i,j,0)
            else:
                M.setElement(i,j,terms[i][j+1])
                
    print(M)
    s = Matrix.solveLinearSystem(M, [0]*(N+1))
    print(s)"""
    """if type(v)==Rat.RationalFunction:
        terms.append
    else:
        for t in v._coefficients:
            lc = """
    
    
def logHelperMultipleInList(term, terms):
    for t in terms:
        c = term.isConstantMultipleOf(t[0])
        if id(c)!=id(False):
            return (t,c) # c = term/t[0]
    return False
if __name__ == '__main__':
    #tests
    one = Rational(1,1)
    #fe1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([one,one])/Pol.Polynomial([0,one]),"T_1")
    fe1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,one]),"T_1") # log(x)
    ft = FE.FieldTower(fe1)
    fe2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,one],fieldTower=ft),"T_2") # log(log(x))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fe1,fe2])
    u = Pol.Polynomial([0,0, one])*Pol.Polynomial([0,one],fieldTower=ft)
    logIsTranscendental(u, FE.fieldTower)