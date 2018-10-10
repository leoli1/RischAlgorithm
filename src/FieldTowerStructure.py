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
import MultivariatePolynomial

from Number import Rational,ONE


def logIsTranscendental(u, fieldTower):
    """
    returns: log(u) is transcendental over fieldTower (True/False), if False -> returns list [k,k1,k2....] such that u^k * product(u_j^k_j) is constant
    let log(u_j) the logarithmic extensions of fieldTower.
    Then: log(u) is algebraic over fieldTower iff there are k,k_j el Z, k!=0 where the product u^k * product(u_j^k_j) is constant
    """
    # 0 = k*(u'/u) + sum(k_j*(u_j'/u))
    
    logs = [x.argFunction.reduceToLowestPossibleFieldTower() for x in fieldTower.getLogarithms()]
    
    N = len(logs)
    if N==0:
        if u.isConstant():
            return False
        else:
            return (True, [1])
            
    
    vs = [u.differentiate()*reduce(lambda x,y:x*y, [logs[j] for j in range(N)])]+  [logs[i].differentiate()*u*reduce(lambda x,y:x*y, [1]+[logs[j] for j in range(N) if j!=i]) for i in range(N)]


    vs_emp = [MultivariatePolynomial.ExtendedMultivariatePolynomial.fromNormalPolynomial(v) for v in vs]
    multiIndices = []
    k_coeffs = []
    for i in range(N+1):
        v_emp = vs_emp[i]
        for term in v_emp.terms:
            try:
                index = multiIndices.index(term.multiIndex)
                k_coeffs[index][i] += term.coefficient
            except ValueError:
                multiIndices.append(term.multiIndex)
                k_coeffs.append([0]*i + [term.coefficient]+[0]*(N-i))
    
    M = Matrix.QuadMatrix(N+1)
    for i in range(N+1):
         for j in range(N+1):
            if i>=len(k_coeffs):
                M.setElement(i,j,0)
            else:
                M.setElement(i,j,k_coeffs[i][j])
                
    s = Matrix.solveLinearSystem(M, [0]*(N+1))
    if s==None:
        raise Exception()
    else:
        unique = s[1]
        if s[0][0]==0 and unique:
            return True
        elif s[0][0]!=0:
            return (False, s[0])
        else:
            raise NotImplementedError("no unique solution of linear system, particular solution has k=0, TODO: is there any other solution with k!=0?")
        
def logarithmIsInFieldTower(log, fieldTower):
    return logIsTranscendental(log.argFunction, fieldTower)==False
    
    
if __name__ == '__main__':
    #tests
    fe1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,ONE]),FE.Variable('T_1')) # log(x)
    ft = FE.FieldTower(fe1)
    fe2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,ONE],fieldTower=ft),FE.Variable('T_2')) # log(log(x))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fe1,fe2])
    u = Pol.Polynomial([0,0, ONE])*Pol.Polynomial([0,ONE],fieldTower=ft) # x^2*log(x)
    print(logIsTranscendental(u, FE.fieldTower)) # log(x^2*log(x)) = 2*log(x)+log(log(x))
    
    u = Pol.Polynomial([0,0,ONE])+Pol.Polynomial([0,ONE],fieldTower=ft) # x^2+log(x)
    print(logIsTranscendental(u, FE.fieldTower))