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
import MultivariatePolynomial
import Integral as Int
from Number import Rational,ONE
from itertools import product

def getLinearLogCombination(logs, log):
    """
    IMPORTANT: this function ignores the (possible) coefficients of logs if they are of type Integral.LogFunction
    logs/log can be of type Integral.LogFunction or FieldExtension.FieldExtension
    tests whether log is a (rational) linear combination of the logs
    return False if there is no such combination and (True, [...]) if there is one and [...] are the coefficients
    e.g. log(x/(x+1)) = log(x) + (-1)*log(x+1) -> returns (True, [1,-1])
         log(x^2) = 2*log(x)                   -> returns (True, [2])
         log(x+1) != k*log(x) for all k el Z (or Q, R)
    """
    # step 1: construct a field tower with the given logs -> note: the field tower only contains the logs and not other functions that might be necessary to define these logs
    fieldTower = FE.FieldTower()
    for i in range(len(logs)):
        l = logs[i]
        fieldTower = fieldTower.addFieldExtension(FE.FieldExtension(FE.TRANS_LOG, l.argFunction, FE.Variable("L_{}".format(i))))
    
    # step 2: log is rational linear combination of the logs iff log is algebraic over this field tower
    trans = logIsTranscendental(log.argFunction, fieldTower)
    if trans==True:
        return False
    else: # step 3: 
        rationalCombo = [-trans[1][i]/trans[1][0] for i in range(1,len(logs)+1)]
        return (True,rationalCombo)

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
            return (False,[1])
        else:
            return True
            
    
    vs = [u.differentiate().reduceToLowestPossibleFieldTower()*reduce(lambda x,y:x*y, [logs[j] for j in range(N)])]+  [logs[i].differentiate().reduceToLowestPossibleFieldTower()*u*reduce(lambda x,y:x*y, [1]+[logs[j] for j in range(N) if j!=i]) for i in range(N)]
    
    p = ONE
    for v in vs:
        if type(v)==Rat.RationalFunction:
            p *= v.denominator
            for c in v.numerator._coefficients: #### TODO!
                if type(c)==Rat.RationalFunction:
                    p *= c.denominator
                    
    p = p.reduceToLowestPossibleFieldTower()
    vs = [v.reduceToLowestPossibleFieldTower()*p for v in vs]
    #vs = 

    vs_emp = [MultivariatePolynomial.ExtendedMultivariatePolynomial.fromNormalPolynomial(v) for v in vs] # convert the vs to ExtendedMultivariatePolynomials which will make combining alike terms easier
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
    
    #M = Matrix.QuadMatrix(N+1)
    M = Matrix.Matrix(data=k_coeffs)
                
    s = Matrix.solveLinearSystem(M, [0]*(len(k_coeffs)))
    if s==None:
        raise Exception()
    else:
        unique = s[1]
        if s[0][0]==0 and unique:
            return True
        elif s[0][0]!=0:
            return (False, s[0])
        else:
            f = s[2]
            for i in product(range(1,4),repeat=s[3]):
                if f(i)[0]!=0:
                    return (False,f(i))
            return True
            #raise NotImplementedError("no unique solution of linear system, particular solution has k=0, TODO: is there any other solution with k!=0?")
        
def logarithmIsInFieldTower(log, fieldTower):
    return logIsTranscendental(log.argFunction, fieldTower)!=True
    
    
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
    print(getLinearLogCombination([fe1,fe2], FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,0,ONE]), FE.Variable('A'))))