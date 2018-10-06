'''
Created on 06.10.2018

@author: Leonard
'''
from __future__ import division
from Matrix import QuadMatrix
import Polynomial as Pol
import RationalFunction as Rat

def solveRischEquation(p,q, fieldTower):
    if fieldTower.towerHeight>0:
        raise NotImplementedError()
    if type(p) is Rat.RationalFunction:
        pn = p.asPolynomial()
        if pn!=None:
            p = pn
    if type(q) is Rat.RationalFunction:
        qn = q.asPolynomial()
        if qn!=None:
            q = qn
    if not type(p) is Pol.Polynomial or not type(q) is Pol.Polynomial:
        raise NotImplementedError()
    return solveComplexPolynomials(p, q)
    
def solveComplexPolynomials(p,q):
    """
    Solves y'+py=q with p,q,y el C[x]
    """
    deg = q.degree-p.degree
    if deg<0:
        return None
    #Ansatz: y=a_l *x^l+a_(l-1)*x^(l-1)+...+a_0, l=deg
    return __helperSolve(p, q, deg)
    
def __helperSolve(p,q,degree):
    if degree==0:
        c = q.getLeadingCoefficient()/p.getLeadingCoefficient()
        pn = c*p
        if pn==q:
            return Pol.Polynomial([c])
        else:
            return None
        
    a = p.getLeadingCoefficient()
    b = q.getLeadingCoefficient()
    y_lc = b/a
    mon = Pol.Monomial(degree, y_lc)
    qn = q + (-mon.differentiate()) + (-p * mon)
    
    deg = qn.degree-p.degree
    if deg<0:
        return None
    #Ansatz: y=a_l *x^l+a_(l-1)*x^(l-1)+...+a_0, l=deg
    s = __helperSolve(p, qn, deg)

    #s = __helperSolve(p, qn, degree-1)
    if s==None:
        return None
    else:
        return s + mon
    
    

if __name__ == '__main__':
    x = Pol.Polynomial([0,1])
    x.replaceNumbersWithRationals()
    one = Pol.Polynomial([1])
    one.replaceNumbersWithRationals()
    print(solveComplexPolynomials(2*x, x)) # 0.5
    print(solveComplexPolynomials(one, x)) # x-1
    print(solveComplexPolynomials(x, x**3)) # x^2-2
    print(solveComplexPolynomials(x**2,x**5)) # x^3-3