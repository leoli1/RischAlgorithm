'''
Created on 06.10.2018

@author: Leonard
'''
from __future__ import division
import Polynomial as Pol
import RationalFunction as Rat
import FieldExtension as FE
from Number import Rational

def solveRischEquation(p,q, fieldTower):
    if fieldTower.towerHeight>1:
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
    
  
# general solution, not finished (yet)
def getADBEG(p,q):
    A = p.numerator
    D = p.denominator
    B = q.numerator
    E = q.denominator
    
    G = Pol.PolyGCD(D,E)
    return (A,D,B,E,G)

def getDenominator(f,g):
    """
    Finds the denominator of the solution y of y'+py=q
    f=A/D, g=B/E, y=Q/T
    G = gcd(D,E)
    T = gcd(E,dE/DT)/gcd(G,dG/dT)
    """
    (A,D,B,E,G) = getADBEG(f, g)
    
    #dE_dT = E.differentiateWRTtoPolyVar()
    #dG_dT = G.differentiateWRTtoPolyVar()
    
    T = E.getLogGCD()/G.getLogGCD()#Pol.PolyGCD(E, dE_dT) / Pol.PolyGCD(G,dG_dT)
    
    return T

def getNumeratorEquation(f,g):
    """
    returns the coefficients of the equation for the numerator of y
    DTQ' + (AT-DT')Q = BD(T^2)/E
    returns (a,b,c) with aQ' + bQ = c
    """
    (A,D,B,E,G) = getADBEG(f, g)
    T = getDenominator(f, g)
    
    (s,r) = Pol.PolyDiv(D*T*T, E)
    if r==0 or r.isZero():
        raise NoRischDESolutionError("E doesn't divide D *T^2")
    return (D*T, A*T-D*T.differentiate(), B*s)

def getFinalNumeratorEquation(f,g,fieldTower):
    (a,b,c) = getNumeratorEquation(f, g)
    k = a.lowestDegree
    T_pk = Pol.Monomial(k, Rational(1,1), fieldTower=fieldTower)
    a = a/T_pk
    b = b/T_pk
    c = c/T_pk
    return (a,b,c)
    #fieldExtension = fieldTower.getLastExtension()
    #if fieldExtension==None or fieldExtension.extensionType!=FE.TRANS_EXP:
    #    return (a,b,c)
    
def find_m(a,b,c):
    """
    finds the smallest number m such that q=T^m * Q el F[Q] (i.e. so that q doesn't contain negative powers, only for exponential extensions)
    """
    pass


class NoRischDESolutionError(Exception):
    pass
if __name__ == '__main__':
    x = Pol.Polynomial([0,1])
    x.replaceNumbersWithRationals()
    one = Pol.Polynomial([1])
    one.replaceNumbersWithRationals()
    print(solveComplexPolynomials(2*x, x)) # 0.5
    print(solveComplexPolynomials(one, x)) # x-1
    print(solveComplexPolynomials(x, x**3)) # x^2-2
    print(solveComplexPolynomials(x**2,x**5)) # x^3-3