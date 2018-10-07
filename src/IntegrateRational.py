'''
Created on 07.10.2018

@author: Leonard
'''
from __future__ import division
import Polynomial as Pol
import RationalFunction as Rat
import FieldExtension as FE
import Integral as Int
import RootSum as RS
from Utils import *

from Matrix import Resultant
from Number import sqrt


def IntegratePolynomialPart(poly):
    """
    integrate poly el C[x]
    """
    if poly==0:
        return Int.Integral([],[])
    if poly.fieldTower.towerHeight!=0:
        raise Exception()
    
    Log("Integrate polynomial part: {}.".format(poly.printFull()))
    
    intPoly = Pol.Polynomial()
    for i in range(poly.degree+1):
        coeff = poly.getCoefficient(i)*(1/(i+1))
        intPoly.setCoefficient(coeff, i+1,callUpdates=False)
        
    intPoly.updateCoefficientsAll()
    
    Log("Integral of polynomial part: {}.".format(intPoly.printFull()))
    
    return Int.Integral(poly_rationals=[intPoly])
    #return Int.Integral(str(intPoly))

def IntegrateRationalFunction(func): # field=0, func element C(x)
    """
    integrate func el C(x)
    """
    if func==0:
        return Int.Integral()
    if isPoly(func):
        return IntegratePolynomialPart(func)
    if func.numerator.degree>=func.denominator.degree:
        (poly,rem) = Pol.PolyDiv(func.numerator,func.denominator)
        if rem==0:
            return IntegratePolynomialPart(poly)+IntegrateRationalFunction(0)
        else:
            return IntegratePolynomialPart(poly)+IntegrateRationalFunction(rem/func.denominator)
        
    Log("Integrate rational part: {}.".format(func.printFull()))
    sqrFree = func.denominator.isSquareFree()
    Log("Denominator is squarefree: {}".format(sqrFree))
    if sqrFree:
        Log("Denominator is squarefree -> use Rothstein/Trager method.")
        a = func.numerator
        b = func.denominator
        bp = b.differentiate()
        bpm = bp*(-1)
        
        zExtension = FE.FieldExtension(FE.TRANSCENDENTAL_SYMBOL,1,"z")
        newFieldTower = FE.FieldTower(zExtension)
        
        poly = Pol.Polynomial([a,bpm],fieldTower=newFieldTower)
        Log("Calculate resultant: res_x[{},{}]".format(poly,b)) 
        
        coeffsA = []
        Adeg = max(a.degree,bp.degree)   
        for i in range(Adeg+1):
            polZ = Pol.Polynomial([a.getCoefficient(i),bpm.getCoefficient(i)],fieldTower=newFieldTower)
            coeffsA.append(polZ)
            
        b_coeffs = b.getCoefficients()
        coeffsB = []
        for bc in b_coeffs:
            coeffsB.append(Pol.Polynomial([bc],fieldTower=newFieldTower))
        res = Resultant(coeffsA, coeffsB) # res_T (a-z*b',b)
        
        Log("Resultant: {}".format(res))
        
        if res!=0:
            primitivePart = res.makeMonic()
        if res==0 or res.isZero():
            primitivePart = 0
                
        sqrfreeFact = primitivePart.factorSquareFree()
        Log("Squarefree factorization of resultant: {}".format(Pol.printFactorization(sqrfreeFact)))
        integral = Int.Integral()
        zExtension = FE.FieldExtension(FE.TRANSCENDENTAL_SYMBOL,1,"z")
        newFieldTower = FE.FieldTower(zExtension)
        
        
        # gcd_x(a-z*b',b)
        coeffsA = []
        Adeg = max(a.degree,bp.degree)   
        for i in range(Adeg+1):
            polZ = Pol.Polynomial([a.getCoefficient(i),bp.getCoefficient(i)*(-1)],fieldTower=newFieldTower)
            coeffsA.append(polZ)
            
        b_coeffs = b.getCoefficients()
        coeffsB = []
        for bc in b_coeffs:
            coeffsB.append(Pol.Polynomial([bc],fieldTower=newFieldTower))
        #v = Pol._PolyGCD(coeffsA, coeffsB)
        
        for fac in sqrfreeFact:
            if fac[0]==1:
                continue
            d = fac[0].degree
            if d<=2:
                roots = fac[0].getRoots()
                for c in roots:
                    v = Pol.PolyGCD(a+(-c)*bp, b)
                    log = Int.LogFunction(v,c)
                    integral += Int.Integral(logs=[log])
            else:
                ratRoots = fac[0].getRationalRoots()
                if len(ratRoots)>0:
                    div = 1
                    for root in ratRoots:
                        p = Pol.Polynomial([-root,1],fieldTower=newFieldTower)
                        div *= p
                        sqrfreeFact.append((p,1))
                        
                    (s,r) = Pol.PolyDiv(fac[0],div)
                    if r!=0:
                        raise Exception(str(r))
                    sqrfreeFact.append((s,1)) # s has no rational roots, but it will be tested again -> TODO
                else:
                    v = Pol.Polynomial( _PolyGCDWithAlgebraicParameter(coeffsA, coeffsB, fac[0]),callUpdateCoeffs =False)
                    rootSum = RS.RootSum(fac[0],"z*{}".format(logExpression(v)), exprVar="z")
                    integral += Int.Integral(rootSums=[rootSum])
                    
        return integral
                
        #roots = primitivePart.getRoots()
        """if func.denominator.degree>2:
            poly = func.denominator
            coeff_rat = func.numerator/func.denominator.differentiate()
            coeff_str = coeff_rat.strCustomVar("y")
            rootSum = RS.RootSum(poly,"({}){}".format(coeff_str, logExpression("x-y")),exprVar="y")
            return Int.Integral([],[],[rootSum])
        """
    else:
        Log("Use Hermite-Reduction to make the denominator squarefree")
        return HermiteReduction(func)
    
def HermiteReduction(rational):
    """
    hermite reduction for rational el C(x), with rational=p/q, dep(p)<deg(q)
    """
    sqrFreeFactorization = rational.denominator.factorSquareFree()
    partialFractions = rational.PartialFraction(sqrFreeFactorization)
    integral = Int.Integral()
    for frac in partialFractions:
        j = frac[2]
        q_i = frac[1]
        r_ij = frac[0] # rational = sum(i=1...n, j=1...i : r_ij/q_i^j
        if j>1: # denominator = q_i^j not square free -> reduce power j -> make it squarefree
            (s,t) = Pol.extendedEuclidGenF(q_i, q_i.differentiate(), r_ij) # finds s,t with s*q_i+t*q_i'=r_ij
            tPrime = 0 if isNumber(t) else t.differentiate() # t'
            p1 = Int.Integral(poly_rationals=[(-1)*t/(j-1)/(q_i**(j-1))]) # integral(r_ij/q_i^j) = p1+p2 = -t/(j-1)/q_i^(j-1) + integral([s+t'/(j-1)]/[q_i^(j-1)])
            #print(t,q_i)
            num = s+tPrime/(j-1)
            if not num.isZero():
                r = num/(q_i**(j-1))
                p2 = IntegrateRationalFunction(r)
                integral += p1+p2
            else:
                integral+= p1
        else: # denominator = q_i is squarefree
            integral += IntegrateRationalFunction(r_ij/q_i)
    return integral


def _PolyDivWithAlgebraicParameter(coeffsA, coeffsB,poly):
    degA = len(coeffsA)-1
    degB = len(coeffsB)-1
    if degA<degB:
        return ([0],coeffsA)
    power = degA-degB
    coeff = coeffsA[-1]/coeffsB[-1]
    if isNumber(coeffsB[-1]):
        coeff = coeffsA[-1]/coeffsB[-1]
    else:
        if isPoly(coeffsB[-1]):
            coeff = coeffsA[-1]*coeffsB[-1].asRational().Inverse()
        else:
            coeff = coeffsA[-1]*coeffsB[-1].Inverse()
    monomial = [0]*power+[coeff] 
    sub = mulListWithList(monomial,mulObjectToListElements(-1,coeffsB))
    newA = addListToList(sub, coeffsA)
    newA = newA[0:len(newA)-1]
    
    zPolys = []
    for coeff in newA:
        if not coeff.isConstant() or not coeff.isZero():
            zPolys.append(coeff)
    print(zPolys)
    isZero = False
    if len(zPolys)==1:
        zPoly = zPolys[0]
        if type(zPoly)==Rat.RationalFunction:
            zPoly = zPoly.numerator
        (s,r) = Pol.PolyDiv(poly, zPoly)
        if r==0 or r.isZero():
            isZero = True
            
    if listIsZero(newA) or isZero:
        return (monomial,0)
    else:
        (quot, remainder)  = _PolyDivWithAlgebraicParameter(newA, coeffsB, poly)
        return (addListToList(quot,monomial),remainder)
def _PolyGCDWithAlgebraicParameter(coeffsA, coeffsB,poly):
    """
    calcs gcd of A,B that depend on an algebraic parameter z that is a zero of poly.
    """
    if listIsZero(coeffsA):
        return coeffsB
    if listIsZero(coeffsB):
        return coeffsA
    degA = len(coeffsA)-1
    degB = len(coeffsB)-1
    (A,B) = (coeffsA,coeffsB) if degA>=degB else (coeffsB,coeffsA)
    (s,r) = _PolyDivWithAlgebraicParameter(A, B, poly)
    r = 0 if r==0 else listStripZeros(r)
    if r==0 or listIsZero(r):
        return B
    gcd = _PolyGCDWithAlgebraicParameter(B,r,poly)
    return gcd

def logExpression(arg):
    return "log({})".format(str(arg))
