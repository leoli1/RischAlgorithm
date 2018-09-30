'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import Polynomial as Pol
import RationalFunction as Rat
import FieldExtension as FE
#import LogFunction as LF
import Integral as Int
import RootSum as RS
from math import sqrt as msqrt
from Utils import isNumber,isPoly

from Matrix import Resultant


logExtensionsInIntegral = []

def Integrate(func, fieldTower):
    (quot, rem) = Pol.PolyDiv(func.numerator, func.denominator)
    rat = Rat.RationalFunction(rem,func.denominator,fieldTower=func.getFieldTower())
    #print("{} = [{}]+[{}]".format(str(func),str(quot),str(rat)))
    lastExtension = fieldTower.getLastExtension()
    try:
        if lastExtension == None:
            out1 = IntegratePolynomialPart(quot)
            out2 = IntegrateRationalFunction(rat)
        elif lastExtension.extensionType==FE.TRANS_LOG:
            out1 = IntegratePolynomialPartLogExt(quot, fieldTower)
            out2 = IntegrateRationalPartLogExt(rat, fieldTower)
        elif lastExtension.extensionType==FE.TRANS_EXP:
            out1 = IntegratePolynomialPartExpExt(quot, fieldTower)
            out2 = IntegrateRationalPartExpExt(rat, fieldTower)
    except Int.IntegralNotElementaryException as e:
        if type(e) == Int.IntegralNotElementaryException:
            return Int.Integral("Integral is not elementary")
        else:
            print(e)
            raise e
        out1 = None
        out2 = None
    
    if out1==None or out2==None: return None
    return out1+out2

def IntegratePolynomialPart(poly):
    if poly==0:
        return Int.Integral("")
    if poly.field!=0:
        raise Exception()
    intPoly = Pol.Polynomial()
    for i in range(poly.degree+1):
        coeff = poly.getCoefficient(i)/(i+1)
        intPoly.setCoefficient(coeff, i+1)
    return Int.Integral(str(intPoly))

def IntegrateRationalFunction(func): # field=0, func element C(x)
    if isPoly(func):
        return IntegratePolynomialPart(func)
    if func.numerator.degree>=func.denominator.degree:
        (poly,rat) = Pol.PolyDiv(func.numerator,func.denominator)
        return IntegratePolynomialPart(poly)+IntegrateRationalFunction(rat/func.denominator)
    if func.denominator.isSquareFree():
        #partial fraction
        if func.denominator.degree>2:
            poly = func.denominator
            coeff_rat = func.numerator/func.denominator.differentiate()
            coeff_str = coeff_rat.strCustomVar("y")
            #logExtensions = []
            #for i in range(func.denominator.degree):
            rootSum = RS.RootSum(poly,"{}{}".format(coeff_str, logExpression("x-y")),exprVar="y")
            return Int.Integral(str(rootSum))
        elif func.denominator.degree==1:
            arg = func.denominator
            coeff = func.numerator/func.denominator.getLeadingCoefficient()
            return Int.Integral("{}*{}".format(str(coeff), logExpression(arg)))
        elif func.denominator.degree==2:
            f = func.numerator
            g = func.denominator
            a = g.getLeadingCoefficient()
            p = g.getCoefficient(1)/a
            q = g.getCoefficient(0)/a
            disc = sqrt((p/2)**2-q)
            if disc==0:
                raise Exception
            zero0 = -p/2+disc
            zero1 = -p/2-disc
            A = f.evaluate(zero0)/g.differentiate().evaluate(zero0)
            B = f.evaluate(zero1)/g.differentiate().evaluate(zero1)
            denomA = Pol.Polynomial([-zero0,1])
            denomB = Pol.Polynomial([-zero1,1])
            one = Pol.Polynomial([1])
            a = Rat.RationalFunction(one,denomA)
            b = Rat.RationalFunction(one,denomB)
            return IntegrateRationalFunction(a*A)+IntegrateRationalFunction(b*B)
    else:
        return HermiteReduction(func)#raise NotImplementedError
    
def HermiteReduction(rational):
    sqrFreeFactorization = rational.denominator.factorSquareFree()
    partialFractions = rational.PartialFraction(sqrFreeFactorization)
    integral = Int.Integral("")
    for frac in partialFractions:
        j = frac[2]
        q_i = frac[1]
        r_ij = frac[0]
        if j>1:
            (s,t) = Pol.extendedEuclidGenF(q_i, q_i.differentiate(), r_ij) # s*q_i+t*q_i'=r_ij
            tPrime = 0 if isNumber(t) else t.differentiate()
            p1 = Int.Integral(str((-1)*t/(j-1)/(q_i**(j-1))))
            #print(t,q_i)
            num = s+tPrime/(j-1)
            if not num.isZero():
                r = num/(q_i**(j-1))
                p2 = IntegrateRationalFunction(r)
                integral += p1+p2
            else:
                integral+= p1
        else:
            integral += IntegrateRationalFunction(r_ij/q_i)
    return integral
    
def IntegratePolynomialPartLogExt(func, fieldTower):
    if func == 0 or func.isZero():
        return Int.Integral("")
    raise NotImplementedError
def IntegrateRationalPartLogExt(func, fieldTower): # Hermite Recuction
    if func.numerator==0:
        return Int.Integral("")
    func.makeDenominatorMonic()
    sqrFreeFactorization = func.denominator.factorSquareFree()
    partialFractions = func.PartialFraction(sqrFreeFactorization)
    integral = Int.Integral("")
    for frac in partialFractions:
        j = frac[2]
        q_i = frac[1]
        r_ij = frac[0]
        if j>1: 
            (s,t) = Pol.extendedEuclidGenF(q_i, q_i.differentiate(), r_ij)
            tPrime = 0 if isNumber(t) else t.differentiate()
            p1 = Int.Integral(str((-1)*t/(j-1)/(q_i**(j-1))))
            num = s+tPrime/(j-1)
            if not num.isZero():
                r = num/(q_i**(j-1))
                p2 = IntegrateRationalPartLogExt(r, fieldTower)
                integral += p1+p2
            else:
                integral += p1
        else:
            a = r_ij
            b = q_i
            bp = q_i.differentiate()
            
            zExtension = FE.FieldExtension(FE.TRANSCENDENTAL_SYMBOL,1,"z")
            newFieldTower = fieldTower.getStrippedTower(fieldTower.towerHeight)
            newFieldTower.addFieldExtension(zExtension)
            
            coeffsA = []
            Adeg = max(a.degree,bp.degree)   
            for i in range(Adeg+1):
                polZ = Pol.Polynomial([a.getCoefficient(i),bp.getCoefficient(i)*(-1)],fieldTower=newFieldTower)
                coeffsA.append(polZ)
                
            b_coeffs = b.getCoefficients()
            coeffsB = []
            for bc in b_coeffs:
                coeffsB.append(Pol.Polynomial([bc],fieldTower=newFieldTower))
            res = Resultant(coeffsA, coeffsB) # res_T (a-z*b',b)
            if res!=0:
                print(res)
                primitivePart = res.makeMonic()
                constantRoots = primitivePart.hasOnlyConstantCoefficients()
            if res==0 or res.isZero():
                constantRoots = False
                primitivePart = 0
                    
            print("Only constant coefficients in primitive part {}: {}".format(primitivePart,constantRoots))
            
            if not constantRoots:
                raise Int.IntegralNotElementaryException("Integral is not elementary")
            
            roots = primitivePart.getRoots()
            vfunc = []
            for c in roots:
                vfunc.append(Pol.PolyGCD(a+bp*(-c), b))
            intPart = Int.Integral("")
            for (c,v) in zip(roots,vfunc):
                intPart += Int.Integral("{}{}".format(c,logExpression(v)))
            integral += intPart
    return integral
            
    #raise NotImplementedError
def IntegratePolynomialPartExpExt(func, fieldTower):
    if func == 0 or func.isZero():
        return Int.Integral("")
    raise NotImplementedError
def IntegrateRationalPartExpExt(func, fieldTower):
    if func.numerator ==0:
        return 0
    raise NotImplementedError

def sqrt(x):
    if x<0:
        return complex(0,1)*(msqrt(-x))
    return msqrt(x)

def logExpression(arg):
    return "log({})".format(str(arg))
if __name__ == '__main__':
    from Parse import parseField0PolyFromStr,parseField0RatFromStr
    """FE.fieldTower = FE.FieldTower()
    polA = Pol.Polynomial([1])
    polB = Pol.Polynomial([2,-3,1])
    polC = Pol.Polynomial([-1,2])
    print(Integrate(Rat.RationalFunction(polA,polC), FE.FieldTower()))
    ratA = Rat.RationalFunction(polA,polB)
    print(ratA)
    print(Integrate(ratA, FE.FieldTower()))
    print("----")
    polD = parseField0PolyFromStr("3*x**2+x+1")
    polE = parseField0PolyFromStr("x**4+x**2+1")
    ratB = Rat.RationalFunction(polD,polE)
    print(ratB)
    print(Integrate(ratB,FE.FieldTower()))
    print("----")
    polA = parseField0PolyFromStr("2*x**4+x**2+x+(-1)")
    polB = parseField0PolyFromStr("(-1)+x")*(parseField0PolyFromStr("x+1")**3)#*parseField0PolyFromStr("x**2+x-1")**3
    ratA = Rat.RationalFunction(polA,polB)
    print(ratA)
    print(Integrate(ratA, FE.FieldTower()))
    print("----")"""
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,1]),"T_1") # field extension with log(x)
    """FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    one = parseField0RatFromStr("1/1")
    ratX = parseField0RatFromStr("x/1")
    
    numerator = Pol.Polynomial([one],field=1)
    denom = Pol.Polynomial([0,ratX],field=1)
    rat = Rat.RationalFunction(numerator,denom,field=1) # 1/(x log(x))
    print(rat.printFull())
    print(Integrate(rat,FE.fieldTower).printFull())
    
    numerator = Pol.Polynomial([one],field=1)
    denom = Pol.Polynomial([ratX,one],field=1)
    rat = Rat.RationalFunction(numerator,denom,field=1) # 1/(x + log(x))
    print(rat.printFull())
    print(Integrate(rat,FE.fieldTower).printFull())
    
    numerator = Pol.Polynomial([one],field=1)
    denom = Pol.Polynomial([0,one,one],field=1)
    rat = Rat.RationalFunction(numerator,denom,field=1) # 1/(log(x)^2 + log(x))
    print(rat.printFull())
    print(Integrate(rat,FE.fieldTower).printFull())"""
    
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,1],fieldTower=FE.FieldTower(fieldExtension1)),"T_2") # field extension with log(log(x))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1,fieldExtension2])
    numerator = Pol.Polynomial([1],fieldTower=FE.fieldTower)
    #ratXf1 = Rat.RationalFunction(Pol.Polynomial([ratX]),field=1)
    polX = Pol.Polynomial([parseField0PolyFromStr("x")],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,polX*Pol.Polynomial([0,1],fieldTower=FE.fieldTower.getStrippedTower(1))],fieldTower=FE.fieldTower)
    rat = Rat.RationalFunction(numerator,denom,fieldTower=FE.fieldTower) # 1/(x log(x) log(log(x)))
    print(rat.printFull())
    print(Integrate(rat,FE.fieldTower).printFull())
    
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,0,1]),"T") # field extension with log(x^2)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    """numerator = Pol.Polynomial([one],field=1)
    denom = Pol.Polynomial([0,one],field=1)
    rat = Rat.RationalFunction(numerator,denom,field=1) # 1/(log(x^2))
    print(rat.printFull())
    print(Integrate(rat,FE.fieldTower).printFull())"""