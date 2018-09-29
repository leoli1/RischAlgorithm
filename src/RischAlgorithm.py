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

def Integrate(func, fieldTower):
    (quot, rem) = Pol.PolyDiv(func.numerator, func.denominator)
    rat = Rat.RationalFunction(rem,func.denominator)
    print("{} = [{}]+[{}]".format(str(func),str(quot),str(rat)))
    lastExtension = fieldTower.getLastExtension()
    if lastExtension == None:
        return IntegrateRationalFunction(func)
    elif lastExtension.extensionType==FE.TRANS_LOG:
        out1 = IntegratePolynomialPartLogExt(quot, fieldTower)
        out2 = IntegrateRationalPartLogExt(rat, fieldTower)
    elif lastExtension.extensionType==FE.TRANS_EXP:
        out1 = IntegratePolynomialPartExpExt(quot, fieldTower)
        out2 = IntegrateRationalPartExpExt(rat, fieldTower)
    
    if out1==None or out2==None: return None
    return out1+out2

def IntegrateRationalFunction(func): # field=0, func element C(x)
    if func.numerator.degree>=func.denominator.degree:
        raise NotImplementedError()
    if func.denominator.isSquareFree():
        #partial fraction
        if func.denominator.degree>2:
            
            poly = func.denominator
            coeff_rat = func.numerator/func.denominator.differentiate()
            coeff_str = coeff_rat.strCustomVar("y")
            rootSum = RS.RootSum(poly,"{}{}".format(coeff_str, logExpression("x-y")),exprVar="y")
            return Int.Integral(str(rootSum))
        elif func.denominator.degree==1:
            arg = func.denominator
            coeff = func.numerator/func.denominator.getLeadingCoefficient()
            return Int.Integral("{}*{}".format(str(coeff), logExpression(arg)))
            #return LF.LogFunction(func.denominator)*(func.numerator/func.denominator.getLeadingCoefficient())
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
            #return (LF.LogFunction(denomA)*A)+(LF.LogFunction(denomB)*B)
    else:
        raise NotImplementedError
def IntegratePolynomialPartLogExt(func, fieldTower):
    raise NotImplementedError
def IntegrateRationalPartLogExt(func, fieldTower):
    if func.numerator==0:
        return 0
    raise NotImplementedError
def IntegratePolynomialPartExpExt(func, fieldTower):
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
    from Parse import parseField0PolyFromStr
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