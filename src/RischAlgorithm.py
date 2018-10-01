'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import Polynomial as Pol
import RationalFunction as Rat
import FieldExtension as FE
import Integral as Int
import RootSum as RS
from math import sqrt as msqrt
from Utils import isNumber,isPoly

from Matrix import Resultant


logExtensionsInIntegral = []

def Integrate(func, fieldTower):
    if isPoly(func):
        quot = func
        rat = 0
    else:
        (quot, rem) = Pol.PolyDiv(func.numerator, func.denominator)
        rat = Rat.RationalFunction(rem,func.denominator,fieldTower=func.getFieldTower())
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
            return None#Int.Integral("Integral is not elementary")
        else:
            print(e)
            raise e
        out1 = None
        out2 = None
    
    if out1==None or out2==None: return None
    return out1+out2

def IntegratePolynomialPart(poly):
    if poly==0:
        return Int.Integral([],[])
    if poly.fieldTower.towerHeight!=0:
        raise Exception()
    intPoly = Pol.Polynomial()
    for i in range(poly.degree+1):
        coeff = poly.getCoefficient(i)/(i+1)
        intPoly.setCoefficient(coeff, i+1)
    return Int.Integral(poly_rationals=[intPoly])
    #return Int.Integral(str(intPoly))

def IntegrateRationalFunction(func): # field=0, func element C(x)
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
    if func.denominator.isSquareFree():
        if func.denominator.degree>2:
            poly = func.denominator
            coeff_rat = func.numerator/func.denominator.differentiate()
            coeff_str = coeff_rat.strCustomVar("y")
            rootSum = RS.RootSum(poly,"({}){}".format(coeff_str, logExpression("x-y")),exprVar="y")
            return Int.Integral([],[],[rootSum])
        
        elif func.denominator.degree==1:
            arg = func.denominator
            coeff = func.numerator/func.denominator.getLeadingCoefficient()
            logExpr = Int.LogFunction(arg,coeff)
            return Int.Integral([],[logExpr])
            #return Int.Integral("{}*{}".format(str(coeff), logExpression(arg)))
        elif func.denominator.degree==2:
            f = func.numerator
            g = func.denominator
            a = g.getLeadingCoefficient()
            p = g.getCoefficient(1)/a
            q = g.getCoefficient(0)/a
            disc = sqrt((p/2)**2-q)
            if disc==0:
                raise Exception("denom was tested square free, but discriminant is 0")
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
        return HermiteReduction(func)
    
def HermiteReduction(rational):
    sqrFreeFactorization = rational.denominator.factorSquareFree()
    partialFractions = rational.PartialFraction(sqrFreeFactorization)
    integral = Int.Integral()
    for frac in partialFractions:
        j = frac[2]
        q_i = frac[1]
        r_ij = frac[0]
        if j>1:
            (s,t) = Pol.extendedEuclidGenF(q_i, q_i.differentiate(), r_ij) # s*q_i+t*q_i'=r_ij
            tPrime = 0 if isNumber(t) else t.differentiate()
            p1 = Int.Integral(poly_rationals=[(-1)*t/(j-1)/(q_i**(j-1))])
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
    
def IntegratePolynomialPartLogExt(poly, fieldTower):
    if poly == 0 or poly.isZero():
        return Int.Integral()
    if poly.deg0():
        return Integrate(poly.getCoefficient(0), fieldTower.prevTower())
    p = []
    for c in poly.getCoefficients():
        p.append(c)
    l = len(p)-1
    q = [0]*(len(p)+1)
    b = [0]*(len(p)+1)
    d = [0]*(len(p)+1)
    
    #PL = Integrate(p[l], p[l].getFieldTower())
    
    u = fieldTower.getLastExtension().characteristicFunction
    log_diff_u = u.differentiate()/u
    
    # 3 conditions
    #IntegratePolynomialPartLogExtCheckIntegralConditions(PL, fieldTower)
    """if PL==None: # integral of pl is elementary
        raise Int.IntegralNotElementaryException()
    logs = PL.getLogExpressionsInFieldTower(fieldTower.prevTower())
    if len(logs)>1:
        raise Int.IntegralNotElementaryException()
    if len(logs)!=0:
        if logs[1].argFunction == u:
            print(logs[1].argFunction, u)
            raise Int.IntegralNotElementaryException"""
    #logs = PL.getLogExpressionsInFieldTower(fieldTower.prevTower())
    #otherLogs = [log for log in PL.logExpressions if log not in logs]
    #cl = 0 if len(logs)==0 else logs[1].factor
    #d[l] = Int.Integral(poly_rationals=PL.poly_rational_partExpressions,logs=otherLogs, rootSums=PL.rootSums).asFunction() # Integral(pl)=cl*T+dl
    #b[l+1] = cl/(l+1)
    #q[l+1] = d[l+1]+b[l+1]
    
    for i in range(l,-1,-1):
        df = d[i+1]
       # if d[i+1]!=0:
        #    df = d[i+1].asFunction()
        integrand = p[i]+(-1)*(i+1)*df*log_diff_u # chechk if field towers are =
        Pi = Integrate(integrand,integrand.getFieldTower())
        IntegratePolynomialPartLogExtCheckIntegralConditions(Pi, fieldTower)
        prev = fieldTower.prevTower()
        logs = Pi.getLogExpressionsInFieldTower(prev)
        otherLogs = [log for log in Pi.logExpressions if log not in logs]
        ci = 0 if len(logs)==0 else logs[1].factor
        d[i] = Int.Integral(poly_rationals=Pi.poly_rational_partExpressions,logs=otherLogs,rootSums=Pi.rootSums).asFunction()
        b[i+1] = ci/(i+1)
        q[i+1] = d[i+1]+b[i+1]
    
    df = d[1]
 #   if d[1]!=0:
  #      df = d[1].asFunction()   
    integrand = p[0]+(-1)*df*log_diff_u
    P0 = Integrate(integrand,integrand.getFieldTower())
    if P0==None:
        raise Int.IntegralNotElementaryException()
    b1 = 0
    for log in P0.logExpressions:
        if log.argFunction==u:
            b1 = log.factor
            
    if b1==0:
        q_0 = P0.asFunction()
    else:
        q_0 = P0+Int.Integral(logs=[Int.LogFunction(u,-b1)]).asFunction()
    q[0] = q_0
    
    integralPoly = Pol.Polynomial(fieldTower=fieldTower)
    for i in range(l+1,-1,-1):
        integralPoly.setCoefficient(q[i], i)
        
    #print(q)
    #print(b)
    #print(d)
    return Int.Integral(poly_rationals=[integralPoly])

        
    
    #raise NotImplementedError
def IntegratePolynomialPartLogExtCheckIntegralConditions(integral,fieldTower):
    if integral==None: # integral of pl is elementary
        raise Int.IntegralNotElementaryException()
    logs = integral.getLogExpressionsInFieldTower(fieldTower.prevTower())
    if len(logs)>1:
        raise Int.IntegralNotElementaryException()
    if len(logs)!=0:
        if logs[1].argFunction == fieldTower.getLastExtension().characteristicFunction:
            print(logs[1].argFunction, fieldTower.getLastExtension().characteristicFunction)
            raise Int.IntegralNotElementaryException
def IntegrateRationalPartLogExt(func, fieldTower): # Hermite Recuction
    if func.numerator==0 or func.numerator.isZero():
        return Int.Integral()
    func.makeDenominatorMonic()
    sqrFreeFactorization = func.denominator.factorSquareFree()
    partialFractions = func.PartialFraction(sqrFreeFactorization)
    integral = Int.Integral()
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
                primitivePart = res.makeMonic()
                constantRoots = primitivePart.hasOnlyConstantCoefficients()
            if res==0 or res.isZero():
                constantRoots = False
                primitivePart = 0
                    
            #print("Only constant coefficients in primitive part {}: {}".format(primitivePart,constantRoots))
            
            if not constantRoots:
                raise Int.IntegralNotElementaryException("Integral is not elementary")
            
            roots = primitivePart.getRoots()
            vfunc = []
            for c in roots:
                vfunc.append(Pol.PolyGCD(a+bp*(-c), b))
            intPart = Int.Integral()
            for (c,v) in zip(roots,vfunc):
                logExpr = Int.LogFunction(v,c)
                intPart += Int.Integral(logs=[logExpr])
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

def printIntegral(func,fieldTower):
    integral = Integrate(func, fieldTower)
    
    return "Integral is not elementary" if integral==None else integral.printFull()
def integratetest(expected, got):
    return "integratetest: should be {}, got {}".format(expected, got)

if __name__ == '__main__':
    from Parse import parseField0PolyFromStr,parseField0RatFromStr
    FE.fieldTower = FE.FieldTower()
    
    polA = Pol.Polynomial([1])
    polB = Pol.Polynomial([2,-3,1])
    polC = Pol.Polynomial([-1,2])
    print(integratetest("0.5*log(x+(-0.5))",Integrate(Rat.RationalFunction(polA,polC),FE.FieldTower()))) #Integral of 1/(2x-1)

    ratA = Rat.RationalFunction(polA,polB) # 1/(x^2-3x+2)
    print(integratetest("log(x+(-2.0))+(-1.0)log(x+(-1.0))",Integrate(ratA, FE.FieldTower())))
    print("----")
    polD = parseField0PolyFromStr("3*x**2+x+1")
    polE = parseField0PolyFromStr("x**4+x**2+1")
    ratB = Rat.RationalFunction(polD,polE)

    print(integratetest("RootSum(w | w**4+w**2+1.0=0, ([0.75w**2+0.25w+0.25]/[w**3+0.5w])log(x-w))",Integrate(ratB,FE.FieldTower())))#print(Integrate(ratB,FE.FieldTower()))
    print("----")
    polA = parseField0PolyFromStr("2*x**4+x**2+x+(-1)")
    polB = parseField0PolyFromStr("(-1)+x")*(parseField0PolyFromStr("x+1")**3)#*parseField0PolyFromStr("x**2+x-1")**3
    ratA = Rat.RationalFunction(polA,polB)

    print(integratetest("2.0x+[4.25x]/[x+1.0]+[(-0.25)x]/[x**2+2.0x+1.0]+[(-0.25)x]/[x+1.0]+0.375*log(x+(-1.0))+(-4.375)log(x+1.0)",printIntegral(ratA, FE.FieldTower())))
    print("----")
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,1]),"T_1") # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    one = parseField0RatFromStr("1/1")
    ratX = parseField0RatFromStr("x/1")
    
    pol = Pol.Polynomial([0,1],fieldTower=FE.fieldTower)  # log(x)
    print(integratetest("(x)log(x)+((-1)x)", printIntegral(pol, FE.fieldTower)))
    # print(Integrate(pol, FE.fieldTower).printFull())
   
    pol = Pol.Polynomial([1,1,1],fieldTower=FE.fieldTower) # log(x)^2 + log(x) + 1
    print(integratetest("(x)log(x)**2+((-1.0)x)log(x)+(2.0x)",printIntegral(pol,FE.fieldTower)))
    #print(printIntegral(pol,FE.fieldTower))
    
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,ratX],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom,fieldTower=FE.fieldTower.getStrippedTower(1)) # 1/(x log(x))

    print(integratetest("log(log(x))", printIntegral(rat, FE.fieldTower)))
    
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([ratX,one],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom,fieldTower=FE.fieldTower.getStrippedTower(1)) # 1/(x + log(x))

    print(integratetest("Integral is not elementary", printIntegral(rat,FE.fieldTower)))
    
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,one,one],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom,fieldTower=FE.fieldTower.getStrippedTower(1)) # 1/(log(x)^2 + log(x))

    print(integratetest("Integral is not elementary", printIntegral(rat,FE.fieldTower)))
    
    
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,1],fieldTower=FE.FieldTower(fieldExtension1)),"T_2") # field extension with log(log(x))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1,fieldExtension2])
    numerator = Pol.Polynomial([1],fieldTower=FE.fieldTower)
    polX = Pol.Polynomial([parseField0PolyFromStr("x")],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,polX*Pol.Polynomial([0,1],fieldTower=FE.fieldTower.getStrippedTower(1))],fieldTower=FE.fieldTower)
    rat = Rat.RationalFunction(numerator,denom,fieldTower=FE.fieldTower) # 1/(x log(x) log(log(x)))
    print(integratetest("log(log(log(x)))", printIntegral(rat,FE.fieldTower)))
    
    
    pol = Pol.Polynomial([0,1],fieldTower=FE.fieldTower) # log(log(x)
   # print(pol.printFull())
    print(integratetest("Integral not elementary", printIntegral(pol, FE.fieldTower)))
  #  print(Integrate(pol, FE.fieldTower))
    ratX = Rat.RationalFunction(1,Pol.Polynomial([0,1]))
    pol = Pol.Polynomial([0,ratX],fieldTower=FE.fieldTower) # (1/x) log(log(x))
    print(pol)
    print(Integrate(pol,FE.fieldTower))
    
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,0,1]),"T") # field extension with log(x^2)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    rat = Rat.RationalFunction(2,Pol.Polynomial([0,parseField0PolyFromStr("x")],fieldTower=FE.fieldTower),fieldTower=FE.fieldTower)#1/(x log(x^2))
    #print(rat.printFull())
    print(integratetest("log(log(x**2))",printIntegral(rat,FE.fieldTower)))#printIntegral(rat, FE.fieldTower))
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,one],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom,fieldTower=FE.fieldTower.getStrippedTower(1)) # 1/(log(x^2))
    #print(rat.printFull())
    print(integratetest("Integral is not elementary", printIntegral(rat, FE.fieldTower)))
   # print(printIntegral(rat,FE.fieldTower))