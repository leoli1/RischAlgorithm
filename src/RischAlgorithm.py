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
from Utils import *

from Matrix import Resultant
from Number import sqrt,Rational
import RischEquation
from IntegrateRational import *
import ExtendedPolynomial as ExtPol


logExtensionsInIntegral = []

def Integrate(func):
    if isNumber(func):
        func = Pol.Polynomial([func])
    else:
        func = func.reduceToLowestPossibleFieldTower()
    if isPoly(func):
        quot = func
        rat = 0
    else:
        (quot, rem) = Pol.PolyDiv(func.numerator, func.denominator)
        rat = Rat.RationalFunction(rem,func.denominator)
    
    if func!=0:
        Log("---------------\nBegin: Integrate {}.".format(str(func)))
        Log("Split up in polynomial/rational part: [{}] + [{}]".format(quot,rat))
    fieldTower = func.fieldTower
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
    except Int.IntegralNotElementaryError as e:
        if type(e) == Int.IntegralNotElementaryError:
            return None#Int.Integral("Integral is not elementary")
        else:
            print(e)
            raise e
        out1 = None
        out2 = None
    
    if out1==None or out2==None: return None
    integral = out1+out2
    if func !=0:
        Log("Finished integrating {} : {}\n---------------".format(str(func),str(integral)))
    return integral
    
def IntegratePolynomialPartLogExt(poly, fieldTower):
    """
    fieldTower = C(x,T_1,...T_N)
    poly el C(x,T_1,...T_(N-1))[T_N], TN is a logarithm over C(x,T_1,...T_(N-1))
    poly = p_l * T_N^l + p_(l-1) * T_N^(l-1) + ... + p_0
    """
    #print(poly)
    #print(poly.degree)
    #print(poly.fieldTower)
    if poly == 0 or poly.isZero():
        return Int.Integral()
    red = poly.reduceToLowestPossibleFieldTower()
    if red.fieldTower.towerHeight<fieldTower.towerHeight:
        return Integrate(red)
    p = []
    for c in poly.getCoefficients():
        p.append(c)
    l = len(p)-1
    q = [0]*(len(p)+1)
    b = [0]*(len(p)+1)
    d = [0]*(len(p)+1)
    
    
    u = fieldTower.getLastExtension().argFunction # last fieldExtension with T=log(u)
    log_diff_u = u.logDifferentiate()#u.differentiate()/u # u'/u
    
    one = Rational(1,1)
    zero = Rational(0,1)
    
    # p_l *T^l + p_(l-1)*T^(l-1) + ... + p_0 = (q_(l+1)*T^(l+1) + q_l*T^l + ... + q_0)' + (c_1 *v_1'/v_i + ... + c_m *v_m'/v_m)
    # (q_i*T^i)' = q_i'*T^i + q_i*i*T'*T^(i-1)
    # -> 0 = q_(l+1)'
    #  p_l = (l+1)*q_(l+1)*T' + q_l'
    # ...
    #  p_0 = q_1*T' + (qe_0)', qe_0 = q_0+(c_1 *log(v_1) + ... + c_m *log(v_m))
    for i in range(l,-1,-1):
        # q_(i+1) = d_(i+1) + b_(i+1)
        # integral(p_i-(i+1)*d_(i+1)*T') = l*b_(i+1)*T' + q_i
        #i = Rational.fromFloat(i)
        integrand = p[i]+(-one)*(i+1)*d[i+1]*log_diff_u
        int_reduced = integrand.reduceToLowestPossibleFieldTower()
        if int_reduced!=None:
            integrand = int_reduced
        P_i = Integrate(integrand)
        if i>0:
            IntegratePolynomialPartLogExtCheckIntegralConditions(P_i, fieldTower)
        else:
            if P_i==None:
                raise Int.IntegralNotElementaryError()
            
        prev = fieldTower.prevTower()
        #logs = P_i.getNewLogExpressionsInFieldTower(prev,fieldTower)
        logTerm = P_i.getLogExpression(u)
        otherLogs = [log for log in P_i.logExpressions if not log==logTerm]#in logs]
        ci = zero if logTerm==None else logTerm.factor # integral = P_i = c_i*T+d_i, d_i = P_i\logs
        d[i] = Int.Integral(poly_rationals=P_i.poly_rational_partExpressions,logs=otherLogs,rootSums=P_i.rootSums).asFunction()
        b[i+1] = ci/(i+1)
        q[i+1] = d[i+1]+b[i+1]
    
    #integrand = p[0]+(-1)*d[1]*log_diff_u
    #P0 = Integrate(integrand)
    #if P0==None:
    #    raise Int.IntegralNotElementaryError()
    #b1 = 0
    #for log in P0.logExpressions:
    #    if log.argFunction==u:
    #        b1 = log.factor
            
    if b[1]==0:
        q_0 = P_i.asFunction()
    else:
        q_0 = (P_i+Int.Integral(logs=[Int.LogFunction(u,(-1)*b[1])])).asFunction()
    q[0] = q_0
    
    integralPoly = Pol.Polynomial(fieldTower=fieldTower)
    for i in range(l+1,-1,-1):
        integralPoly.setCoefficient(q[i], i,callUpdates=False)
    integralPoly.updateCoefficientsAll()
        
    return Int.Integral(poly_rationals=[integralPoly])

        
def IntegratePolynomialPartLogExtCheckIntegralConditions(integral,fieldTower):
    if integral==None: # integral of pl is elementary
        raise Int.IntegralNotElementaryError()
    logs = integral.getNewLogExpressionsInFieldTower(fieldTower.prevTower(),fieldTower)
    if len(logs)>1: # at most one log extension of C(x,T_1,...,T_(N-1))
        if len(logs)==2:
            x = logs[0].factor if isNumber(logs[0].factor) else logs[0].factor.getConstant()
            y = logs[1].factor if isNumber(logs[1].factor) else logs[1].factor.getConstant()
            if abs(x)==abs(y):
                a = logs[0].argFunction if sign(x)==1 else logs[0].argFunction.Inverse()
                b = logs[1].argFunction if sign(y)==1 else logs[1].argFunction.Inverse()
                logs = [Int.LogFunction(a*b,abs(x))]
                integral.logExpressions = logs
        else:
            raise NotImplementedError()#Int.IntegralNotElementaryError()
    if len(logs)!=0: # if a log extension appears, it must be T_N
        if logs[0].argFunction != fieldTower.getLastExtension().argFunction:
            raise Int.IntegralNotElementaryError()
        
def IntegrateRationalPartLogExt(func, fieldTower): # Hermite Reduction
    if func.numerator==0 or func.numerator.isZero():
        return Int.Integral()
    red = func.reduceToLowestPossibleFieldTower()
    if red.fieldTower.towerHeight<fieldTower.towerHeight:
        return Integrate(red)
    func = func.makeDenominatorMonic()
    sqrFreeFactorization = func.denominator.factorSquareFree()
    partialFractions = func.BasicPartialFraction(sqrFreeFactorization)
    integral = Int.Integral()
    integratedPart = 0
    toIntegratePart = 0 # squarefree part
    #0.5/(x^2+x+0.25) *(1/log(x+0.5)) + -x/(x^2+x+0.25) * 1/(log(x+0.5)^2)
    for frac in partialFractions:
        j = frac[2]
        q_i = frac[1]
        r_ij = frac[0]
        integrand = r_ij/(q_i**j)
        while j>1:
            j = Rational.fromFloat(j)
            (s,t) = Pol.extendedEuclidGenF(q_i, q_i.differentiate(), r_ij)
            #p1 = Int.Integral(poly_rationals=[(-1)*t/(j-1)/(q_i**(j-1))])
            integratedPart += (-1)*t*(1/(j-1))/(q_i**(j-1))
            tPrime = 0 if isNumber(t) else t.differentiate()
            num = s+tPrime*(1/(j-1))
            integrand = num/(q_i**(j-1))
            r_ij = num
            j -= 1
        toIntegratePart += integrand
    
    integral += Int.Integral(poly_rationals=[integratedPart])
    if toIntegratePart==0:
        return integral
    
    r = toIntegratePart.reduceToLowestPossibleFieldTower()
    if r.fieldTower!=toIntegratePart.fieldTower:
        integral += Integrate(r)
    else:   
        a = toIntegratePart.numerator
        b = toIntegratePart.denominator
        bp = b.differentiate()
        
        zvar = FE.Variable('z')
        zExtension = FE.FieldExtension(FE.TRANSCENDENTAL_SYMBOL,1,zvar)
        newFieldTower = fieldTower.getStrippedTower(fieldTower.towerHeight-1)
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
            raise Int.IntegralNotElementaryError("Integral is not elementary")
        
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


def IntegratePolynomialPartExpExt(func, fieldTower):
    if func == 0 or func.isZero():
        return Int.Integral()
    
    #red = func.reduceToLowestPossibleFieldTower()
    #if red.fieldTower.towerHeight<fieldTower.towerHeight:
    #    return Integrate(red)
    #q = [0]*(func.degree+1)
    up = fieldTower.getLastExtension().argFunction.differentiate().reduceToLowestPossibleFieldTower()
    integralPoly = ExtPol.ExtendedPolynomial(fieldTower=fieldTower)
    for i in range(1,func.degree+1):
        o = func.getCoefficient(i).reduceToLowestPossibleFieldTower()
        ft = o.fieldTower if o.fieldTower.towerHeight>up.fieldTower.towerHeight else up.fieldTower
        s = RischEquation.solveRischEquation(i*up, o, ft)
        if s==None:
            raise Int.IntegralNotElementaryError()
        integralPoly.setCoefficient(s, i,callUpdates=False)
    if type(func) is ExtPol.ExtendedPolynomial:
        for i in range(-func.princDegree,0):
            o = func.getCoefficient(i).reduceToLowestPossibleFieldTower()
            ft = o.fieldTower if o.fieldTower.towerHeight>up.fieldTower.towerHeight else up.fieldTower
            s = RischEquation.solveRischEquation(i*up, o, ft)
            if s==None:
                raise Int.IntegralNotElementaryError()
            integralPoly.setCoefficient(s, i,callUpdates=False)
    integralPoly.updateCoefficientsAll()
    p = integralPoly.asPolynomial()
    if p!=None:
        integralPoly = p
    q0 = Integrate(func.getCoefficient(0))
    integral = Int.Integral(poly_rationals=[integralPoly])+q0
    return integral
        
    #raise NotImplementedError()
def IntegrateRationalPartExpExt(func, fieldTower):
    if func==0 or func.numerator==0 or func.numerator.isZero():
        return Int.Integral()

    func = func.makeDenominatorMonic()
    l = func.denominator.lowestDegree
    if l!=0:
        powerPol = Pol.Polynomial([0,Number.Rational(1,1)],fieldTower=fieldTower)**l # T^l
        qs = func.denominator/powerPol
        (rs,w) = Pol.extendedEuclidGenF(powerPol, qs, func.numerator)
        w_ext = ExtPol.extPolyFromPol_power(w, l) # w/T^l
        #rat = rs/qs
        #if type(rat) is Pol.Polynomial:
        #    ratInt = Integrate()
        return IntegratePolynomialPartExpExt(w_ext,fieldTower) + Integrate(rs/qs)
        #raise NotImplementedError()
    
    sqrFreeFactorization = func.denominator.factorSquareFree()
    partialFractions = func.BasicPartialFraction(sqrFreeFactorization)
    integral = Int.Integral()
    integratedPart = 0
    toIntegratePart = 0 # squarefree part

    for frac in partialFractions:
        j = frac[2]
        q_i = frac[1]
        r_ij = frac[0]
        integrand = r_ij/(q_i**j)
        while j>1:
            j = Rational.fromFloat(j)
            (s,t) = Pol.extendedEuclidGenF(q_i, q_i.differentiate(), r_ij)
            #p1 = Int.Integral(poly_rationals=[(-1)*t/(j-1)/(q_i**(j-1))])
            integratedPart += (-1)*t*(1/(j-1))/(q_i**(j-1))
            tPrime = 0 if isNumber(t) else t.differentiate()
            num = s+tPrime*(1/(j-1))
            integrand = num/(q_i**(j-1))
            r_ij = num
            j -= 1
        toIntegratePart += integrand
    
    integral += Int.Integral(poly_rationals=[integratedPart])
    if toIntegratePart==0:
        return integral
    
    r = toIntegratePart.reduceToLowestPossibleFieldTower()
    if r.fieldTower!=toIntegratePart.fieldTower:
        integral += Integrate(r)
    else:   
        a = toIntegratePart.numerator
        b = toIntegratePart.denominator
        bp = b.differentiate()
        
        zExtension = FE.FieldExtension(FE.TRANSCENDENTAL_SYMBOL,1,"z")
        newFieldTower = fieldTower.getStrippedTower(fieldTower.towerHeight-1)
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
            primitivePart = res.makeMonic().reduceToLowestPossibleFieldTower()
            constantRoots = primitivePart.hasOnlyConstantCoefficients()
        if res==0 or res.isZero():
            constantRoots = False # TODO
            primitivePart = 0
                
        #print("Only constant coefficients in primitive part {}: {}".format(primitivePart,constantRoots))
        
        if not constantRoots:
            raise Int.IntegralNotElementaryError("Integral is not elementary")
        
        
        roots = primitivePart.getRoots()
        vfunc = []
        for c in roots:
            x = a+bp*(-c)
            vfunc.append(Pol.PolyGCD(x, b))
        intPart = Int.Integral()
        for (c,v) in zip(roots,vfunc):
            logExpr = Int.LogFunction(v,c)
            intPart += Int.Integral(logs=[logExpr])
            vd = v.degree
            k = c*vd*(-1)
            polyPart = k*fieldTower.getLastExtension().argFunction
            intPart += Int.Integral(poly_rationals=[polyPart])
        integral += intPart
    
    return integral
    raise NotImplementedError()


def printIntegral(func,fieldTower=None):
    integral = Integrate(func)#, fieldTower)
    if integral!=None:
        integral.simplify()
        int_str = integral.printFull()
        if len(int_str)==0:
            return "C"
        else:
            return int_str + "  + C"
    else:
        return "Integral not elementary"
    
def integratetest(expected, got):
    return "integratetest: should be {}, got {}".format(expected, got)

if __name__ == '__main__':
    from Parse import parseField0PolyFromStr,parseField0RatFromStr,parseExpressionFromStr
    import Utils
    #global log_algorithm
    Utils.log_algorithm = False
    FE.fieldTower = FE.FieldTower()
    
    polA = Pol.Polynomial([1])
    polB = Pol.Polynomial([2,-3,1])
    polC = Pol.Polynomial([-1,2])
    print(integratetest("0.5*log(x+(-0.5))",Integrate(Rat.RationalFunction(polA,polC))))#,FE.FieldTower()))) #Integral of 1/(2x-1)

    ratA = Rat.RationalFunction(polA,polB) # 1/(x^2-3x+2)
    print(integratetest("log(x+(-2.0))+(-1.0)log(x+(-1.0))",Integrate(ratA)))
    print("----")
    polD = parseExpressionFromStr("3*x**2+x+1",FE.FieldTower())
    polE = parseExpressionFromStr("x**4+x**2+1",FE.FieldTower())
    ratB = Rat.RationalFunction(polD,polE)

    print(integratetest("RootSum(w | w**4+w**2+1.0=0, ([0.75w**2+0.25w+0.25]/[w**3+0.5w])log(x-w))",Integrate(ratB)))#print(Integrate(ratB,FE.FieldTower()))
    print("----")
    polA = parseField0PolyFromStr("2*x**4+x**2+x+(-1)")
    polB = parseField0PolyFromStr("(-1)+x")*(parseField0PolyFromStr("x+1")**3)#*parseField0PolyFromStr("x**2+x-1")**3
    ratA = Rat.RationalFunction(polA,polB)

    print(integratetest("2.0x+[4.25x]/[x+1.0]+[(-0.25)x]/[x**2+2.0x+1.0]+[(-0.25)x]/[x+1.0]+0.375*log(x+(-1.0))+(-4.375)log(x+1.0)",printIntegral(ratA)))
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
    rat = Rat.RationalFunction(numerator,denom) # 1/(x log(x))

    print(integratetest("log(log(x))", printIntegral(rat, FE.fieldTower)))
    
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([ratX,one],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom) # 1/(x + log(x))

    print(integratetest("Integral is not elementary", printIntegral(rat,FE.fieldTower)))
    
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,one,one],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom) # 1/(log(x)^2 + log(x))

    print(integratetest("Integral is not elementary", printIntegral(rat,FE.fieldTower)))
    
    c = Rat.RationalFunction(1,Pol.Polynomial([0,1])) # 1/x
    p = Pol.Polynomial([0,c],FE.fieldTower) # 1/x log(x)
    #print(p.printFull())
    print(integratetest("0.5log(x)**2",printIntegral(p,FE.fieldTower)))
    
    c = Rat.RationalFunction(1,Pol.Polynomial([1,1])) # 1/(x+1)
    p = Pol.Polynomial([0,c],FE.fieldTower) # 1/(x+1) log(x)
    #print(p.printFull())
    print(integratetest("Integral is not elementary", printIntegral(p,FE.fieldTower)))
    
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,1],fieldTower=FE.FieldTower(fieldExtension1)),"T_2") # field extension with log(log(x))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1,fieldExtension2])
    
    
    numerator = Pol.Polynomial([1],fieldTower=FE.fieldTower)
    polX = Pol.Polynomial([parseField0PolyFromStr("x")],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,polX*Pol.Polynomial([0,1],fieldTower=FE.fieldTower.getStrippedTower(1))],fieldTower=FE.fieldTower)
    rat = Rat.RationalFunction(numerator,denom) # 1/(x log(x) log(log(x)))
    print(integratetest("log(log(log(x)))", printIntegral(rat,FE.fieldTower)))
    
    
    pol = Pol.Polynomial([0,1],fieldTower=FE.fieldTower) # log(log(x)
    # print(pol.printFull())
    print(integratetest("Integral not elementary", printIntegral(pol, FE.fieldTower)))
    #  print(Integrate(pol, FE.fieldTower))
    ratX = Rat.RationalFunction(1,Pol.Polynomial([0,1]))
    pol = Pol.Polynomial([0,ratX],fieldTower=FE.fieldTower) # (1/x) log(log(x))
    print(integratetest("log(x)*log(log(x))+(-1)*log(x)",printIntegral(pol, FE.fieldTower)))
    
    X = Pol.Polynomial([0,1],fieldTower=FE.FieldTower())
    X2 = X*X
    log = Pol.Polynomial([0,1],fieldTower=FE.fieldTower.getStrippedTower(1))
    log2 = log*log
    log3 = log2*log
    loglogCoeff = 2*X2*log3+X2*log2+(-1)
    num = Pol.Polynomial([log2,loglogCoeff],fieldTower=FE.fieldTower)
    denom = Pol.Polynomial([0,X*log2],fieldTower=FE.fieldTower)
    rat = Rat.RationalFunction(num,denom)
    print(integratetest("Integral is not elementary", printIntegral(rat, FE.fieldTower)))
    
    loglogCoeff = 2*X2*log3+X2*log2+(-1)
    num = Pol.Polynomial([log,loglogCoeff],fieldTower=FE.fieldTower) # (2x^2*log(x)^3+x^2*log(x)^2-1) *log(log(x)) + log(x)
    denom = Pol.Polynomial([0,X*log2],fieldTower=FE.fieldTower) # x*log(x)^2*log(log(x))
    rat = Rat.RationalFunction(num,denom) 
    #print(rat)
    print(integratetest("(x**2)log(x)+[(-1.0)log(x)+1.0]/[log(x)]+log(log(log(x)))",printIntegral(rat,FE.fieldTower)))
    # print(printIntegral(rat,FE.fieldTower))
    
    
    
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,0,1]),"T") # field extension with log(x^2)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    rat = Rat.RationalFunction(2,Pol.Polynomial([0,parseField0PolyFromStr("x")],fieldTower=FE.fieldTower))#1/(x log(x^2))
    #print(rat.printFull())
    print(integratetest("log(log(x**2))",printIntegral(rat,FE.fieldTower)))#printIntegral(rat, FE.fieldTower))
    numerator = Pol.Polynomial([one],fieldTower=FE.fieldTower.getStrippedTower(1))
    denom = Pol.Polynomial([0,one],fieldTower=FE.fieldTower.getStrippedTower(1))
    rat = Rat.RationalFunction(numerator,denom) # 1/(log(x^2))
    #print(rat.printFull())
    print(integratetest("Integral is not elementary", printIntegral(rat, FE.fieldTower)))
    # print(printIntegral(rat,FE.fieldTower))