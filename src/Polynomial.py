'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
import RationalFunction as Rat

from Utils import isNumber

class Polynomial(object):
    '''
    class for polynomials
    '''


    def __init__(self, coefficients = None, field=FE.BASE_FIELD):
        '''
        field: the field(-extension) of its coefficients
        '''
        self.field = field
        if coefficients == None:
            coefficients = [0]
        if type(coefficients)!=list:
            raise Exception()
        self._coefficients = coefficients
        
        self.degree = len(self._coefficients)-1
        self.updateDegree()

    def setCoefficients(self, newCoeffs):
        self._coefficients = newCoeffs
        self.updateDegree()
        
    def setCoefficient(self, coeff, power):
        if power>self.degree:
            self._coefficients += [0]*(power-self.degree)
            self.degree = power

        self._coefficients[power] = coeff
        self.updateDegree()
        
    def getCoefficient(self, power):
        return 0 if power>self.degree else self._coefficients[power]
    def getLeadingCoefficient(self):
        return self.getCoefficient(self.degree)
    def getFieldTower(self):
        return FE.fieldTower.getStrippedTower(self.field)
    def getFieldExtension(self):
        if self.field==FE.BASE_FIELD:
            return None
        t = self.getFieldTower()
        return t.getFieldExtension(t.towerHeight-1)
    
    def updateDegree(self):
        '''
        tests if highest power terms have coefficients equal to 0 and eventually recalculates the degree
        '''
        newDeg = 0
        for i in range(self.degree, 0, -1):
            if not self.coeffIsZero(i):
                newDeg = i
                break
        self.degree = newDeg
        self._coefficients = self._coefficients[0:newDeg+1]
        
    def coeffIsZero(self, power):
        c = self.getCoefficient(power)
        if c==0:
            return True
        elif self.field == FE.BASE_FIELD:
            return False 
        return c.isZero()
    def isZero(self):
        for i in range(self.degree+1):
            if not self.coeffIsZero(i):
                return False
        return True
    def deg0(self):
        self.updateDegree()
        return self.degree==0
    def isConstant(self):
        if self.deg0():
            if self.field==0:
                return True
            else:
                return self.getCoefficient(0).isConstant()
        else:
            return False
        
    def isSquareFree(self): # polynomial p is square-free iff gcd(p,p')=1
        gcd = PolyGCD(self,self.differentiateWRTtoPolyVar())#PolyGCD(self,self.differentiate())
        return gcd==1
    def factorSquareFree(self):
        """
        factors self = f = a_1*(a_2^2)*(a_3^3)...
        where the a_i are square-free and gcd(a_i,a_j)=1 for i!=j
        returns [(a_1,1),(a_2,2), ...]
        """
        #raise NotImplementedError()
        if self.isSquareFree():
            return [(self,1)]
        c = PolyGCD(self, self.differentiateWRTtoPolyVar()) # = a_2*(a_3)^2*(a_4)^3...
        w = self/c # = a_1*a_2*a_3...
        y = PolyGCD(c, w) # = a_2*a_3...
        a1 = w/y
        rest = c.factorSquareFree()
        an = [(a1,1)]+rest
        for i in range(1,len(an)):
            an[i] = (an[i][0],an[i][1]+1)
        return an
        
    
    def evaluate(self, val):
        if self.field!=FE.BASE_FIELD:
            raise Exception()
        sum = 0
        for i in range(self.degree+1):
            sum += self.getCoefficient(i)*(val**i)
        return sum
    def getVariable(self):
        return "x" if self.field==FE.BASE_FIELD else self.getFieldExtension().variable
    def differentiate(self):
        dPoly = Polynomial(field = self.field)
        for i in range(self.degree,-1,-1):
            p = self.getCoefficient(i)
            fieldExtension = self.getFieldExtension()
            if p==0:
                continue
            if self.field == FE.BASE_FIELD:
                if i>0:
                    dPoly.setCoefficient(p*i, i-1)
            elif fieldExtension.extensionType==FE.TRANS_EXP: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'*T^i since T'=u'T
                dPoly += Monomial(i, p.differentiate(), field=self.field)
                dPoly += i*Monomial(i,p*fieldExtension.characteristicFunction.differentiate(),field=self.field)
            elif fieldExtension.extensionType==FE.TRANS_LOG:
                dPoly += Monomial(i,p.differentiate(),field=self.field)
                if i>0:
                    log_diff_u = Rat.RationalFunction(fieldExtension.characteristicFunction.differentiate(),fieldExtension.characteristicFunction,field=self.field-1)
                    dPoly += i*Monomial(i-1,p*log_diff_u, field=self.field)
        return dPoly
    def differentiateWRTtoPolyVar(self):
        """
        differentiates function with respect to the field extension variable of the polynomial:
        p = 2*T**2+T+1
        -> differentiateWRTtoPolyVar(p) = (d/dT)p = 4T+1
        """
        dPoly = Polynomial(field=self.field)
        for i in range(self.degree,0,-1):
            p = self.getCoefficient(i)
            dPoly.setCoefficient(p*i, i-1)
        return dPoly
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other ==0:
            return self
        if type(other)==Rat.RationalFunction:
            return other.__add__(self)
        if (self.field!=other.field):
            raise Exception("Polynomials have to have coefficients in the same field in order to add them")
        tDeg = max(self.degree,other.degree)
        tPoly = Polynomial(field=self.field)
        for i in range(tDeg,-1,-1):
            a = self.getCoefficient(i)
            b = other.getCoefficient(i)
            tPoly.setCoefficient(a+b, i)
        return tPoly
    def __sub__(self, other):
        return self.__add__(other*(-1))
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        if other == 1:
            return self
        if isNumber(other):
            return self.__mul__(Polynomial(coefficients=[other],field=self.field))
        elif type(other) == Rat.RationalFunction:
            return other.__mul__(self)
        if (self.field!=other.field):
            raise Exception("Polynomials have to have coefficients in the same field in order to multiply them")
        tDeg = self.degree + other.degree
        tPoly = Polynomial(field=self.field)
        # calcs new coefficients using Cauchy-Product
        for i in range(tDeg,-1,-1):
            newCoeff = 0
            
            for j in range(0,i+1):
                a = self.getCoefficient(j)
                b = other.getCoefficient(i-j)
                newCoeff += a*b
            tPoly.setCoefficient(newCoeff, i)
        return tPoly
    
    def __truediv__(self, other):
        if isNumber(other):
            return self.__mul__(1/other)
        (quot,rem)= PolyDiv(self, other)
        if rem==0:
            return quot
        
        return Rat.RationalFunction(rem,other)+quot
    def __str__(self):
        return self.strCustomVar(self.getVariable())
    def strCustomVar(self, variable):
        out = ""#FE.VARIABLES[self.field]
        for i in range(self.degree,-1,-1):
            if self.coeffIsZero(i):
                continue
            coeff = ""
            coeff_v = self.getCoefficient(i)
            if coeff_v != 1 or i==0:
                if isNumber(coeff_v) and coeff_v>0:
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant():
                    coeff = str(coeff_v)
                else:
                    coeff = "({})".format(str(self.getCoefficient(i)))
                
            var = ""
            if i>1:
                var = "{}**{}".format(variable,i)
            elif i==1:
                var = variable
            out += coeff +"{}+".format(var)
        out = out.strip("+")
        if len(out)==0:
            return "0"
        return out.strip("+")
    def printFull(self):
        """
        writes out field extension variables
        """
        out = ""
        fE = self.getFieldExtension()
        for i in range(self.degree,-1,-1):
            if self.coeffIsZero(i):
                continue
            coeff = self.getCoefficient(i)
            coeffstr = str(coeff) if isNumber(coeff) else coeff.printFull()
            out += "({})".format(coeffstr)
            if i>0:
                if fE==None:
                    out += "x**{}".format(i)
                else:
                    out += fE.strFunc()
                out += "+"
        if len(out)==0:
            return "0"
        return out.strip("+")
    
def PolyDiv(polA, polB):
    if (polA.field!=polB.field):
        raise Exception("Polynomials have to have coefficients in the same field in order to apply PolyDiv")
    if polA.deg0() and polB.deg0():
        return (Polynomial([polA.getCoefficient(0)/polB.getCoefficient(0)],field=polA.field), 0)
    if (polA.degree<polB.degree):
        return (0,polA)
    power = polA.degree-polB.degree
    coeff = polA.getCoefficient(polA.degree)/polB.getCoefficient(polB.degree)
    monomial = Monomial(power,coeff,field=polB.field)
    newA = polA-(monomial*polB)
    if newA.isZero():
        return (monomial,0)
    else:
        (quot,remainder) = PolyDiv(newA,polB)
        return (monomial+quot,remainder)
    
def PolyGCD(polA,polB):
    (A,B) = (polA,polB) if polA.degree>=polB.degree else (polB,polA)
    if B.isZero():
        return A
    (s,r) = PolyDiv(A,B) # A=s*B+r
    if r == 0 or r.isZero():
        if B.deg0():# and type(B.getCoefficient(0)) in numbers:
            return 1#B*(1/B.getCoefficient(0)) # normalize
        return B
    gcd = PolyGCD(B,r)
    if gcd!=1:
        lcoeff = gcd.getLeadingCoefficient()
        if isNumber(lcoeff):
            return gcd*Polynomial([1/gcd.getLeadingCoefficient()],field=gcd.field)
        else:
            return gcd*Polynomial([gcd.getLeadingCoefficient().Inverse()],field=gcd.field)
    else:
        return 1
    
def extendedEuclid(p,q):
    """
    xp+yq=1
    """
    if (PolyGCD(p,q))!=1:
        raise Exception()
    (P,Q) = (p,q) if p.degree>=q.degree else (q,p)
    (s,r) = PolyDiv(P,Q)
    if r.deg0():
        if r.field==FE.BASE_FIELD:
            r_inv_poly = Polynomial([1/r.getCoefficient(0)])
            return (r_inv_poly,r_inv_poly*(-1)*s)
        else:
            r_inv_poly = Polynomial([r.getCoefficient(0).Inverse()],field=r.field)#Rat.RationalFunction(1,r,field=r.field)
            return (r_inv_poly,r_inv_poly*(-1)*s)
    (x,y) = extendedEuclid(Q, r)
    return (y,x+(-1)*s*y)
    
def extendedEuclidGenF(p,q,f):
    """
    xp+yq=f
    """
    if PolyGCD(p, q)!=1:
        raise Exception()
    (s,t) = extendedEuclid(p, q)
    
def PartialFraction(f,p,q):
    """
    f/(pq)=x/p+y/q
    gcd(p,q)=1
    """
    if f.field!=p.field or f.field!=q.field:
        raise Exception()
    if PolyGCD(p,q)!=1:
        raise Exception()
    
def Monomial(degree, coeff, field = FE.BASE_FIELD):
    return Polynomial([0]*degree+[coeff], field=field)

def printFactorization(fact):
    out = ""
    for i in range(len(fact)):
        if fact[i][0]!=1:
            out += "({})".format(fact[i][0])
            if fact[i][1]!=1:
                out += "**{}".format(fact[i][1])
            out += "*"
    return out.strip("*")
            

if __name__ == '__main__': #tests
    from Parse import parseField0PolyFromStr,parseField0RatFromStr
    fieldExtension1 = FE.FieldExtension(FE.TRANS_EXP,Polynomial([0,1]),"T_{{1}}") # field extension with e^x=exp(x)
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Polynomial([Polynomial([0,1])], field=1),"T_{{2}}") # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1,fieldExtension2])
    #FE.updateVariables()
    print("---------- General poly tests")
    polA = Polynomial(coefficients=[1,0,2,4])
    print(polA)
    polB = Polynomial(coefficients=[0,1,2])
    print(polB)
    polC = polA+polB
    print("{} + {} = {}".format(polA,polB,polC))
    polD = polA*polB
    print("[{}] * [{}] = {}".format(polA,polB,polD))
    
    polF_T = Polynomial(coefficients=[polA,polB],field=1) # polB*(T_1)**1+polA = polB*(e**x)+polA
    print(polF_T)
    polG_T = polF_T*polF_T
    print("[{}] * [{}] = {}".format(polF_T,polF_T,polG_T))
    
    polH = Polynomial(coefficients=[0,Polynomial([1])], field=1)
    polG = Polynomial(coefficients=[polH,polH],field=2)
    print(polH)
    print(polG)
    print(polG.printFull())
    
    print("---------- Diff poly tests")
    print("[{}]'={}".format(polA,polA.differentiate()))
    print("[{}]'={}".format(polG,polG.differentiate()))
    
    print("---------- Poly div test:")
    polA = Polynomial(coefficients=[11,4,0,3])
    polB = Polynomial(coefficients=[2,-3,1])
    (quot,rem) = PolyDiv(polA, polB)
    print("[{}] / [{}] = {} R {}".format(polA,polB,quot,rem))
    
    polA = Polynomial(coefficients=[18,9,-11,12])
    polB = Polynomial(coefficients=[3,4])
    (quot,rem) = PolyDiv(polA, polB)
    print("[{}] / [{}] = {} R {}".format(polA,polB,quot,rem)) # (12x^3-11x^2+9x+18)/(4x+3)
    
    polA = parseField0PolyFromStr("x**2+")
    
    (quot,rem) = PolyDiv(polG_T, polF_T)
    print("[{}] / [{}] = {} R {}".format(polG_T,polF_T,quot,rem))
    
    print("---------- GCD Tests")
    polA = Polynomial(coefficients=[-4,0,1])#x^2-4
    polB = Polynomial(coefficients=[2,1])#x+2
    print("gcd({},{})={}".format(polA,polB,PolyGCD(polA, polB)))
    
    polA = Polynomial(coefficients=[1,0,1])
    polB = Polynomial(coefficients=[-1,1])
    print("gcd({},{})={}".format(polA,polB,PolyGCD(polA, polB)))
    
    polA = Polynomial(coefficients=[1,1,2,0,1,0,0,1])#1+x+2x^2+x^4+x^7
    polB = Polynomial(coefficients=[1,1,1,1,1,1])#1+x+x^2+x^3+x^4+x^5
    print("gcd({},{})={}".format(polA,polB,PolyGCD(polA, polB)))
    
    print("---------- ExtendedEuclid Test")
    polA = parseField0PolyFromStr("x**3+x**2+(-1)*x+2")
    polB = parseField0PolyFromStr("x**2+1")
    (x,y) = extendedEuclid(polA,polB)
    print(str(x),str(y))
    
    print("---------- Factor square free test")
    polA = parseField0PolyFromStr("x**3+4*x**2+5*x+2")
    fact = polA.factorSquareFree()
    print("{} = {}".format(str(polA),printFactorization(fact)))
    
    
    # http://www.wolframalpha.com/input/?i=expand+(T-x)(T%2B(x%2B1)%2F(x%2B2))%5E2+at+T%3D0
    ratA = parseField0RatFromStr("x+2*x**2+x**3/4+4*x+x**2")
    ratB = parseField0RatFromStr("1+(-2)*x+(-5)*x**2+(-2)*x**3/4+4*x+x**2")
    ratC = parseField0RatFromStr("2+(-1)*x**2/2+x")
    ratD = parseField0RatFromStr("1/1")
    
    polA = Polynomial([ratA*(-1),ratB,ratC,ratD],field=1)
    print(polA)
    fact = polA.factorSquareFree()
    print("{} = {}".format(str(polA),printFactorization(fact)))
        