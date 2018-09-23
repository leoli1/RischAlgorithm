'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
import RationalFunction as Rat

numbers = [int,float,complex]

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
    def isConstant(self):
        self.updateDegree()
        return self.degree==0
    
    def differentiate(self):
        dPoly = Polynomial(field = self.field)
        for i in range(self.degree,-1,-1):
            p = self.getCoefficient(i)
            if self.field == FE.BASE_FIELD:
                if i>0:
                    dPoly.setCoefficient(p*i, i-1)
            elif FE.fieldExtension.extensionType==FE.TRANS_EXP: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'*T^i since T'=u'T
                dPoly += Monomial(i, p.differentiate(), field=FE.BASEFUNCTION_FIELD)
                dPoly += i*Monomial(i,p*FE.fieldExtension.characteristicFunction.differentiate(),field=FE.BASEFUNCTION_FIELD)
            elif FE.fieldExtension.extensionType==FE.TRANS_LOG:
                dPoly += Monomial(i,p.differentiate(),field=FE.BASEFUNCTION_FIELD)
                if i>0:
                    dPoly += i*Monomial(i-1,p*Rat.RationalFunction(FE.fieldExtension.characteristicFunction.differentiate(),FE.fieldExtension.characteristicFunction))
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
        if type(other) in numbers:
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
        (quot,rem)= PolyDiv(self, other)
        if rem==0:
            return quot
        
        return Rat.RationalFunction(rem,other)+quot
    def __str__(self):
        out = ""
        var = FE.VARIABLES[self.field]
        for i in range(self.degree,-1,-1):
            if self.coeffIsZero(i):
                continue
            out += "("+str(self.getCoefficient(i))+")"+("["+var+"**"+str(i)+"]" if i!=0 else "") +"+"
        out = out.strip("+")
        if len(out)==0:
            return "0"
        return out.strip("+")
def PolyDiv(polA, polB):
    if (polA.field!=polB.field):
        raise Exception("Polynomials have to have coefficients in the same field in order to apply PolyDiv")
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
    (s,r) = PolyDiv(A,B) # A=s*B+r
    if r == 0 or r.isZero():
        if B.isConstant() and type(B.getCoefficient(0)) in numbers:
            return B*(1/B.getCoefficient(0)) # normalize
        return B
    return PolyGCD(B,r)
    
def Monomial(degree, coeff, field = FE.BASE_FIELD):
    return Polynomial([0]*degree+[coeff], field=field)
if __name__ == '__main__': #tests
    print("---------- General poly tests")
    polA = Polynomial(coefficients=[1,0,2,4])
    print(polA)
    polB = Polynomial(coefficients=[0,1,2])
    print(polB)
    polC = polA+polB
    print("{} + {} = {}".format(polA,polB,polC))
    polD = polA*polB
    print("[{}] * [{}] = {}".format(polA,polB,polD))
    
    polF_T = Polynomial(coefficients=[polA,polB],field=FE.BASEFUNCTION_FIELD)
    print(polF_T)
    polG_T = polF_T*polF_T
    print("[{}] * [{}] = {}".format(polF_T,polF_T,polG_T))
    
    print("---------- Diff poly tests")
    print("[{}]'={}".format(polA,polA.differentiate()))
    
    print("---------- Poly div test:")
    polA = Polynomial(coefficients=[11,4,0,3])
    polB = Polynomial(coefficients=[2,-3,1])
    (quot,rem) = PolyDiv(polA, polB)
    print("[{}] / [{}] = {} R {}".format(polA,polB,quot,rem))
    
    polA = Polynomial(coefficients=[18,9,-11,12])
    polB = Polynomial(coefficients=[3,4])
    (quot,rem) = PolyDiv(polA, polB)
    print("[{}] / [{}] = {} R {}".format(polA,polB,quot,rem)) # (12x^3-11x^2+9x+18)/(4x+3)
    
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
    
        