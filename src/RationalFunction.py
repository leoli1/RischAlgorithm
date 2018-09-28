'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
#from Polynomial import *
import FieldExtension as FE
import Polynomial as Pol
from Utils import isNumber

class RationalFunction(object):
    '''
    classdocs
    '''


    def __init__(self, numerator, denominator, field = FE.BASE_FIELD):
        '''
        Constructor
        '''
        self.field = field

        self.numerator = numerator
        self.denominator = denominator
        if not isNumber(numerator) and not isNumber(denominator):
            if self.numerator.field!=self.denominator.field:
                raise Exception()
        if not isNumber(numerator):
            if self.numerator.field!=self.field:
                raise Exception()
        if not isNumber(denominator):
            if self.denominator.field!=self.field:
                raise Exception()
        if self.denominator!=1 and self.numerator!=1 and self.numerator!=0:
            self.removeCommonFactors()
        
    def removeCommonFactors(self):
        '''
        cancels common factors, so that: numerator/denominator = p/q with gcd(p,q)=1
        '''
        gcd = Pol.PolyGCD(self.numerator, self.denominator)
        #print(gcd)
        if not gcd==1:
            self.numerator = Pol.PolyDiv(self.numerator, gcd)[0]
            self.denominator = Pol.PolyDiv(self.denominator,gcd)[0]
            #print(PolyDiv(self.denominator,gcd)[0])
    def isZero(self):
        return self.numerator.isZero()
    
    def differentiate(self):
        p = self.numerator
        q = self.denominator
        dp = p.differentiate()
        dq = q.differentiate()
        return RationalFunction(dp*q-p*dq,q*q,field=self.field) # (p/q)'=(p'q-pq')/(q^2)
    
    def PartialFraction(self, denomFactorization):
        
        raise NotImplementedError
    """def PartialFraction(self):
        if self.denominator.isSquareFree():
            if self.denominator.degree>2:
                raise NotImplementedError()
            else:
                raise NotImplementedError()
              #  c = 
        else:
            raise NotImplementedError()"""
    
    def Inverse(self):# f->1/f
        return RationalFunction(self.denominator,self.numerator,field=self.field)
    
    def MakeDenominatorMonic(self):
        lcoeff = self.denominator.getLeadingCoefficient()
        lcoeff_poly = Pol.Polynomial(coefficients=[lcoeff],field=self.field)
        self.numerator = self.numerator/lcoeff_poly
        self.denominator = self.denominator/lcoeff_poly
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other): # a/b+c/d = (ad+cb)/(bd)
        if other==0:
            return self
        if type(other) == Pol.Polynomial or _isPoly(other):
            return self.__add__(RationalFunction(other,1,field=self.field))
        num = self.numerator*other.denominator+self.denominator*other.numerator
        den = self.denominator*other.denominator
        return RationalFunction(num,den,field=self.field)
    def __mul__(self, other):
        if other == 1:
            return self
        if type(other)==Pol.Polynomial or _isPoly(other):
            return self.__mul__(RationalFunction(other,1,field=self.field))
        num = self.numerator * other.numerator
        denom = self.denominator*other.denominator
        return RationalFunction(num,denom,field=self.field)
    def __truediv__(self, other):
        if other==1:
            return self
        if type(other) == Pol.Polynomial or _isPoly(other):#_isPoly(other):
            return self.__mul__(RationalFunction(1,other,field=self.field))
        return self.__mul__(other.Inverse())
    def __str__(self):
        out = "["+str(self.numerator)+"]"
        if self.numerator==0:
            return out
        d = str(self.denominator)
        if d=="(1)":
            return out
        return out+"/["+d+"]"
    def printFull(self):
        numStr = str(self.numerator) if isNumber(self.numerator) else self.numerator.printFull()
        denomStr = str(self.denominator) if isNumber(self.denominator) else self.denominator.printFull()
        return "[{}]/[{}]".format(numStr,denomStr)

def _isPoly(x):
    try:
        x.degree
        return True
    except:
        return False
if __name__=='__main__':
    polA = Pol.Polynomial(coefficients=[1,2,1])
    polB = Pol.Polynomial(coefficients=[2,1])
    ratA = RationalFunction(polA,polB)
    print(ratA)
    polC = Pol.Polynomial(coefficients=[-1,0,1])
    polD = Pol.Polynomial(coefficients=[1,1])
    ratB = RationalFunction(polC,polD)
    print(ratB)
    fieldExtension1 = FE.FieldExtension(FE.TRANS_EXP,Pol.Polynomial([0,1]),"T_{{1}}") # field extension with e^x=exp(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    
    polE = Pol.Polynomial(coefficients=[polA,polC],field=1)
    polF = Pol.Polynomial(coefficients=[polB,polD],field=1)
    ratC = RationalFunction(polE,polF,field=1)
    print(ratC)
    ratC.MakeDenominatorMonic()
    print(ratC)
    
    print("[{}]'={}".format(ratA,ratA.differentiate()))
    polE = Pol.Polynomial(coefficients=[])