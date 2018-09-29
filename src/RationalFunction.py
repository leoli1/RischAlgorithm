'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
#from Polynomial import *
import FieldExtension as FE
import Polynomial as Pol
from Utils import isNumber
from test.test_binop import isnum

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
            if (isNumber(self.numerator) or not self.numerator.isZero()) and not ( isNumber(self.numerator) or isNumber(self.denominator)):
                self.removeCommonFactors()
        
        
        
        # TODO: clean up if's
        if (isNumber(numerator) or numerator.isConstant()) and (isNumber(denominator) or denominator.isConstant()):
            if isNumber(numerator) and isNumber(denominator):
                self.numerator = numerator/denominator
                print(self.numerator)
                self.denominator = 1
            elif isNumber(numerator) and denominator.isConstant():
                #print(str(numerator),str(denominator),"a")
                self.numerator = Pol.Polynomial([numerator/denominator.getConstant()],field=field)#(denominator/numerator).asRational().Inverse().numerator()
                self.denominator = 1
            elif isNumber(denominator) and numerator.isConstant():
                self.numerator = Pol.Polynomial([numerator.getConstant()/denominator],field=field)
                self.denominator = 1
            elif numerator.isConstant() and denominator.isConstant():
                self.numerator = Pol.Polynomial([numerator.getConstant()/denominator.getConstant()],field=field)
                self.denominator = 1
            
        if isNumber(self.denominator):
            self.numerator = self.numerator/self.denominator
            self.denominator = 1
        elif self.denominator.isConstant():
            self.numerator = self.numerator/self.denominator.getConstant()
            self.denominator = 1
            
        self.MakeDenominatorMonic()
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
        if isNumber(p):
            dp = 0
        else:
            dp = p.differentiate()
        if isNumber(q):
            dq = 0
        else:
            dq = q.differentiate()
        return RationalFunction(dp*q-p*dq,q*q,field=self.field) # (p/q)'=(p'q-pq')/(q^2)
    
    def PartialFraction(self, denomFactorization):
        
        self.MakeDenominatorMonic()
        return Pol.PartialFractionWithPowerFactorization(self.numerator, denomFactorization)
    
    def Inverse(self):# f->1/f
        return RationalFunction(self.denominator,self.numerator,field=self.field)
    
    def MakeDenominatorMonic(self):
        if isNumber(self.denominator):
            return
        lcoeff = self.denominator.getLeadingCoefficient()
        lcoeff_poly = Pol.Polynomial(coefficients=[lcoeff],field=self.field)
        if isNumber(self.numerator):
            self.numerator = Pol.Polynomial([self.numerator],field=self.field)/lcoeff_poly
        else:
            self.numerator = self.numerator/lcoeff_poly
        self.denominator = self.denominator/lcoeff_poly
        
    def isConstant(self):
        if isNumber(self.numerator):
            return isNumber(self.denominator) or self.denominator.isConstant()
        if self.numerator.isConstant():
            return isNumber(self.denominator) or self.denominator.isConstant()
        
    def getConstant(self):
        if not self.isConstant():
            return None
        if isNumber(self.numerator):
            if isNumber(self.denominator):
                return self.numerator/self.denominator
            else:
                return self.numerator/self.denominator.getConstant()
        else:
            if isNumber(self.denominator):
                return self.numerator.getConstant()/self.denominator
            else:
                return self.numerator.getConstant()/self.denominator.getConstant()
            
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
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __mul__(self, other):
        if other == 1:
            return self
        if type(other)==Pol.Polynomial or _isPoly(other):
            return self.__mul__(RationalFunction(other,1,field=self.field))
        if isNumber(other):
            num = self.numerator*other
            return RationalFunction(num,self.denominator, field=self.field)
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
        if self.numerator!=1 and not self.numerator.isConstant():
            out = "["+str(self.numerator)+"]"
        else:
            out = str(self.numerator.getLeadingCoefficient())
        if self.numerator==0:
            return "0"
        d = str(self.denominator)
        if d=="(1)" or self.denominator==1:
            return out
        return out+"/["+d+"]"
    
    def __repr__(self):
        return self.__str__()
    
    def strCustomVar(self, variable):
        out = "["+str(self.numerator.strCustomVar(variable))+"]"
        if self.numerator==0:
            return out
        d = str(self.denominator.strCustomVar(variable))
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
    
    polA = Pol.Polynomial([4])
    polB = Pol.Polynomial([3])
    ratA = RationalFunction(polA,polB)
    print(ratA)