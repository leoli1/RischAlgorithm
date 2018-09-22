'''
Created on 22.09.2018

@author: Leonard
'''
#from Polynomial import *
from FieldExtension import *
import Polynomial as Pol

class RationalFunction(object):
    '''
    classdocs
    '''


    def __init__(self, numerator, denominator, field = BASE_FIELD):
        '''
        Constructor
        '''
        self.field = field

        self.numerator = numerator
        self.denominator = denominator
        
        self.removeCommonFactors()
        
    def removeCommonFactors(self):
        '''
        cancels common factors, so that: numerator/denominator = p/q with gcd(p,q)=1
        '''
        gcd = Pol.PolyGCD(self.numerator, self.denominator)
        if not gcd.isConstant():
            self.numerator = Pol.PolyDiv(self.numerator, gcd)[0]
            self.denominator = Pol.PolyDiv(self.denominator,gcd)[0]
            #print(PolyDiv(self.denominator,gcd)[0])
    def isZero(self):
        return self.numerator.isZero()
    
    def Inverse(self):
        return RationalFunction(self.denominator,self.numerator,field=self.field)
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
        if _isPoly(other):
            return self.__mul__(RationalFunction(1,other,field=self.field))
        return self.__mul__(other.Inverse())
    def __str__(self):
        out = "["+str(self.numerator)+"]"
        d = str(self.denominator)
        if d=="(1)":
            return out
        return out+"/["+d+"]"

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