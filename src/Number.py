'''
Created on 03.10.2018

@author: Leonard
'''
from __future__ import division
import Utils


class Rational(object):
    """
    Represents a rational Number p/q el Q
    """
    def __init__(self, p,q):
        self._p = p
        self._q = q
        if self._q==0:
            raise ValueError("Denominator can't be zero.")
        if self._q<0:
            self._p *= -1
            self._q *= -1
        
        self.removeCommonFactors()    
        self._p = int(self._p)
        self._q = int(self._q)
        
    @staticmethod
    def fromFloat(f):
        n = 0
        p = f
        while not Utils.isInt(p):
            p = p*10
            n+=1
        return Rational(int(p),10**n)
            
    def removeCommonFactors(self):
        gcd = NumberGCD(abs(self._p),abs(self._q))
        if gcd!=1:
            self._p /= gcd
            self._q /= gcd
        
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other==0:
            return self
        if type(other)==int or type(other)==float:
            return self.__add__(Rational.fromFloat(other))
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._q+self._q*other._p,self._q*other._q)
    def __rsub__(self, other):
        print("asdf")
    def __sub__(self, other):
        if type(other)!=Rational:
            return other.__add__(self)
        return self.__add__(-other)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        if type(other)==int or type(other)==float:
            return self.__mul__(Rational.fromFloat(other))
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._p,self._q*other._q)
    def __rtruediv__(self, other):
        return other*Rational(self._q,self._p)
    def __truediv__(self, other):
        if type(other)==int or type(other)==float:
            return self.__truediv__(Rational.fromFloat(other))
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._q,self._q*other._p)
    def __neg__(self):
        return Rational(-self._p,self._q)
    
    def __lt__(self, other):
        return float(self) < float(other)
    def __gt__(self,other):
        return float(self)>float(other)
    def __eq__(self, other):
        if other==None:
            return False
        return Utils.numberIsZero(self+(-other))
    def __ne__(self, other):
        return not self.__eq__(other)
    def __float__(self):
        return self._p/self._q
    def __int__(self):
        return self._p//self._q
    def __abs__(self):
        return Rational(abs(self._p),self._q)
    def __repr__(self):
        return str(self)
    def __str__(self):
        pstr = "({})".format(str(self._p)) if self._p<0 else str(self._p)
        if self._q==1:
            return pstr
        return pstr+"/"+str(self._q)
    
def NumberGCD(a,b):
    (A,B) = (a,b) if a>=b else (b,a)
    if B==0:
        return A
    r = A%B  # A=s*B+r
    if r==0:
        return B
    else:
        return NumberGCD(B,r)
    

if __name__=='__main__':
    a = Rational(12,8)
    b = Rational(5,4)
    print("{}+{}={}".format(a,b,a+b))
    c = Rational.fromFloat(1.25)
    print(c)
    