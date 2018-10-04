'''
Created on 03.10.2018

@author: Leonard
'''
from __future__ import division
from Utils import *


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
        
    @staticmethod
    def fromFloat(f):
        n = 0
        p = f
        while not isInt(p):
            p = p*10
            n+=1
        return Rational(int(p),10**n)
            
    def removeCommonFactors(self):
        gcd = NumberGCD(abs(self._p),abs(self._q))
        self._p /= gcd
        self._q /= gcd
        
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other==0:
            return self
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._q+self._q*other._p,self._q*other._q)
    def __rsub__(self, other):
        print("asdf")
    def __sub__(self, other):
        if type(other)!=Rational:
            return other.__add__(self)
        return self.__add__(-other)
    def __mul__(self, other):
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._p,self._q*other._q)
    def __truediv__(self, other):
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._q,self._q*other._p)
    def __neg__(self):
        return Rational(-self._p,self._q)
    
    def __abs__(self):
        return Rational(abs(self._p),self._q)
    def __repr__(self):
        return str(self)
    def __str__(self):
        if self._q==1:
            return str(self._p)
        return str(self._p)+"/"+str(self._q)
    
def NumberGCD(a,b):
    (A,B) = (a,b) if a>=b else (b,a)
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
    