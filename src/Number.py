'''
Created on 03.10.2018

@author: Leonard
'''
from __future__ import division


class AlgebraicNumber(object):
    def __init__(self, poly):
        self.poly = poly
        if not Utils.isPoly(self.poly):
            raise TypeError("poly argument is no Polynomial")
        if self.poly.fieldTower.towerHeight==0:
            raise TypeError("poly field tower must be C(x)")
        

class SqrRootPolynomial(object):
    """
    Represents a number/expression a+b*sqrt(q) el Q[sqrt(q)] with q el Q
    """
    
    def __init__(self, radicand,a=0,b=0):
        self.radicand = radicand
        if not type(radicand) is Rational:
            self.radicand = Rational.fromFloat(radicand)
        self.a = ZERO if a==0 else a
        self.b = ZERO if b==0 else b # coefficient of sqrt(radicand)
        
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if type(other)==Rational:
            return SqrRootPolynomial(self.radicand, a=self.a+other, b=self.b)
        if not type(other) is SqrRootPolynomial:
            return self.__add__(Rational.fromFloat(other))
        if self.radicand!=other.radicand:
            raise ValueError()
        return SqrRootPolynomial(self.radicand,a=self.a+other.a,b=self.b+other.b)
    def __sub__(self, other):
        return self.__add__(-other)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other): # (a+bX)(p+qX) = (ap+bqX**2)+(bp+qa)X
        if type(other)==Rational:
            return SqrRootPolynomial(self.radicand, a=self.a*other, b=self.b*other)
        if not type(other) is SqrRootPolynomial:
            if Utils.isNumber(other):
                return self.__mul__(Rational.fromFloat(other))
            else:
                return other.__mul__(self)
        if self.radicand!=other.radicand:
            raise ValueError()
        return SqrRootPolynomial(self.radicand,a=self.a*other.a+self.b*other.b*self.radicand,b=self.b*other.a+self.a*other.b)
    def __rtruediv__(self, other):#other/self
        return self.Inverse().__mul__(other)
    def __truediv__(self, other): # self/other = (a+bX)/(p+qX) = [(a+bX)(p-qX)]/[(p+qX)(p-qX)] = (a+bX)(p-qX)/(p^2-q^2X^2)
        if type(other)==int or type(other)==float:
            return self * (1/other)
        if type(other)==Rational:
            return SqrRootPolynomial(self.radicand, a=self.a/other, b=self.b/other)
        if self.radicand!=other.radicand:
            raise ValueError()
        return self*other.Inverse()
    def __pow__(self, other):
        if not Utils.isInt(other):
            return TypeError()
        if other==0:
            return Rational(1,1)
        if self==0:
            return 0
        if other>0:
            return self*(self**(other-1))
        else:
            return 1/self * (self**(other+1))
    
    def Inverse(self):
        return self.conjugate()*(self.a**2-self.b**2*self.radicand).Inverse()
    def __neg__(self):
        return SqrRootPolynomial(self.radicand,a=-self.a,b=-self.b)
    def conjugate(self):
        return SqrRootPolynomial(self.radicand,a=self.a,b=-self.b)
    
    def __eq__(self, other):
        if other==None:
            return False
        if type(other)==Rational:
            return self.__eq__(SqrRootPolynomial(self.radicand,a=other))
        elif type(other)!=SqrRootPolynomial:
            return self.__eq__(Rational.fromFloat(other))
        return self.a==other.a and self.b==other.b
    def __ne__(self, other):
        return not self.__eq__(other)
    def __str__(self):
        rstr = str(self.radicand)
        if self.radicand._q==1:
            rstr = rstr.replace("(", "").replace(")","")
        if self.b==0:
            bStr = ""
        elif self.b==1:
            bStr= "sqrt({})".format(rstr)
        else:
            #if self.b._q==1:
            bStr= "{}sqrt({})".format(str(self.b),rstr)
            #else: 
            #    bStr= "({})sqrt({})".format(str(self.b),rstr)
        if self.a==0:
           return bStr
        else:
            if bStr=="":
                return str(self.a)
            else:
                return str(self.a)+"+"+bStr
    def __repr__(self):
        return self.__str__()
class Rational(AlgebraicNumber):
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
        if type(f) is Rational:
            return f
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
            return other.__add__(self)
        return Rational(self._p*other._q+self._q*other._p,self._q*other._q)
    def __rsub__(self, other): # other -self
        if other==0:
            return -self
        raise NotImplementedError()
        
    def __sub__(self, other):
        if type(other)==int or type(other)==float:
            return self.__add__(Rational.fromFloat(-other))
        #if type(other)!=Rational:
        #    return other.__add__(self)
        return self.__add__(-other)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        if type(other)==int or type(other)==float:
            return self.__mul__(Rational.fromFloat(other))
        if type(other)!=Rational:
            return other.__mul__(self)
        return Rational(self._p*other._p,self._q*other._q)
    def __rtruediv__(self, other):# other/self
        return other*self.Inverse()
    def __truediv__(self, other): # self/other
        if type(other)==int or type(other)==float:
            return self.__truediv__(Rational.fromFloat(other))
        if type(other) is SqrRootPolynomial:
            return self*other.Inverse()
        return Rational(self._p*other._q,self._q*other._p)
    def __neg__(self):
        return Rational(-self._p,self._q) 
    def Inverse(self):
        return Rational(self._q,self._p)
    def __pow__(self, other):
        if not Utils.isInt(other):
            return TypeError()
        if other==0:
            return Rational(1,1)
        if self==0:
            return 0
        if other>0:
            return self*(self**(other-1))
        else:
            return 1/self * (self**(other+1))
    def __lt__(self, other):
        return float(self) < float(other)
    def __gt__(self,other):
        return float(self)>float(other)
    def __eq__(self, other):
        if other==None:
            return False
        if type(other) is SqrRootPolynomial:
            return SqrRootPolynomial.__eq__(self)
        if type(other)!=Rational:
            return self.__eq__(Rational.fromFloat(other))
        return self._p==other._p and (self._q==other._q or self._p==0)
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
def NumberGCDList(l):
    if len(l)==2:
        return NumberGCD(l[0], l[1])
    return NumberGCD(l[0],NumberGCDList(l[1:len(l)]))

def NumberLCM(a,b):
    return a*b/NumberGCD(a, b)
def NumberLCMList(l):
    if len(l)==2:
        return NumberLCM(l[0], l[1])
    return NumberLCM(l[0],NumberLCMList(l[1:len(l)]))

ONE = Rational(1,1)
ZERO = Rational(0,1)

def sqrt(x):
    if type(x)!=Rational:
        return sqrt(Rational.fromFloat(x))
    if x._p<0:
        ai = False
    else:
        a = x._p**0.5
        ai = a.is_integer()
    b = x._q**0.5
    bi = b.is_integer()
    if ai and bi:
        return Rational(a,b)
    if ai and not bi: # sqrt(x)=a/sqrt(x._q)
       f = Rational(a,1)
       r = Rational(1,x._q) 
    elif not ai and bi:
        f = Rational(1,b)
        r = Rational(x._p,1)
    else:
        f = ONE
        r = x
    return SqrRootPolynomial(r,0,f)

def getAllDivisors(n):
    r = int(n**0.5)
    step = 2 if n%2 else 1
    a = []
    b = []
    for i in xrange(1,r+1,step):
        if n%i==0:
            a.append(i)
            if i!=r:
                b.append(n//i)
    return a+list(reversed(b))

import Utils

if __name__=='__main__':
    a = Rational(12,8)
    b = Rational(5,4)
    print("{}+{}={}".format(a,b,a+b))
    c = Rational.fromFloat(1.25)
    print(c)
    print(sqrt(4))
    print(5+2*sqrt(5))
    print(sqrt(4/5))
    print((1+2*sqrt(-3))/sqrt(-3))
    print(getAllDivisors(120))
    print(NumberLCMList([2,3,4]))