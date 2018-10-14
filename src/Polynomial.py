'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
import RationalFunction as Rat

from Utils import *
from FieldExtension import fieldTower
from Parse import parseExpressionFromStr
import Number
import ExtendedPolynomial as ExtPol


def call_counter(func):
    def helper(*args, **kwargs):
        helper.calls += 1
        return func(*args, **kwargs)
    helper.calls = 0
    helper.__name__= func.__name__
    return helper

class Polynomial(object):
    def __init__(self, coefficients=None,variable=None):
        
        self.__derivative = None
        self.__logDerivative = None # u'/u
        self.__logGCD = None # GCD(u,u')
        
        
        if coefficients == None:
            self._coefficients = [Number.ZERO]
        else:
            self._coefficients = coefficients
            
        if variable==None:
            variable = FE.getVariable('x')
        self.variable = variable
        
        if type(variable)!=None and type(variable) != FE.Variable:
            raise TypeError()
        
        if coefficients!=None:
            self.updateCoefficientsAll()
            
        for c in self._coefficients:
            if not isNumber(c):
                if c.variable==self.variable:
                    raise TypeError()
                
        # self.replaceNumbersWithRationals()
    # ========================================== Coefficients stuff =========================================
    def getCoefficient(self, power):
        if power>self.degree:
            return Number.ZERO
        return self._coefficients[power]
    def setCoefficient(self, power, coeff, callUpdates=True):
        if type(power)!=int: raise TypeError()
        if power>self.degree:
            self._coefficients += [Number.ZERO]*(power-self.degree)
        self._coefficients[power] = coeff
        if callUpdates:
            self.updateCoefficientsAll()
    
    def getLeadingCoefficient(self):
        return self.getCoefficient(self.degree)
    
    def getCoefficients(self):
        return self._coefficients
    def coeffIsZero(self, power):
        return objEqualsNumber(self.getCoefficient(power),0)
    
    def updateCoefficientsAll(self):
        self.updateCoefficients()

    def updateCoefficients(self):
        """
        Removes leading coefficients equal to Zero
        """
        newDeg = 0
        for i in range(len(self._coefficients)-1,-1,-1):
            if not self.coeffIsZero(i):
                newDeg = i
                break
        self._coefficients = self._coefficients[0:newDeg+1]
        
    def simplifyCoefficients(self):
        pass#TODO #raise NotImplementedError()
    
    def simplified(self):
        p = self.reduceToLowestPossibleFieldTower()
        for i in range(self.degree+1):
            c = self.getCoefficient(i)
            if not isNumber(c):
                self._coefficients[i] = c.simplified()
                
        return p
    
    def isConstant(self):
        if self.deg0():
            if isNumber(self.getCoefficient(0)):
                return True
            else:
                return self.getCoefficient(0).isConstant()
        else:
            return False
            
    def getConstant(self):
        if self.deg0():
            if isNumber(self.getCoefficient(0)):
                return self.getCoefficient(0)
            else:
                return self.getCoefficient(0).getConstant()
        else:
            return None
        
    def isZero(self):
        return self.getConstant()==0
    
    def getConstantCoefficients(self):
        coeffs = []
        for c in self.getCoefficients():
            if isNumber(c):
                coeffs.append(c)
            else:
                cc = c.getConstant()
                if cc==None:
                    return None
                coeffs.append(c.getConstant())
        return coeffs
    
    def hasOnlyConstantCoefficients(self):
        self.simplifyCoefficients() # TODO
        for coeff in self.getCoefficients():
            if isNumber(coeff):
                continue
            if not coeff.isConstant():
                return False
            
        return True
            
    @property
    def degree(self):
        return len(self._coefficients)-1
    
    @property
    def lowestDegree(self):
        for i in range(self.degree+1):
            if not self.coeffIsZero(i):
                return i
        return 0
    
    def deg0(self):
        return self.degree==0
    
    def replaceNumbersWithRationals(self):
        for i in range(len(self._coefficients)):
            c = self._coefficients[i]
            if isNumber(c):
                self._coefficients[i] = Number.Rational.fromFloat(c)
            else:
                self._coefficients[i].replaceNumbersWithRationals()
    
    # ========================================== Mixed stuff =========================================
    def asRational(self):
        return Rat.RationalFunction(self,1)
    def Inverse(self):
        return self.asRational().Inverse()
    def isMonic(self):
        return objEqualsNumber(self.getLeadingCoefficient(), 1)
    def makeMonic(self):
        if self.isZero():
            return self
        x = self.getLeadingCoefficient()
        poly = Polynomial([x],variable=self.variable)
        return self/poly
    
    def evaluate(self, val):
        return sum(self.getCoefficient(i)*(val**i) for i in range(self.degree+1))
    def completeEvaluate(self, val):
        sum = 0
        for i in range(self.degree+1):
            if isNumber(i):
                c = self.getCoefficient(i)
            else:
                c = self.getCoefficient(i).completeEvalute(val)
            sum += c*(val**i)
        return sum
    
    def hasRationalRoots(self):
        return len(self.getAllRationalRoots())>0
    
    def getRationalRoots(self):
        coeffs = self.getConstantCoefficients()
        lcm = Number.NumberLCMList([(x._q if type(x)==Number.Rational else (1 if isInt(x) else Number.Rational.fromFloat(x)._q)) for x in coeffs])
        newCoeffs = mulObjectToListElements(lcm, coeffs)
        lc = newCoeffs[len(coeffs)-1]
        ac = newCoeffs[0]
        ps = Number.getAllDivisors(int(abs(ac)))
        qs = Number.getAllDivisors(int(abs(lc)))
        roots = []
        for p in ps+map(lambda x: -x, ps):
            for q in qs:
                x = Number.Rational(p,q)
                v = self.evaluate(x)
                if v==0:
                    if x not in roots:
                        roots.append(x)
        return roots
    
    def getRoots(self):
        coeffs = self.getConstantCoefficients()
        if coeffs==None:
            raise NotImplementedError()
        if self.degree>2:
            raise NotImplementedError()
        
        if self.degree==1:
            return [-coeffs[0]/coeffs[1]]
        elif self.degree==2:
            a = coeffs[2]
            p = coeffs[1]/a
            q = coeffs[0]/a
            rad = (p/2)**2-q
            disc = Number.sqrt(rad)
            zero0 = -p/2+disc
            zero1 = -p/2-disc
            return [zero0,zero1]
        
        
     # ========================================== Field Tower/Extension =========================================
            
    @property
    def fieldTower(self):
        return self.getFieldTower()
    def getFieldTower(self):
        return self.variable.fieldExtension.fieldTower
    def getFieldExtension(self):
        return self.getFieldTower().getLastExtension()

    """def increaseToFieldTower(self, targetFieldTower):
        if targetFieldTower.isExtendedTowerOf(self.getFieldTower()):
            return Polynomial([self],targetFieldTower)
        else:
            raise Exception("")"""
    def reduceToFieldTower(self, targetFieldTower):
        if targetFieldTower==self.fieldTower:
            return self
        if not self.fieldTower.isExtendedTowerOf(targetFieldTower):
            raise Exception()
        if self.deg0():
            return self.getCoefficient(0)
        return None
    
    def reduceToLowestPossibleFieldTower(self):
        if self.fieldTower.towerHeight==1:
            return self
        r = self.reduceToFieldTower(self.fieldTower.prevTower())
        if r==None:
            return self
        elif isNumber(r):
            return self
        else:
            return r.reduceToLowestPossibleFieldTower()
    def isLowestFieldTower(self):
        r = self.reduceToFieldTowr(self.fieldTower.prevTower())
        return r==None

    # ========================================== SquareFree stuff =========================================
    def getLogGCD(self):
        """
        returns the 'logarithmic' gcd: gcd(u,du/dT) # where T is the variable of self
        """
        if id(self.__logGCD)==id(None):
            self.__logGCD = PolyGCD(self, self.differentiateWRTtoPolyVar())
        return self.__logGCD
    
    def isSquareFree(self): # polynomial p is square-free iff gcd(p,(d/dT)p)=1
        gcd = self.getLogGCD()#PolyGCD(self,self.differentiate())
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
        c = self.getLogGCD() # = a_2*(a_3)^2*(a_4)^3...
        w = self/c # = a_1*a_2*a_3...
        y = PolyGCD(c, w) # = a_2*a_3...
        a1 = w/y
        if isNumber(a1):
            assert a1==1
        else:
            if a1.getConstant()==1:
                a1 = 1
        rest = c.factorSquareFree()
        an = [(a1,1)]+rest
        for i in range(1,len(an)):
            an[i] = (an[i][0],an[i][1]+1)
        return an
    
    # ========================================== Differentiate stuff =========================================
    def differentiate(self):
        if id(self.__derivative)!=id(None):
            return self.__derivative
        
        fieldTower = self.getFieldTower()
        dPoly = Polynomial(variable=self.variable)
        fieldExtension = self.getFieldExtension()
        if fieldTower.towerHeight>1:
            u = fieldExtension.argFunction
            a = u.reduceToFieldTower(u.fieldTower.prevTower())
            if a!=None:
                u = a

            uP = u.differentiate()
            log_diff_u = uP/u
        for i in range(self.degree,-1,-1):
            p = self.getCoefficient(i)
            if p==0:
                continue
            if i==0 and fieldExtension!=None:
                if isNumber(p):
                    continue
                else:
                    dPoly += p.differentiate()
                    continue
                
            if fieldTower.towerHeight==1:
                if i>0:
                    dPoly.setCoefficient(i-1,p*i,callUpdates=False)
            elif fieldExtension.extensionType==FE.TRANS_EXP: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'*T^i since T'=u'T
                dPoly += Monomial(i, p.differentiate(), variable=self.variable)
                dPoly += i*Monomial(i,p*uP,variable=self.variable)
            elif fieldExtension.extensionType==FE.TRANS_LOG: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'/u*T^(i-1)
                dPoly += Monomial(i,p.differentiate(),variable=self.variable)
                dPoly += i*Monomial(i-1, p* log_diff_u,variable=self.variable)
                

        dPoly.updateCoefficientsAll()
        if type(dPoly)==Rat.RationalFunction:
            red = dPoly.reduceToLowestPossibleFieldTower()
            if red.fieldTower.towerHeight<self.fieldTower.towerHeight:
                self.__derivative = Polynomial([red],variable=self.variable)
            else:
                if dPoly.denominator.getConstant()==1:
                    self.__derivative = dPoly.numerator
                else:
                    self.__derivative = dPoly.numerator/dPoly.denominator
        else:
            self.__derivative = dPoly
        return self.__derivative
    
    def differentiateWRTtoPolyVar(self):
        """
        differentiates function with respect to the field extension variable of the polynomial:
        p = 2*T**2+T+1
        -> differentiateWRTtoPolyVar(p) = (d/dT)p = 4T+1
        """
        dPoly = Polynomial(variable=self.variable)
        for i in range(self.degree,0,-1):
            p = self.getCoefficient(i)
            dPoly.setCoefficient(i-1,p*i,callUpdates=False)
        dPoly.updateCoefficientsAll()
        return dPoly
    
    def logDifferentiate(self):
        if id(self.__logDerivative)!=id(None):
            return self.__logDerivative
        self.__logDerivative = self.differentiate()/self
        return self.__logDerivative
    
    # ========================================== Arithmetic stuff =========================================
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if id(other)==id(0):
            return self
        if type(other)==Rat.RationalFunction:
            return other.__add__(self)
        if isNumber(other):
            deg0Part = self.getCoefficient(0)
            tPoly = Polynomial(variable=self.variable)
            for i in range(self.degree,0,-1):
                tPoly.setCoefficient(i, self.getCoefficient(i), callUpdates=False)
            tPoly.setCoefficient(0,deg0Part+other,callUpdates=False)
            return tPoly
        if type(other)==ExtPol.ExtendedPolynomial:
            return other.__add__(self)
        
        if self.variable != other.variable:
            if self.fieldTower.isExtendedTowerOf(other.fieldTower):
                return self.__add__(Polynomial([other],variable=self.variable))
            elif other.fieldTower.isExtendedTowerOf(self.fieldTower):
                return other.__add__(self)
            else:
                raise Exception("Polynomials have too different field towers")
        
        tDeg = max(self.degree,other.degree)
        tPoly = Polynomial(variable=self.variable)
        for i in range(tDeg,-1,-1):
            a = self.getCoefficient(i)
            b = other.getCoefficient(i)
            tPoly.setCoefficient(i,a+b,callUpdates=False)
        tPoly.updateCoefficientsAll()
        return tPoly
    
    def __rsub__(self, other): # other-self
        return (-self).__add__(other)
    def __sub__(self, other): # self -other
        return self.__add__((-1)*other)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        if isNumber(other):
            if other==1:
                return self
            elif other==0:
                return Number.ZERO
            tPoly = Polynomial(variable=self.variable)
            for i in range(self.degree,-1,-1):
                tPoly.setCoefficient(i, self.getCoefficient(i)*other,callUpdates=False)
            return tPoly

        if type(other) == Rat.RationalFunction:
            if other.numerator.isZero():
                return Number.ZERO
            return other.__mul__(self)
        if self.variable != other.variable:
            if self.fieldTower.isExtendedTowerOf(other.fieldTower):
                return self.__mul__(Polynomial([other],variable=self.variable))
            elif other.fieldTower.isExtendedTowerOf(self.fieldTower):
                return other.__mul__(self)
            else:
                raise Exception("Polynomials have too different field towers")

        
        tDeg = self.degree+other.degree
        tPoly = Polynomial(variable=self.variable)
        # calcs new coefficients using Cauchy-Product
        for i in range(tDeg,-1,-1):
            newCoeff = 0
            for j in range(0,i+1):
                a = self.getCoefficient(j)
                b = other.getCoefficient(i-j)
                newCoeff += a*b
            tPoly.setCoefficient(i,newCoeff,callUpdates=False)
            
        #tPoly.updateCoefficientsAll()            
        return tPoly
    
    def __rtruediv__(self, other): # other/self
        return Rat.RationalFunction(other,self)
    def __truediv__(self, other): # self/other
        if self.isZero():
            return self
        if isNumber(other):
            return self.__mul__(1/other)
        if type(other) is Rat.RationalFunction:
            if other.fieldTower.towerHeight<self.fieldTower.towerHeight:
                return self.__mul__(Polynomial([other.Inverse()],variable=self.variable))
            else:
                return self.__mul__(other.Inverse())
        (quot,rem)= PolyDiv(self, other)
        if rem==0:
            return quot
        
        x = Rat.RationalFunction(self,other)
        if x.denominator==1:
            return x.numerator
        return x
    
    def isMultipleOf(self, other):
        (s,r) = PolyDiv(self, other)
        return r==0 or r.isZero()
    
    def isConstantMultipleOf(self, other):
        c = (self/other).getConstant()
        return False if c==None else c
        
    def __pow__(self, other):
        if not isNumber(other) or int(other)!=other or other<0:
            raise Exception("power has to be an integer >= 0")
        if other==0:
            return 1
        return (self**(other-1))*self
    
    def __neg__(self):
        return (-1)*self
    
    # ========================================== Object comparison =========================================
    def __eq__(self, other):
        if other==None:
            return False
        return polyEqual(self, other)
    def __ne__(self, other):
        return not self.__eq__(other)

     # ========================================== String output =========================================
    def __str__(self):
        return self.strCustomVar(self.variable)
    
    def __repr__(self):
        return self.__str__()

    def strCustomVar(self, variable):
        out = ""#FE.VARIABLES[self.field]
        for i in range(len(self._coefficients)-1,-1,-1):
            if self.coeffIsZero(i):
                continue
            coeff = ""
            coeff_v = self.getCoefficient(i)
            if coeff_v != 1 or i==0:
                if isNumber(coeff_v):# and coeff_v>0:
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant():# and coeff_v.getConstant()>0:
                    if coeff_v.getConstant()!=1 or i==0:
                        c = coeff_v.getConstant()
                        if type(c)==Number.SqrRootPolynomial and c.a!=0:
                            coeff = "("+str(c)+")"
                        else:
                            coeff = str(coeff_v.getConstant())
                #elif isNumber(coeff_v):
                #    coeff = "({})".format(str(coeff_v))
                else:
                    c = str(coeff_v)
                    if (coeff_v.isConstant() and c.startswith("(") and c.endswith(")")) or i==0:
                        coeff = "{}".format(c)
                    else:
                        coeff = "({})".format(c)
                
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
    
    def printFull(self,reverse=False):
        """
        writes out field extension variables
        """
        out = ""
        fE = self.getFieldExtension()
        r = range(len(self._coefficients)) if reverse else range(len(self._coefficients)-1,-1,-1)
        for i in r:
            if self.coeffIsZero(i):
                continue
            coeff = ""
            coeff_v = self.getCoefficient(i)
            if coeff_v != 1 or i==0:
                if isNumber(coeff_v):# and (type(coeff_v)!=complex and coeff_v>0):
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant():# and (type(coeff_v)!=complex and coeff_v.getConstant()>0):
                    if coeff_v.getConstant()!=1 or i==0:
                        coeff = str(coeff_v.getConstant())
                else:
                    c = coeff_v.printFull()
                    if (coeff_v.isConstant() and c.startswith("(") and c.endswith(")")) or i==0:
                        coeff = "{}".format(c)
                    else:
                        coeff = "({})".format(c)

            power = ""
            if i>1:
                power = "**{}".format(i)
            elif i==1:
                power = ""
            var = ""
            if i>0:
                if fE.argFunction==None:
                    var += fE.variable.stringRepr
                else:
                    var += fE.strFunc()
            
            out += coeff +"{}{}+".format(var,power)
        if len(out)==0:
            return "0"
        return out.strip("+")
    
    

def polyEqual(A,B):
    if isNumber(A) and isNumber(B):
        return A==B
    if isNumber(B):
        return A.getConstant()==B
    if isNumber(A):
        return B.getConstant()==A
    
    return (A+(-1)*B).isZero()

def PolyDiv(polA, polB):
    """
    Polynomial division, returns s,r where polA = s*polB + r, with deg(r)<deg(polB)
    """
    if isNumber(polB) or isNumber(polA):
        raise Exception()
    if polA.fieldTower!=polB.fieldTower:
        if polA.fieldTower.isExtendedTowerOf(polB.fieldTower):
            return PolyDiv(polA,Polynomial([polB],variable=polA.variable))
        elif polB.fieldTower.isExtendedTowerOf(polA.fieldTower):
            return PolyDiv(Polynomial([polA],variable=polB.variable), polB)
        else:
            raise Exception()
      
    (s,r) = _PolyDiv(polA._coefficients, polB._coefficients)
    s = 0 if s==0 else Polynomial(s,variable=polA.variable)
    r = 0 if r==0 else  Polynomial(r,variable=polA.variable)
    return (s,r)
    """
     
    '''if polB==1:
        return (polA,0)
    if (polA.fieldTower!=polB.fieldTower):
        if polA.fieldTower.isExtendedTowerOf(polB.fieldTower):
            return PolyDiv(polA,Polynomial([polB],fieldTower=polA.fieldTower))
        elif polB.fieldTower.isExtendedTowerOf(polA.fieldTower):
            if polB.isLowestFieldTower():
                return (0,polB)
            else:
                return PolyDiv(polA,polB.reduceToLowestPossibleFieldTower()) 
        else:
            raise Exception()
        #raise Exception("Polynomials have to have coefficients in the same field in order to apply PolyDiv")
    if polA.deg0() and polB.deg0():
        if isNumber(polB.getCoefficient(0)):
            return (Polynomial([polA.getCoefficient(0)/polB.getCoefficient(0)],fieldTower=polA.fieldTower),0)#Rat.RationalFunction(polA.getCoefficient(0),polB.getCoefficient(0),polA.field-1)],field=polA.field), 0)
        else:
            if type(polB.getCoefficient(0))==Polynomial or isPoly(polB.getCoefficient(0)):
                return (Polynomial([polA.getCoefficient(0)*polB.getCoefficient(0).asRational().Inverse()],fieldTower=polA.getFieldTower()),0) 
            else:
                return (Polynomial([polA.getCoefficient(0)*polB.getCoefficient(0).Inverse()],fieldTower=polA.getFieldTower()),0) 
    '''
    if (polA.degree<polB.degree):
        return (0,polA)
    power = polA.degree-polB.degree
    if isNumber(polB.getLeadingCoefficient()):
        coeff = polA.getLeadingCoefficient()/polB.getLeadingCoefficient()
    else:
        if isPoly(polB.getLeadingCoefficient()):
            coeff = polA.getLeadingCoefficient()*polB.getLeadingCoefficient().asRational().Inverse()
        else:
            coeff = polA.getLeadingCoefficient()*polB.getLeadingCoefficient().Inverse()
            
    monomial = Monomial(power,coeff,fieldTower=polB.getFieldTower())
    sub = (-1)*(monomial*polB)
    newA = polA+sub # 
    
    if newA.degree==polA.degree: # calculation errors
        newA.setCoefficient(0,newA.degree)
        
    if newA.isZero():
        return (monomial,0)
    
    else:
        (quot,remainder) = PolyDiv(newA,polB)
        return (monomial+quot,remainder)"""
def _PolyDiv(coeffsA, coeffsB):
    degA = len(coeffsA)-1
    degB = len(coeffsB)-1
    if degA<degB:
        return ([0],coeffsA)
    power = degA-degB
    coeff = coeffsA[-1]/coeffsB[-1]
    if isNumber(coeffsB[-1]):
        coeff = coeffsA[-1]/coeffsB[-1]
    else:
        if isPoly(coeffsB[-1]):
            coeff = coeffsA[-1]*coeffsB[-1].asRational().Inverse()
        else:
            coeff = coeffsA[-1]*coeffsB[-1].Inverse()
    monomial = [0]*power+[coeff] 
    sub = mulListWithList(monomial,mulObjectToListElements(-1,coeffsB))
    newA = addListToList(sub, coeffsA)
    newA = newA[0:len(newA)-1]
    if listIsZero(newA):
        return (monomial,0)
    else:
        (quot, remainder)  = _PolyDiv(newA, coeffsB)
        return (addListToList(quot,monomial),remainder)
@call_counter     
def PolyGCD(polA,polB):
    """
    calculates gcd(polA,polB) (monic)
    """
    g = _PolyGCD(polA._coefficients, polB._coefficients)
    if len(g)==1:
        return 1
    else:
        p = Polynomial(g,variable=polA.variable)
        return p.makeMonic()
    """if polB==0:
        return polA
    if polA==0:
        return polB
    (A,B) = (polA,polB) if polA.degree>=polB.degree else (polB,polA)
    if B.isZero():
        return A
    (s,r) = PolyDiv(A,B) # A=s*B+r
    if (isNumber(r) and numberIsZero(r)) or (not isNumber(r) and r.isZero()):
        if B.deg0():
            return 1
        else:
            return B.makeMonic()

    gcd = PolyGCD(B,r)
    if gcd!=1:
        if not gcd.isMonic():
            return gcd.makeMonic()
        else:
            return gcd
    else:
        return 1"""
def _PolyGCD(coeffsA, coeffsB):
    if listIsZero(coeffsA):
        return coeffsB
    if listIsZero(coeffsB):
        return coeffsA
    degA = len(coeffsA)-1
    degB = len(coeffsB)-1
    (A,B) = (coeffsA,coeffsB) if degA>=degB else (coeffsB,coeffsA)
    (s,r) = _PolyDiv(A, B)
    r = 0 if r==0 else listStripZeros(r)
    if id(r)!=id(0):
        lc = r[-1]
        r = mulObjectToListElements(1/lc, r)
    if r==0 or listIsZero(r):
        return B
    gcd = _PolyGCD(B,r)
    return gcd

def __extendedEuclid(p,q):
    """
    finds x,y with
    xp+yq=gcd(p,q), deg(x)<deg(q)-deg(gcd(p,q)),deg(y)<deg(p)-deg(gcd(p,q))
    """
    (P,Q) = (p,q)# if p.degree>=q.degree else (q,p)
    (s,r) = PolyDiv(P,Q)
    rm = 0 if r==0 else r.makeMonic()
    if polyEqual(PolyGCD(p,q),rm):
        if r.getFieldTower().towerHeight==0:#####update !!!!!!!!!!!! TODO
            r_inv_poly = Polynomial([1/r.getLeadingCoefficient()])
            return (r_inv_poly,r_inv_poly*(-1)*s)
        else:
            r_inv_poly = Polynomial([r.getLeadingCoefficient().Inverse()],variable=r.variable)#Rat.RationalFunction(1,r,field=r.field)
            return (r_inv_poly,r_inv_poly*(-1)*s)
    (x,y) = extendedEuclid(Q, r)
    return (y,x+(-1)*s*y)
    
def extendedEuclid(p,q):
    """
    finds x,y with
    xp+yq=gcd(p,q), deg(x)<deg(q)-deg(gcd(p,q)),deg(y)<deg(p)-deg(gcd(p,q))
    """
    (P,Q,a) = (p,q,0) if p.degree>=q.degree else (q,p,1)
    if Q.deg0():
        return (1,((-1)*P+1)/Q)
    (x,y) = __extendedEuclid(P, Q)
    return (x,y) if a==0 else (y,x)

def extendedEuclidGenF(p,q,f):
    """
    finds x,y with
    xp+yq=f
    gcd(p,q)|f
    """
    gcd = PolyGCD(p,q)
    if gcd!=1 and PolyDiv(f,gcd)[1]!=0:
        raise Exception("gcd(p,q) doesn't divide f")
    (s,t) = extendedEuclid(p, q)
    if f.deg0():
        return (s*f,t*f)
    sig1 = s*(f/gcd)
    tau1 = t*(f/gcd)
    (a,r) = PolyDiv(sig1, q/gcd)
    sig = r
    tau = tau1+a*(p/gcd)
    return (sig,tau)


def BasicPartialFraction(f,p,q):
    """
    finds x,y with f/(pq)=x/p+y/q
    gcd(p,q)=1
    returns (x,y)
    """
    if f.getFieldTower()!=p.getFieldTower() or f.getFieldTower()!=q.getFieldTower():
        raise Exception()
    if PolyGCD(p,q)!=1:
        raise Exception("p, q shouldn't have a common divisor.")
    (sig,tau) = extendedEuclidGenF(p, q, f)
    return (tau,sig)


def PartialFractionWithPowerFactorization(numerator,factorization):
    """
    factorization = squarefree factorization of the denominator = a1* a2^2 * a3^3 * ...
    numerator/factorization = r11/a1 + r21/a2 + r22/a2^2 + r31/a3 + ..., deg(rij)<deg(ai) 
    returns [(r11,a1,1),(r21,a2,1),(r22,a2,2),...]
    """
    if len(factorization)==1:
        return PartialFractionPower(numerator, factorization[0][0], factorization[0][1])
    
    otherFactors = 1
    for i in range(1,len(factorization)):
        otherFactors *= factorization[i][0]**factorization[i][1]
        
    if factorization[0][0]!=1:
        (r1,r2) = BasicPartialFraction(numerator, factorization[0][0], otherFactors)
        pfrac = PartialFractionPower(r1, factorization[0][0], factorization[0][1])
        return pfrac + PartialFractionWithPowerFactorization(r2, factorization[1:len(factorization)])
    else:
        return PartialFractionWithPowerFactorization(numerator, factorization[1:len(factorization)])
    
    
    
def PartialFractionPower(p,q,n):
    """
    deg(p)<deg(q)
    returns [(r1,q,1),(r2,q,2),(r3,q,3),...] with
    p/(q^n) = r1/q+r2/q^2+...+rn/q^n
    deg(ri)<deg(q)
    """
    if p==0:
        return []
    if n==1:
        if p.degree>=q.degree:
            raise Exception("deg(p) should be < deg(q)")
        return [(p,q,1)]
    (L,J) = PolyDiv(p,q)
    return PartialFractionPower(L, q, n-1)+[(J,q,n)]

def Monomial(degree, coeff,variable=None):
    """
    T = variable of last extension of fieldTower
    returns coeff * (T)^degree
    """
    var = variable
    if var==None:
        var = FE.BASEVARIABLE
    return Polynomial([0]*degree+[coeff], variable=var)
    


def printFactorization(fact):
    """
    returns string representation of factorization
    fac = [(a1,1),(a2,2),...]
    returns "{a1}*{a2}**2*{a3}**3..."
    """
    out = ""
    for i in range(len(fact)):
        if fact[i][0]!=1:
            out += "({})".format(fact[i][0])
            if fact[i][1]!=1:
                out += "**{}".format(fact[i][1])
            out += "*"
    return out.strip("*")
            
def printPartialFraction(pfrac):
    """
    returns string representation of partial fraction decomposition
    """
    out = ""
    for frac in pfrac:
        if (isNumber(frac[0]) and numberIsZero(frac[0])) or frac[0].isZero():
            continue
        power = ""
        if frac[2]>1:
            power = "**{}".format(frac[2])
        out += "[{}]/[{}]{}+".format(frac[0],frac[1],power)
    return out.strip("+")

def polytest(expected, got):
    return "polytest: should be {}, got {}".format(expected, got)

if __name__ == '__main__': #tests
    from Parse import parseField0PolyFromStr,parseField0RatFromStr,parseExpressionFromStr
    
    [T1,T2] = FE.Variables(['T_1','T_2'])
    
    
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Polynomial([0,Number.ONE]),T1) # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Polynomial([0,Number.ONE],variable=T1),T2) # field extension with log(log(x))
    FE.fieldTower = FE.fieldTower.addFieldExtension(fieldExtension2)
    
    p1 = Polynomial(coefficients=[0,0,1,2,3,4])
    p1.replaceNumbersWithRationals()
    print(p1)
    
    p2 = Polynomial(coefficients=[Number.ONE, p1],variable=T1)
    print(p2)
    
    print(p1.fieldTower)
    print(p2.fieldTower)
    print(p1+p2)
    
    print("Field Tower: {}".format(FE.fieldTower))
    
    polA = Polynomial([1,1])
    print("polytest: should be x+1, got {}".format(polA))
    polB = Polynomial([1,1],variable=T1)
    print("polytest: should be T_1+1, got {}".format(polB))
    #print("polytest: should be rationalfunction, got {}".format(type(polB.getCoefficient(0))))
    #print("polytest: should be 0, got {}".format(polB.getCoefficient(0).getFieldTower().towerHeight))
    ratA = Rat.RationalFunction(1,polA)
    print("polytest: should be 1/[x+1], got {}".format(ratA))
    
    polD = Polynomial([ratA,3],variable=T1)
    print("polytest: should be 3T_1+(1/[x+1]), got {}".format(polD))
    #print("polytest: should be 2x 0, got {}, {}".format(polD.getCoefficient(0).getFieldTower().towerHeight,polD.getCoefficient(1).getFieldTower().towerHeight))
    
    polC = Polynomial([ratA,1], variable=T2)
    print("polytest: should be T2+(1/[x+1]), got {}".format(polC))
    #print("polytest: should be 2x RationalFunction, got {},{}".format(type(polC.getCoefficient(0)),type(polC.getCoefficient(1))))
    #print("polytest: should be 2, got {}".format(polC.getFieldTower().towerHeight))
    #print("polytest: should be 2x 1, got {},{}".format(polC.getCoefficient(0).getFieldTower().towerHeight,polC.getCoefficient(1).getFieldTower().towerHeight))
    
    polE = Polynomial([0],variable=T2)
    print(polytest("0",polE))
    print(polytest("3",polE.getFieldTower().towerHeight))
    #print(polytest("Rationalfunction",type(polE.getCoefficient(0))))
    
    
    print("======================== Arithmetic stuff tests ==============================")
    polE = polB+polD
    print("polytest: should be 4T1+([x+2]/[x+1]), got {}".format(polE))
    polF = polE+polD
    print("polytest: should be 7T1+([x+3]/[x+1]), got {}".format(polF))
    
    polA = Polynomial([1,1])*Polynomial([-1,1]) # (x+1)(x-1)
    polB = Polynomial([1,1])*Polynomial([1,0,1])# (x+1)(x^2+1)
    (s,r) = PolyDiv(polB, polA)
    print(polytest("(x+1,2x+2)",(s,r)))
    gcd = PolyGCD(polA, polB)
    print(polytest("x+1",gcd))
    
    polA = parseExpressionFromStr("(T_2-T_1) * (T_2+T_1)",FE.fieldTower)
    print(polytest("T_2**2+((-1)*T_1**2)",polA))
    polB = parseExpressionFromStr("(T_2+x**2)*(T_2-T_1)*(T_2-x)",FE.fieldTower)
    print(polytest("T_2 +(-1)T_1", PolyGCD(polA, polB)))
    
    print("======================== euclid stuff tests ==============================")
    polA = parseField0PolyFromStr("x**2+2*x+2")
    polB = parseField0PolyFromStr("2*x+1")
    (s,r) = PolyDiv(polA, polB)
    print("polytest: should be 0.5x+0.75, 1.25, got {}, {}".format(s,r))
    gcd = PolyGCD(polA, polB)
    print("polytest: should be 1, got {}".format(gcd))
    (x,y) = extendedEuclid(polA, polB)
    print("polytest: should be 0.8, (-0.4)x+(-0.6), got {}, {}".format(x,y))
    
    print("============================= Old poly/rat class tests ===================================")
    print("---------- General poly tests") #################################################################################
    polA = Polynomial(coefficients=[1,0,2,4])
    print(polytest("4x**3+2x**2+1",str(polA)))
    polB = Polynomial(coefficients=[0,1,2])
    print(polytest("2x**2+x",str(polB)))
    polC = polA+polB
    #print("{} + {} = {}".format(polA,polB,polC))
    print(polytest("4x**3+4x**2+x+1",str(polC)))
    polD = polA*polB
    # print("[{}] * [{}] = {}".format(polA,polB,polD))
    print(polytest("8x**5+8x**4+2x**3+2x**2+x",str(polD)))
    
    polF_T = Polynomial(coefficients=[polA,polB],variable=T1) # polB*(T_1)**1+polA = 
    print(polytest("(2x**2+x)T1+(4x**3+2x**2+1)",str(polF_T)))#polF_T)
    polG_T = polF_T*polF_T
    print(polytest("(4x**4+4x**3+x**2)T1**2+(16x**5+16x**4+4x**3+4x**2+2x)T1+(16x**6+16x**5+4x**4+8x**3+4x**2+1)",polG_T))#print("[{}] * [{}] = {}".format(polF_T,polF_T,polG_T))
    
    polH = Polynomial(coefficients=[0,Polynomial([1])], variable=T1)
    polG = Polynomial(coefficients=[polH,polH],variable=T2)
    print(polytest("T1",polH))
    print(polytest("(T1)T2+(T1)",polG))
    print(polytest("(log(x))log(log(x))+(log(x))",polG.printFull()))
    
    ratA = parseField0RatFromStr("1/x")
    ratB = parseField0RatFromStr("x**2/x+1")
    polA = Polynomial(coefficients=[ratA,ratB],variable=T1) # x^2/(x+1)*T1+1/x
    print(polytest("T1+([x+1]/[x**3])",polA.makeMonic()))
    
    #polA = Polynomial([0,1]) # x in C(x)
    #polB = Polynomial([0,1],fieldTower=FE.fieldTower.getStrippedTower(1))# T in C(x,T)
    #print(polA,polB)
    #polC = polA.increaseField(1)
    #print(polC,polC.field)
    #print(polA*polB)
    print("---------- Diff poly tests") #################################################################################
    pol = parseExpressionFromStr("(x/(x+0.5)) * T_1 - x/(x+0.5)", FE.fieldTower.getStrippedTower(3))
    print(pol)
    print(pol.differentiate())
    print(polytest("([x**2+2x]/[x**2+2x+1])T1+([x**3+(-1.0)x+(-1.0)]/[x**3+x**2])",polA.differentiate()))#"[{}]'={}".format(polA,polA.differentiate()))
    print(polytest("((1.0/[x]))T2+(2.0/[x])", polG.differentiate()))#[{}]'={}".format(polG,polG.differentiate()))
    
    print("---------- Poly div test:") #################################################################################
    polA = Polynomial(coefficients=[11,4,0,3]) # 3x^3+4x^2+11
    polB = Polynomial(coefficients=[2,-3,1])#x^2-3x+2
    (quot,rem) = PolyDiv(polA, polB)
    print(polytest("3x+9, 25x+(-7)",(quot,rem)))#print("[{}] / [{}] = {} R {}".format(polA,polB,quot,rem))
    
    polA = Polynomial(coefficients=[18,9,-11,12]) # 12x^3-11x^2+9x+18
    polB = Polynomial(coefficients=[3,4])  #4x+3
    (quot,rem) = PolyDiv(polA, polB)
    print(polytest("3x**2+(-5)x+6, 0", (quot,rem)))#"[{}] / [{}] = {} R {}".format(polA,polB,quot,rem)) # (12x^3-11x^2+9x+18)/(4x+3)
    #polA = parseField0PolyFromStr("x**2+")
    polA = Polynomial(coefficients=[parseField0PolyFromStr("(-1)*x**2"),0,1],variable=T1)
    polB = Polynomial(coefficients=[parseField0PolyFromStr("x"),1],variable=T2)
    (quot,rem) = PolyDiv(polA, polB)
    print(polytest("T1+(-1)x, 0",(quot,rem)))
    (quot,rem) = PolyDiv(polG_T, polF_T)
    print(polytest("(2x**2+x)T1+(4x**3+2x**2+1), 0", (quot,rem)))#"[{}] / [{}] = {} R {}".format(polG_T,polF_T,quot,rem))
    
    print("---------- GCD Tests") #################################################################################
    polA = Polynomial(coefficients=[-4,0,1])#x^2-4
    polB = Polynomial(coefficients=[2,1])#x+2
    print(polytest("x+2",PolyGCD(polA, polB)))
    
    polA = Polynomial(coefficients=[-1,0,1])#x^2-1
    polB = Polynomial(coefficients=[1,2,1])#1+2x+x^2
    print(polytest("x+1",PolyGCD(polA, polB)))
    
    polA = Polynomial(coefficients=[1,0,1])
    polB = Polynomial(coefficients=[-1,1])
    print(polytest("1", PolyGCD(polA, polB)))
    
    polA = Polynomial(coefficients=[1,1,2,0,1,0,0,1])#1+x+2x^2+x^4+x^7
    polB = Polynomial(coefficients=[1,1,1,1,1,1])#1+x+x^2+x^3+x^4+x^5
    print(polytest("1",PolyGCD(polA, polB)))
    
    ratA = parseField0RatFromStr("1/1")
    ratB = parseField0RatFromStr("(-1)*x**2/1")
    polA = Polynomial([ratB,0,ratA],variable=T1) # = -x^2+T^2=(T-x)(T+x)
    ratC = parseField0RatFromStr("x**2+1/x")
    polB = Polynomial([ratA,ratC,ratA],variable=T1)# = T^2+(x^2+1)/x*T+1=(T+x)(T+1/x)
    print(polytest("T1+(x)",PolyGCD(polA, polB)))
    
    (x,y) = extendedEuclid(polA,polB)
    print(polytest("([x]/[(-1.0)x**2+(-1.0)]), ([(-1.0)x]/[(-1.0)x**2+(-1.0)])",(x,y)))
    #print("extended euclid: {}, {}, {}".format(x,y,PolyGCD(polA, polB)))
    f = PolyGCD(polA, polB)*Polynomial([Polynomial([1,1]),1],variable=T1) #*(1+x+T+T^2)
    (x,y) = extendedEuclidGenF(polA, polB, f)
    print(polytest("([(-1.0)x**2+(-1.0)x+1.0]/[x**2+1.0]), ([2.0x**2+x]/[x**2+1.0])",(x,y)))

    print(polytest("True",polA.isSquareFree()))
    
    polC = Polynomial([Polynomial([0,1]).asRational(),ratA],variable=T1)
    polD = polC*polA # (T+x)*(T-x)*(T+x)
    print(polytest("False",polD.isSquareFree()))

    print(polytest("(T1+((-1.0)x))*(T1+(x))**2",printFactorization(polD.factorSquareFree())))

    print(polytest("T1**2+((-1)*x**2)",polD/polC))
    
    print("---------- ExtendedEuclid Test") #################################################################################
    polA = parseField0PolyFromStr("x**3+x**2+(-1)*x+2")
    polB = parseField0PolyFromStr("x**2+1")
    (x,y) = extendedEuclid(polA,polB)
    print(polytest("1",PolyGCD(polA, polB)))
    
    print(polytest("0.4x+0.2, (-0.4)x**2+(-0.6)x+0.6",(x,y)))
    f = parseField0PolyFromStr("x**2+(-1)")
    (x,y) = extendedEuclidGenF(polA, polB, f)
    print(polytest("(-0.8)x+(-0.4), 0.8x**2+1.2x+(-0.2)",(x,y)))
    
    print("---------- Factor square free test") #################################################################################
    polA = parseField0PolyFromStr("x**3+4*x**2+5*x+2")
    fact = polA.factorSquareFree()
    print(polytest("(x+2.0)*(x+1.0)**2",printFactorization(fact)))
    #print("{} = {}".format(str(polA),printFactorization(fact)))
    
    ratA = parseField0RatFromStr("x+1/x+2")
    polA = Polynomial([ratA,Number.ONE],variable=T1)**2
    print(polytest("False",polA.isSquareFree()))
    fact = polA.factorSquareFree()
    print(polytest("(T1+([x+1.0]/[x+2.0]))**2",printFactorization(fact)))
    #print("{} = {}".format(polA,printFactorization(fact)))
    assert polA.isSquareFree()==False
    
    ratB = parseExpressionFromStr("(-1)*x/1",FE.fieldTower)
    polB = Polynomial([ratB,Number.ONE],variable=T1)*polA # = (-x+T1)*((x+1)/(x+2)+T1)**2
    print(polytest("False",polB.isSquareFree()))
    fact = polB.factorSquareFree()
    print(polytest("(T1+((-1.0)x))*(T1+([x+1.0]/[x+2.0]))**2",printFactorization(fact)))
    #print("{} = {}".format(polB,printFactorization(fact)))
    assert polB.isSquareFree()==False
    
    # http://www.wolframalpha.com/input/?i=expand+(T-x)(T%2B(x%2B1)%2F(x%2B2))%5E2+at+T%3D0
    ratA = parseExpressionFromStr("(x+2*x**2+x**3)/(4+4*x+x**2)")
    ratB = parseExpressionFromStr("(1+(-2)*x+(-5)*x**2+(-2)*x**3)/(4+4*x+x**2)")
    ratC = parseExpressionFromStr("(2+(-1)*x**2)/(2+x)")
    ratD = parseExpressionFromStr("1/1")
    
    polA = Polynomial([ratA*(-1),ratB,ratC,ratD],variable=T1)
    #print(polA)
    fact = polA.factorSquareFree()
    #print(polyEqual(polA, polB))
    assert polA.isSquareFree()==False
    print(polytest("(T1+((-1.0)x))*(T1+([x+1.0]/[x+2.0]))**2",printFactorization(fact)))
    
    polA = parseField0PolyFromStr("x**4+3*x**3+3*x**2+x")
    print(polytest("(x)*(x+1.0)**3",printFactorization(polA.factorSquareFree())))
    # print("{} = {}".format(str(polA),printFactorization(polA.factorSquareFree())))  
    
    print("---------- partial fraction test") #################################################################################
    
    polA = parseField0PolyFromStr("x**5+x**2+1")
    polB = parseField0PolyFromStr("x**2+1")
    pfrac = PartialFractionPower(polA, polB, 4)
    print(polytest("[x]/[x**2+1]**2+[(-2.0)x+1.0]/[x**2+1]**3+[x]/[x**2+1]**4",printPartialFraction(pfrac)))
    #print("[{}]/[{}]^{} = {}".format(polA,polB,4,printPartialFraction(pfrac)))
    
    polA = parseField0PolyFromStr("x**2+1")
    polB = parseField0PolyFromStr("x**4+3*x**3+3*x**2+x")
    factorization = polB.factorSquareFree()
    pfrac = PartialFractionWithPowerFactorization(polA, factorization)
    print(polytest("[1.0]/[x]+[(-1.0)]/[x+1.0]+[(-2.0)]/[x+1.0]**3",printPartialFraction(pfrac)))
    #print("[{}]/[{}] = {}".format(polA,polB,printPartialFraction(pfrac)))
    # a = parseField0PolyFromStr("x")
    # b = parseField0PolyFromStr("x**3+3*x**2+3*x+1")
    # c = parseField0PolyFromStr("x**2+1")
    # (x,y) = extendedEuclidGenF(a, b, c)
    # print(x,y)  
    
    one = parseField0RatFromStr("1/1")
    x = parseField0RatFromStr("x/1")
    polA = Polynomial([(-1)*x,one],variable=T1) # T-x
    polB = Polynomial([x,one],variable=T1) # T+x
    
    num = Polynomial([one,one],variable=T1)#1+T
    denom = polA*(polB**2) # (T-x)(T+x)^2
    
    factorization = denom.factorSquareFree()
    #print("{} = {}".format(denom,printFactorization(factorization)))
    print(polytest("(T1+((-1.0)x))*(T1+(x))**2",printFactorization(factorization)))
    pfrac = PartialFractionWithPowerFactorization(num, factorization)
    print(polytest("[([0.25x+0.25]/[x**2])]/[T1+((-1.0)x)]+[([(-0.25)x+(-0.25)]/[x**2])]/[T1+(x)]+[([0.5x+(-0.5)]/[x])]/[T1+(x)]**2",printPartialFraction(pfrac)))
    # print("[{}]/[{}] = {}".format(num,denom,printPartialFraction(pfrac)))
        