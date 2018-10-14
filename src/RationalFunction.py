'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
import Polynomial as Pol
from Utils import *
import Number

class RationalFunction(object):
    def __init__(self, numerator, denominator):
        
        self.__derivative = None
        self.__logDerivative = None
         
         
        self.numerator = numerator
        self.denominator = denominator
        
        while type(self.numerator) is RationalFunction or type(self.denominator) is RationalFunction:
            if type(self.numerator)==RationalFunction:
                tempNum = self.numerator
                self.numerator = self.numerator.numerator
                self.denominator = self.denominator*tempNum.denominator
            if type(self.denominator)==RationalFunction:
                tempDenom = self.denominator
                self.numerator = self.numerator*self.denominator.denominator
                self.denominator = tempDenom.numerator
        
        if isNumber(self.numerator):
            if isNumber(denominator):
                self.denominator = Pol.Polynomial([self.denominator])
            self.numerator = Pol.Polynomial([self.numerator],variable=self.denominator.variable)#fieldTower=self.denominator.fieldTower)
        elif isNumber(self.denominator):
            self.denominator = Pol.Polynomial([self.denominator],variable=self.numerator.variable)#fieldTower = self.numerator.fieldTower)
            
        numFieldTower = self.numerator.getFieldTower()
        denomFieldTower = self.denominator.getFieldTower()    
        if numFieldTower!=denomFieldTower:
            if numFieldTower.isExtendedTowerOf(denomFieldTower):
                self.denominator = Pol.Polynomial([self.denominator],variable=self.numerator.variable)#fieldTower=numFieldTower)
            elif denomFieldTower.isExtendedTowerOf(numFieldTower):
                self.numerator = Pol.Polynomial([self.numerator],variable=self.denominator.variable)#fieldTower=denomFieldTower)
            else:
                raise Exception()
        
        #self.fieldTower = self.numerator.fieldTower
        self.variable = self.numerator.variable
        
        self.removeCommonFactors()
        
        if self.denominator.isConstant():
            self.numerator = self.numerator/self.denominator.getConstant()
            self.denominator = Pol.Polynomial([1],variable=self.variable)
        elif self.denominator.deg0():
            self.numerator = self.numerator/self.denominator
            self.denominator = Pol.Polynomial([Number.ONE],variable=self.variable)

        self.makeDenominatorMonicInPlace()
        
    @property
    def fieldTower(self):
        return self.variable.fieldExtension.fieldTower
    
    def isConstant(self):
        if isNumber(self.numerator):
            return numberIsZero(self.numerator) or isNumber(self.denominator) or self.denominator.isConstant()
        if self.numerator.isConstant():
            return numberIsZero(self.numerator.getConstant()) or isNumber(self.denominator) or self.denominator.isConstant()
        return False
        
    def getConstant(self):
        if not self.isConstant():
            return None
        if isNumber(self.numerator):
            if numberIsZero(self.numerator):
                return 0
            if isNumber(self.denominator):
                return self.numerator/self.denominator
            else:
                return self.numerator/self.denominator.getConstant()
        else:
            if self.numerator.isZero():
                return 0
            if isNumber(self.denominator):
                return self.numerator.getConstant()/self.denominator
            else:
                return self.numerator.getConstant()/self.denominator.getConstant()
            
    def isZero(self):
        return objEqualsNumber(self.getConstant(),Number.ZERO)

    # ========================================== Field Tower/Extension =========================================
    def getFieldTower(self):
        return self.fieldTower
    def getFieldExtension(self):
        return self.getFieldTower().getLastExtension()
    
    
    def removeCommonFactors(self):
        """
        cancels common factors, so that: numerator/denominator = p/q with gcd(p,q)=1
        """
        gcd = Pol.PolyGCD(self.numerator, self.denominator)
        if not gcd==1:
            self.numerator = Pol.PolyDiv(self.numerator, gcd)[0]
            self.denominator = Pol.PolyDiv(self.denominator,gcd)[0]
            
    def reduceToFieldTower(self, targetFieldTower):
        if isNumber(self.numerator):
            newNum = Pol.Polynomial([self.numerator],variable=targetFieldTower.getLastVariable())
        else:
            newNum = self.numerator.reduceToFieldTower(targetFieldTower)
        if isNumber(self.denominator):
            newDenom = Pol.Polynomial([self.denominator],variable=targetFieldTower.getLastVariable())
        else:
            newDenom = self.denominator.reduceToFieldTower(targetFieldTower)
        if newNum==None or newDenom==None:
            return None
        
        return RationalFunction(newNum,newDenom)#newNum/newDenom#RationalFunction(newNum,newDenom,fieldTower=targetFieldTower)
    
    def reduceToLowestPossibleFieldTower(self):
        if self.fieldTower.towerHeight==1:
            return self
        r = self.reduceToFieldTower(self.fieldTower.prevTower())
        if r==None:
            return self
        else:
            return r.reduceToLowestPossibleFieldTower()
    def isLowestFieldTower(self):
        r = self.reduceToFieldTower(self.fieldTower.prevTower())
        return r==None
    
    # ========================================== Differentiate stuff =========================================
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
            
        num = dp*q+(-1)*p*dq
        denom = q*q
        self.__derivative = RationalFunction(num,denom)
            
        return self.__derivative# (p/q)' = (p'q-pq')/(q^2)
    def logDifferentiate(self):
        if id(self.__logDerivative)!=id(None):
            return self.__logDerivative
        self.__logDerivative = self.differentiate()/self
        return self.__logDerivative
    
    # ========================================== Arithmetic stuff =========================================
    def updateCoefficientsAll(self):
        self.numerator.updateCoefficientsAll()
        self.denominator.updateCoefficientsAll()
        
    def Inverse(self):
        return RationalFunction(self.denominator,self.numerator)
    def makeDenominatorMonic(self):
        if isNumber(self.denominator):
            return
        lcoeff = self.denominator.getLeadingCoefficient()
        lcoeff_poly = Pol.Polynomial(coefficients=[lcoeff],variable=self.variable)
        if isNumber(self.numerator):
            newNumerator = Pol.Polynomial([self.numerator],variable=self.variable)/lcoeff_poly
        else:
            newNumerator = self.numerator/lcoeff_poly
        newDenominator = self.denominator/lcoeff_poly
        return RationalFunction(newNumerator,newDenominator)
    def makeDenominatorMonicInPlace(self):
        if isNumber(self.denominator):
            return
        lcoeff = self.denominator.getLeadingCoefficient()
        lcoeff_poly = Pol.Polynomial(coefficients=[lcoeff],variable=self.variable)
        if isNumber(self.numerator):
            self.numerator = Pol.Polynomial([self.numerator],variable=self.variable)/lcoeff_poly
        else:
            self.numerator = self.numerator/lcoeff_poly
        self.denominator = self.denominator/lcoeff_poly
        
    def asPolynomial(self):
        self.numerator.simplifyCoefficients()
        self.denominator.simplifyCoefficients()
        c = self.denominator.getConstant()
        if c==None:
            denom = self.denominator.reduceToLowestPossibleFieldTower()
            if denom.fieldTower.towerHeight<self.numerator.fieldTower.towerHeight:
                return self.numerator/denom
            (q,r) = Pol.PolyDiv(self.numerator, self.denominator)
            if r==0:
                return q
            return None
        return self.numerator*(1/c)
    
    def simplified(self):
        num = self.numerator.simplified()
        denom = self.denominator.simplified()
        return RationalFunction(num,denom)
        
    def BasicPartialFraction(self, denomFactorization):
        return Pol.PartialFractionWithPowerFactorization(self.numerator, denomFactorization)
    
    def replaceNumbersWithRationals(self):
        if not isNumber(self.numerator):
            self.numerator.replaceNumbersWithRationals()
        if not isNumber(self.denominator):
            self.denominator.replaceNumbersWithRationals()
            
    def completeEvaluate(self, val):
        return self.numerator.completeEvaluate(val)/self.denominator.completeEvaluate(val)
    
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        #if other==0:
        #    return self
        if isNumber(other):
            if other==0:
                return self
            return self.__add__(Pol.Polynomial([other],variable=self.variable))
        if other.isZero():
            return self
        if type(other)==Pol.Polynomial or isPoly(other):
            return self.__add__(other.asRational())
        num = self.numerator*other.denominator+self.denominator*other.numerator
        den = self.denominator*other.denominator
        return RationalFunction(num,den)
    
    def __sub__(self, other):
        return self.__add__((-1)*other)
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        #if other==1:
        #    return self
        #if other==0:
        #    return 0
        if isNumber(other):
            if other==0:
                return 0
            elif other==1:
                return self
            return RationalFunction(self.numerator*other,self.denominator)

        if type(other)==Pol.Polynomial or isPoly(other):
            return self.__mul__(other.asRational())
        num = self.numerator*other.numerator
        den = self.denominator*other.denominator
        return RationalFunction(num,den)
    
    def __neg__(self):
        return (-1)*self
    def __rtruediv__(self, other):
        return self.Inverse()*other
    def __truediv__(self, other):
        if isNumber(other):
            return self.__mul__(1/other)
        if isPoly(other):
            return self.__truediv__(other.asRational())
        newNum = self.numerator*other.denominator
        newDenom = self.denominator*other.numerator
        return RationalFunction(newNum,newDenom)
    
    def isConstantMultipleOf(self, other):
        c = (self/other).getConstant()
        return False if c==None else c
    
    
    def __eq__(self, other):
        if other==None:
            return False
        return (self+(-1)*other).isZero()
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return self.__str__()
    # ========================================== String output =========================================
    def __str__(self):
        out = ""
        if not isNumber(self.numerator) and not self.numerator.isConstant() and self.denominator!=1:
            out = "["+str(self.numerator)+"]"
        elif isNumber(self.numerator):
            out = str(self.numerator)
        elif self.denominator==1:
            out = str(self.numerator)
        else:
            out = str(self.numerator.getLeadingCoefficient())
        if self.numerator==0:
            return "0"
        d = str(self.denominator)
        if self.denominator==1:
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
        out = ""
        if not isNumber(self.numerator) and not self.numerator.isConstant() and self.denominator!=1:
            out = "["+self.numerator.printFull()+"]"
        elif isNumber(self.numerator):
            out = str(self.numerator)
        elif self.denominator==1 or (not isNumber(self.denominator) and objEqualsNumber(self.denominator.getConstant(),1)):
            out = self.numerator.printFull()
        else:# => self.numerator.isConstant() = True
            out = self.numerator.printFull()
        if self.numerator==0:
            return "0"
        if isNumber(self.denominator):
            d = str(self.denominator)
        else:
            d = self.denominator.printFull()
        if self.denominator==1:
            return out
        if len(d)==1:
            return out + "/"+d
        return out+"/["+d+"]"
    

if __name__=='__main__':
    pass
    """FE.fieldTower = FE.FieldTower()
    polA = Pol.Polynomial(coefficients=[1,2,1])
    polB = Pol.Polynomial(coefficients=[2,1])
    ratA = RationalFunction(polA,polB)
    print(ratA)
    polC = Pol.Polynomial(coefficients=[-1,0,1])
    polD = Pol.Polynomial(coefficients=[1,1])
    ratB = RationalFunction(polC,polD)
    print(ratB)
    fieldExtension1 = FE.FieldExtension(FE.TRANS_EXP,Pol.Polynomial([0,1]),"T_1") # field extension with e^x=exp(x)
    fieldExtension2 = FE.FieldExtension(FE.TRANS_EXP,Pol.Polynomial([0,1],field=1),"T_2") # field extension with exp(exp(x))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1,fieldExtension2])
    
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
    x = Pol.Polynomial([0,1])
    polA = Pol.Polynomial([0,x],field=1)
    polB = Pol.Polynomial([0,polA],field=2)
    polC = Pol.Polynomial([polA],field=2)
    print(polB,polC)
    print(polB/polC)
    print(polB.getCoefficients())"""