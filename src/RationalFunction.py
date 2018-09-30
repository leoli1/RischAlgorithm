'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
import Polynomial as Pol
from Utils import *

class RationalFunction(object):
    
    def __init__(self, numerator, denominator, fieldTower=None):
        
        self.numerator = numerator
        self.denominator = denominator
        
        if fieldTower==None:
            self.fieldTower = FE.FieldTower()
        else:
            self.fieldTower = fieldTower
            
        if self.denominator!=1 and self.numerator!=1 and self.numerator!=0:
            if (isNumber(self.numerator) or not self.numerator.isZero()) and not ( isNumber(self.numerator) or isNumber(self.denominator)): # sieht komisch aus, stimmt aber so
                self.removeCommonFactors()
            
        # TODO: clean up if's
        if (isNumber(numerator) or numerator.isConstant()) and (isNumber(denominator) or denominator.isConstant()):
            if isNumber(numerator) and isNumber(denominator):
                self.numerator = numerator/denominator
                self.denominator = 1
            elif isNumber(numerator) and denominator.isConstant():
                self.numerator = Pol.Polynomial([numerator/denominator.getConstant()],fieldTower=self.getFieldTower())#(denominator/numerator).asRational().Inverse().numerator()
                self.denominator = 1
            elif isNumber(denominator) and numerator.isConstant():
                self.numerator = Pol.Polynomial([numerator.getConstant()/denominator],fieldTower=self.getFieldTower())
                self.denominator = 1
            elif numerator.isConstant() and denominator.isConstant():
                self.numerator = Pol.Polynomial([numerator.getConstant()/denominator.getConstant()],fieldTower=self.getFieldTower())
                self.denominator = 1
            
        if isNumber(self.denominator):
            self.numerator = self.numerator/self.denominator
            self.denominator = 1
        elif self.denominator.isConstant():
            self.numerator = self.numerator/self.denominator.getConstant()
            self.denominator = 1
            
        self.makeDenominatorMonic()
        
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
        return objEqualsNumber(self.getConstant(),0)
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
            
        return RationalFunction(dp*q+(-1)*p*dq,q*q,fieldTower=self.getFieldTower())# (p/q)' = (p'q-pq')/(q^2)
    
    # ========================================== Arithmetic stuff =========================================
    def Inverse(self):
        return RationalFunction(self.denominator,self.numerator,fieldTower=self.getFieldTower())
    def makeDenominatorMonic(self):
        if isNumber(self.denominator):
            return
        lcoeff = self.denominator.getLeadingCoefficient()
        lcoeff_poly = Pol.Polynomial(coefficients=[lcoeff],fieldTower=self.getFieldTower())
        if isNumber(self.numerator):
            self.numerator = Pol.Polynomial([self.numerator],fieldTower=self.getFieldTower())/lcoeff_poly
        else:
            self.numerator = self.numerator/lcoeff_poly

            #print(self.denominator,self.denominator.field)
            #print(self.denominator.getCoefficients())
            #print(lcoeff_poly,lcoeff_poly.field)
        self.denominator = self.denominator/lcoeff_poly
        
    def PartialFraction(self, denomFactorization):
        
        self.makeDenominatorMonic()
        return Pol.PartialFractionWithPowerFactorization(self.numerator, denomFactorization)
    
    
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other==0:
            return self
        if isNumber(other):
            return self.__add__(Pol.Polynomial([other],fieldTower=self.getFieldTower()))
            #raise Exception("??")
        if other.isZero():
            return self
        if type(other)==Pol.Polynomial or isPoly(other):
            return self.__add__(other.asRational())
        num = self.numerator*other.denominator+self.denominator*other.numerator
        den = self.denominator*other.denominator
        return RationalFunction(num,den,fieldTower=self.getFieldTower())
    
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        if other==1:
            return self
        if other==0:
            return 0
        if isNumber(other):
            return RationalFunction(self.numerator*other,self.denominator,fieldTower=self.getFieldTower())
            #raise Exception("??")
        if type(other)==Pol.Polynomial or isPoly(other):
            return self.__mul__(other.asRational())
        num = self.numerator*other.numerator
        den = self.denominator*other.denominator
        return RationalFunction(num,den,fieldTower=self.getFieldTower())    
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
        numStr = str(self.numerator) if isNumber(self.numerator) else self.numerator.printFull()
        denomStr = str(self.denominator) if isNumber(self.denominator) else self.denominator.printFull()
        if self.denominator==1:
            return "[{}]".format(numStr)
        return "[{}]/[{}]".format(numStr,denomStr)


'''class RationalFunction(object):


    def __init__(self, numerator, denominator, field = FE.BASE_FIELD,fieldTower=None,output=False):
        self.output = output
        if output:
            print("A")
        self.field = field

        self.numerator = numerator
        self.denominator = denominator
        
        if fieldTower==None:
            self.fieldTower = FE.fieldTower
        else:
            self.fieldTower = fieldTower
            
            
        """if not isNumber(numerator) and not isNumber(denominator):
            if self.numerator.field!=self.denominator.field:
                raise Exception()
        if not isNumber(numerator):
            if self.numerator.field!=self.field:
                raise Exception()
        if not isNumber(denominator):
            if self.denominator.field!=self.field:
                raise Exception()"""
        if self.denominator!=1 and self.numerator!=1 and self.numerator!=0:
            if (isNumber(self.numerator) or not self.numerator.isZero()) and not ( isNumber(self.numerator) or isNumber(self.denominator)):
                self.removeCommonFactors()
        
        
        
        # TODO: clean up if's
        if (isNumber(numerator) or numerator.isConstant()) and (isNumber(denominator) or denominator.isConstant()):
            if isNumber(numerator) and isNumber(denominator):
                self.numerator = numerator/denominator
                self.denominator = 1
            elif isNumber(numerator) and denominator.isConstant():
                #print(str(numerator),str(denominator),"a")
                self.numerator = Pol.Polynomial([numerator/denominator.getConstant()],field=field,fieldTower=self.getFieldTower())#(denominator/numerator).asRational().Inverse().numerator()
                self.denominator = 1
            elif isNumber(denominator) and numerator.isConstant():
                self.numerator = Pol.Polynomial([numerator.getConstant()/denominator],field=field,fieldTower=self.getFieldTower())
                self.denominator = 1
            elif numerator.isConstant() and denominator.isConstant():
                self.numerator = Pol.Polynomial([numerator.getConstant()/denominator.getConstant()],field=field,fieldTower=self.getFieldTower())
                self.denominator = 1
            
        if isNumber(self.denominator):
            self.numerator = self.numerator/self.denominator
            self.denominator = 1
        elif self.denominator.isConstant():
            self.numerator = self.numerator/self.denominator.getConstant()
            self.denominator = 1

        if output:
            print("B")
        self.MakeDenominatorMonic()
        if output:
            print("Z")
    def removeCommonFactors(self):
        """
        cancels common factors, so that: numerator/denominator = p/q with gcd(p,q)=1
        """
        gcd = Pol.PolyGCD(self.numerator, self.denominator)
        #print(gcd)
        if not gcd==1:
            self.numerator = Pol.PolyDiv(self.numerator, gcd)[0]
            self.denominator = Pol.PolyDiv(self.denominator,gcd)[0]
            #print(PolyDiv(self.denominator,gcd)[0])
            
    def isZero(self):
        if isNumber(self.numerator):
            return numberIsZero(self.numerator)
        return self.numerator.isZero()
    
    def getFieldTower(self):
        if self.fieldTower == None:
            self.fieldTower = FE.fieldTower
            if FE.fieldTower==None:
                raise Warning()
        return self.fieldTower
    
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
        return RationalFunction(dp*q+(-1)*p*dq,q*q,field=self.field,fieldTower=self.getFieldTower()) # (p/q)'=(p'q-pq')/(q^2)
    
    def PartialFraction(self, denomFactorization):
        
        self.MakeDenominatorMonic()
        return Pol.PartialFractionWithPowerFactorization(self.numerator, denomFactorization)
    
    def Inverse(self):# f->1/f
        return RationalFunction(self.denominator,self.numerator,field=self.field,fieldTower=self.getFieldTower())
    
    def MakeDenominatorMonic(self):
        if isNumber(self.denominator):
            return
        lcoeff = self.denominator.getLeadingCoefficient()
        if self.output:
            print("K")
        lcoeff_poly = Pol.Polynomial(coefficients=[lcoeff],field=self.field,fieldTower=self.getFieldTower())
        if self.output:
            print("N")
        if isNumber(self.numerator):
            self.numerator = Pol.Polynomial([self.numerator],field=self.field,fieldTower=self.getFieldTower())/lcoeff_poly
        else:
            self.numerator = self.numerator/lcoeff_poly
        if self.output:
            print("M")
            print(self.denominator,self.denominator.field)
            print(self.denominator.getCoefficients())
            print(lcoeff_poly,lcoeff_poly.field)
        self.denominator = self.denominator/lcoeff_poly
        if self.output:
            print("J")
        
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
        if isNumber(other):
            return RationalFunction(self.numerator+other*self.denominator,self.denominator,field=self.field,fieldTower=self.getFieldTower())
        if other.isZero():
            return self
        if type(other) == Pol.Polynomial or isPoly(other):
            return self.__add__(RationalFunction(other,1,field=self.field,fieldTower=self.getFieldTower()))
        num = self.numerator*other.denominator+self.denominator*other.numerator
        den = self.denominator*other.denominator
        return RationalFunction(num,den,field=self.field,fieldTower=self.getFieldTower())
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __mul__(self, other):
        if other == 1:
            return self
        if type(other)==Pol.Polynomial or isPoly(other):
            return self.__mul__(RationalFunction(other,1,field=self.field,fieldTower=self.getFieldTower()))
        if isNumber(other):
            num = self.numerator*other
            return RationalFunction(num,self.denominator, field=self.field,fieldTower=self.getFieldTower())
     #   if self.isConstant():
     #       return other*self.getConstant()
     #   if other.isConstant():
     #       return self*other.getConstant()
        num = self.numerator * other.numerator
        denom = self.denominator*other.denominator
        return RationalFunction(num,denom,field=self.field,fieldTower=self.getFieldTower())
    
    def __truediv__(self, other):
        if other==1:
            return self
        if isNumber(other):
            return self * (1/other)
        if type(other) == Pol.Polynomial or isPoly(other):#_isPoly(other):
            return self.__mul__(RationalFunction(1,other,field=self.field,fieldTower=self.getFieldTower()))
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
        if self.denominator==1:
            return "[{}]".format(numStr)
        return "[{}]/[{}]".format(numStr,denomStr)
'''
    
    
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