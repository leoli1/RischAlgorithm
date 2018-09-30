'''
Created on 22.09.2018

@author: Leonard
'''
from __future__ import division
import FieldExtension as FE
import RationalFunction as Rat

from Utils import *
from FieldExtension import fieldTower

class Polynomial(object):
    
    def __init__(self, coefficients=None, fieldTower=None):
        """
        fieldTower = C(x,T1,T2,...,TN)
        polynomial in C(x,T1,T2,...T(N-1))[TN]
        -> coeffs in C(x,T1,T2,...T(N-1))
        (if n=0, then fieldTower=C(x),polynomial in C[x] -> coeffs in C)
        
        """
        
        if fieldTower==None:
            self.fieldTower = FE.FieldTower()
        else:
            self.fieldTower = fieldTower
            
        if coefficients == None:
            self._coefficients = [0]
        else:
            self._coefficients = coefficients
        
        self.updateCoefficients()
        self.updateCoefficientsFields()
        
        
    # ========================================== Coefficients stuff =========================================
    def setCoefficient(self, coeff, power):
        if power>self.degree:
            self._coefficients += [0]*(power-self.degree)
        self._coefficients[power] = coeff
        
        self.updateCoefficients()
        self.updateCoefficientsFields()
        
    def getCoefficient(self, power):
        if power>=len(self._coefficients):
            return 0
        return self._coefficients[power]
            
    def getLeadingCoefficient(self):
        return self.getCoefficient(self.degree)
    def setCoefficients(self,coeffs):
        if type(coeffs)!=list:
            raise Exception("coefficients parameter must be a list")
        if len(coeffs)<=self.degree:
            self._coefficients = self._coefficients[0:len(coeffs)]
        for i in range(len(coeffs)):
            self._coefficients[i] = coeffs[i]
        
        self.updateCoefficients()
        
    def getCoefficients(self):
        return self._coefficients
    
    def coeffIsZero(self, power):
        return objEqualsNumber(self.getCoefficient(power),0)
    
    def updateCoefficients(self):
        """
        Removes coefficients equal to Zero
        """
        newDeg = 0
        for i in range(len(self._coefficients)-1,-1,-1):
            if not self.coeffIsZero(i):
                newDeg = i
                break
        self._coefficients = self._coefficients[0:newDeg+1]
    def updateCoefficientsFields(self):   
        fieldTower = self.getFieldTower()
        for i in range(len(self._coefficients)):
            coeff = self.getCoefficient(i)
            coeff_tower = FE.FieldTower()
            if not isNumber(coeff):
                coeff_tower = coeff.getFieldTower().copy()
                if coeff.isZero():
                    continue
            if not isNumber(coeff) and not (fieldTower.isExtendedTowerOf(coeff_tower)):
                raise Exception("Coefficients of polynomial {} have to be in sub-towers of the polynomials field tower.".format(str(self)))
            
            newCoeff = coeff
            while coeff_tower.towerHeight< fieldTower.towerHeight-1 or (isNumber(newCoeff) and fieldTower.towerHeight>0):
                newExtension = fieldTower.getFieldExtension(coeff_tower.towerHeight)
                coeff_tower.addFieldExtension(newExtension)
                prev_tower = coeff_tower.getStrippedTower(coeff_tower.towerHeight-1)
                if not isNumber(newCoeff) and newCoeff.fieldTower==prev_tower:
                    pol = newCoeff
                else:
                    pol = Polynomial([newCoeff],fieldTower=prev_tower)
                if not isNumber(newCoeff):
                    newCoeff = Rat.RationalFunction(Polynomial([pol],fieldTower=coeff_tower),1,fieldTower=coeff_tower)#Rat.RationalFunction(Polynomial([newCoeff],fieldTower=prev_tower),1,fieldTower=prev_tower)
                else:
                    newCoeff = Rat.RationalFunction(pol,1,fieldTower=prev_tower) # numbers are 'below' C[x] Polynomials / C(x) Rational first, so they have to be converted to a C(x) first
                    coeff_tower = FE.FieldTower()
            self._coefficients[i] = newCoeff
                        
                        
    def isConstant(self):
        if self.deg0():
            if isNumber(self.getCoefficient(0)):
                return True
            else:
                return self.getCoefficient(0).isConstant()
        else:
            return False
        
    def getConstant(self):
        if not self.isConstant():
            return None
        if isNumber(self.getCoefficient(0)):
            return self.getCoefficient(0)
        return self.getCoefficient(0).getConstant()
    
    def getConstantCoefficients(self):
        if not self.hasOnlyConstantCoefficients():
            return None
        coeffs = []
        for c in self.getCoefficients():
            if (isNumber(c)):
                coeffs.append(c)
            else:
                coeffs.append(c.getConstant())
        return coeffs
    def hasOnlyConstantCoefficients(self):
        for coeff in self.getCoefficients():
            if isNumber(coeff):
                continue
            if not coeff.isConstant():
                return False
        return True
    
    def isZero(self):
        return objEqualsNumber(self.getConstant(),0)
    
    def deg0(self):
        return self.degree==0
    # ========================================== Field Tower/Extension =========================================
            
    def getFieldTower(self):
       # if self.fieldTower == None:
       #     self.fieldTower = FE.FieldTower()
       #     if FE.fieldTower==None:
       #         raise Warning()
        return self.fieldTower
    def getFieldExtension(self):
        return self.getFieldTower().getLastExtension()
    
    
    # ========================================== SquareFree stuff =========================================
    def isSquareFree(self): # polynomial p is square-free iff gcd(p,(d/dPolyVar)p)=1
        gcd = PolyGCD(self,self.differentiateWRTtoPolyVar())#PolyGCD(self,self.differentiate())
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
        d = self.differentiateWRTtoPolyVar()
        c = PolyGCD(self, d) # = a_2*(a_3)^2*(a_4)^3...
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
    # ========================================== Mixed stuff =========================================
    @property
    def degree(self):
        self.updateCoefficients()
        return len(self.getCoefficients())-1#-1 if self.isZero() else len(self.getCoefficients())-1
    
    def asRational(self):
        return Rat.RationalFunction(self,1,fieldTower = self.getFieldTower())
    
    def isMonic(self):
        return objEqualsNumber(self.getLeadingCoefficient(), 1)
    def makeMonic(self):
        if self.isZero():
            return self
        x = self.getLeadingCoefficient()
        return self/Polynomial([x],fieldTower=self.getFieldTower())
    def evaluate(self, val):
        if self.getFieldTower().towerHeight!=0:
            raise Exception()
        sum = 0
        for i in range(self.degree+1):
            sum += self.getCoefficient(i)*(val**i)
        return sum
    
    def getRoots(self):
        coeffs = self.getConstantCoefficients()
        if coeffs==None:
            return NotImplementedError()
        if self.degree>2:
            return NotImplementedError()
        
        if self.degree==1:
            return [-coeffs[0]/coeffs[1]]
        elif self.degree==2:
            raise NotImplementedError()
    # ========================================== Differentiate stuff =========================================
    def differentiate(self):
        fieldTower = self.getFieldTower()
        dPoly = Polynomial(fieldTower=fieldTower)
        for i in range(self.degree,-1,-1):
            p = self.getCoefficient(i)
            fieldExtension = self.getFieldExtension()
            if p==0:
                continue
            if fieldExtension==None:
                if i>0:
                    dPoly.setCoefficient(p*i, i-1)
            elif fieldExtension.extensionType==FE.TRANS_EXP: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'*T^i since T'=u'T
                dPoly += Monomial(i, p.differentiate(), fieldTower=fieldTower)
                dPoly += i*Monomial(i,p*fieldExtension.characteristicFunction.differentiate(),fieldTower=fieldTower)
            elif fieldExtension.extensionType==FE.TRANS_LOG: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'/u*T^(i-1)
                dPoly += Monomial(i,p.differentiate(),fieldTower=fieldTower)
                log_diff_u = Rat.RationalFunction(fieldExtension.characteristicFunction.differentiate(),fieldExtension.characteristicFunction,fieldTower=fieldTower.getStrippedTower(fieldTower.towerHeight-1))
                dPoly += i*Monomial(i-1, p* log_diff_u,fieldTower=fieldTower)
        return dPoly
    
    def differentiateWRTtoPolyVar(self):
        """
        differentiates function with respect to the field extension variable of the polynomial:
        p = 2*T**2+T+1
        -> differentiateWRTtoPolyVar(p) = (d/dT)p = 4T+1
        """
        dPoly = Polynomial(fieldTower=self.getFieldTower())
        for i in range(self.degree,0,-1):
            p = self.getCoefficient(i)
            dPoly.setCoefficient(p*i, i-1)
        return dPoly
    
    # ========================================== Arithmetic stuff =========================================
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other==0:
            return self
        if type(other)==Rat.RationalFunction:
            return other.__add__(self)
        if isNumber(other):
            deg0Part = self.getCoefficient(0)
            tPoly = Polynomial(fieldTower=self.getFieldTower())
            for i in range(self.degree,0,-1):
                tPoly.setCoefficient(self.getCoefficient(i), i)
            tPoly.setCoefficient(deg0Part+other,0)
            return tPoly
            #raise Exception(str(other))
        if (self.fieldTower!=other.fieldTower):
            #if self.field<other.field:
            #    return other+self.increaseField(other.field)
            #else:
            #    return self+other.increaseField(self.field)
            raise Exception("Polynomials have to have coefficients in the same field in order to add them")
        
        tDeg = max(self.degree,other.degree)
        tPoly = Polynomial(fieldTower=self.fieldTower)
        for i in range(tDeg,-1,-1):
            a = self.getCoefficient(i)
            b = other.getCoefficient(i)
            tPoly.setCoefficient(a+b, i)
        return tPoly
    
    def __rmul__(self, other):
        return self.__mul__(other)
    def __mul__(self, other):
        if other==1:
            return self
        if other==0:
            return 0
        if isNumber(other):
            tPoly = Polynomial(fieldTower=self.getFieldTower())
            for i in range(self.degree,-1,-1):
                tPoly.setCoefficient(self.getCoefficient(i)*other, i)
            return tPoly
            #raise Exception("??")
        if type(other) == Rat.RationalFunction:
            return other.__mul__(self)
        if (self.fieldTower!=other.fieldTower):
            raise Exception("Polynomials have to have coefficients in the same field in order to multiply them")
        
        tDeg = self.degree+other.degree
        tPoly = Polynomial(fieldTower=self.getFieldTower())
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
        if isNumber(other):
            return self.__mul__(1/other)
        (quot,rem)= PolyDiv(self, other)
        if rem==0:
            return quot
        
        return Rat.RationalFunction(self,other,field=self.field,fieldTower=self.getFieldTower())#Rat.RationalFunction(rem,other, field=self.field)+quot
        
    def __pow__(self, other):
        if not isNumber(other) or int(other)!=other or other<0:
            raise Exception("power has to be an integer >= 0")
        if other==0:
            return 1
        return (self**(other-1))*self
    
    # ========================================== String output =========================================
    def getVariable(self):
        return "x" if self.getFieldTower().towerHeight==0 else self.getFieldExtension().variable
    def __str__(self):
        return self.strCustomVar(self.getVariable())
    
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
                if isNumber(coeff_v) and coeff_v>0:
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant() and coeff_v.getConstant()>0:
                    if coeff_v.getConstant()!=1 or i==0:
                        coeff = str(coeff_v)
                else:
                    coeff = "({})".format(str(self.getCoefficient(i)))
                
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
    
    def printFull(self):
        """
        writes out field extension variables
        """
        out = ""
        fE = self.getFieldExtension()
        for i in range(len(self._coefficients)-1,-1,-1):
            if self.coeffIsZero(i):
                continue
            coeff = ""
            coeff_v = self.getCoefficient(i)
            if coeff_v != 1 or i==0:
                if isNumber(coeff_v) and coeff_v>0:
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant() and coeff_v.getConstant()>0:
                    if coeff_v.getConstant()!=1 or i==0:
                        coeff = str(coeff_v)
                else:
                    coeff = "({})".format(self.getCoefficient(i).printFull())

            power = ""
            if i>1:
                power = "**{}".format(i)
            elif i==1:
                power = ""
            var = ""
            if i>0:
                if fE==None:
                    var += "x"
                else:
                    var += fE.strFunc()
            
            out += coeff +"{}{}+".format(var,power)
        if len(out)==0:
            return "0"
        return out.strip("+")
'''class Polynomial(object):
    

    def __init__(self, coefficients = None, field=FE.BASE_FIELD,fieldTower=None):
        """
        field: the field(-extension) of its coefficients
        """
        self.field = field
        if coefficients == None:
            coefficients = [0]
        if type(coefficients)!=list:
            raise Exception()
        self._coefficients = coefficients
        
        self.degree = len(self._coefficients)-1
        self.updateDegree()
        
        if fieldTower==None:
            self.fieldTower = FE.fieldTower
        else:
            self.fieldTower = fieldTower

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
    def getCoefficients(self):
        return self._coefficients
    
    def getLeadingCoefficient(self):
        self.updateDegree()
        return self.getCoefficient(self.degree)
    
    def getFieldTower(self):
        if self.fieldTower == None:
            self.fieldTower = FE.fieldTower
            if FE.fieldTower==None:
                raise Warning()
        return self.fieldTower
    def getStrippedFieldTower(self):
        return self.getFieldTower().getStrippedTower(self.field)
    
    def getFieldExtension(self):
        if self.field==FE.BASE_FIELD:
            return None
        t = self.getStrippedFieldTower()
        return t.getFieldExtension(t.towerHeight-1)
    
    def increaseField(self, newField):
        if newField<=self.field:
            return self
        p = Polynomial([self],self.field+1)
        return p.increaseField(newField)
    def updateDegree(self):
        """
        tests if highest power terms have coefficients equal to 0 and eventually recalculates the degree
        """
        newDeg = 0
        for i in range(self.degree, 0, -1):
            if not self.coeffIsZero(i):
                newDeg = i
                break
        self.degree = newDeg
        self._coefficients = self._coefficients[0:newDeg+1]
        
    def coeffIsZero(self, power):
        c = self.getCoefficient(power)
        if isNumber(c):
            if numberIsZero(c):
                return True
            else:
                return False
        if c==0:
            return True
        return c.isZero()
    
    def isZero(self):
        for i in range(self.degree+1):
            if not self.coeffIsZero(i):
                return False
        return True
    
    def deg0(self):
        self.updateDegree()
        return self.degree==0
    
    def isConstant(self):
        if self.deg0():
            if isNumber(self.getCoefficient(0)):
                return True
            else:
                return self.getCoefficient(0).isConstant()
        else:
            return False
    def getConstant(self):
        if not self.isConstant():
            return None
        if isNumber(self.getCoefficient(0)):
            return self.getCoefficient(0)
        return self.getCoefficient(0).getConstant()
    def getConstantCoefficients(self):
        if not self.hasOnlyConstantCoefficients():
            return None
        coeffs = []
        for c in self.getCoefficients():
            if (isNumber(c)):
                coeffs.append(c)
            else:
                coeffs.append(c.getConstant())
        return coeffs
    def hasOnlyConstantCoefficients(self):
        for coeff in self.getCoefficients():
            if isNumber(coeff):
                continue
            if not coeff.isConstant():
                return False
        return True
    
    def getRoots(self):
        coeffs = self.getConstantCoefficients()
        if coeffs==None:
            return NotImplementedError()
        if self.degree>2:
            return NotImplementedError()
        
        if self.degree==1:
            return [-coeffs[0]/coeffs[1]]
        elif self.degree==2:
            raise NotImplementedError()
        
    def isSquareFree(self): # polynomial p is square-free iff gcd(p,p')=1
        gcd = PolyGCD(self,self.differentiateWRTtoPolyVar())#PolyGCD(self,self.differentiate())
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
        d = self.differentiateWRTtoPolyVar()
        c = PolyGCD(self, d) # = a_2*(a_3)^2*(a_4)^3...
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
        
    
    def evaluate(self, val):
        if self.field!=FE.BASE_FIELD:
            raise Exception()
        sum = 0
        for i in range(self.degree+1):
            sum += self.getCoefficient(i)*(val**i)
        return sum
    
    def getVariable(self):
        return "x" if self.field==FE.BASE_FIELD else self.getFieldExtension().variable
    
    def differentiate(self):
        dPoly = Polynomial(field = self.field,fieldTower=self.getFieldTower())
        for i in range(self.degree,-1,-1):
            p = self.getCoefficient(i)
            fieldExtension = self.getFieldExtension()
            if p==0:
                continue
            if self.field == FE.BASE_FIELD:
                if i>0:
                    dPoly.setCoefficient(p*i, i-1)
            elif fieldExtension.extensionType==FE.TRANS_EXP: # (p*T^i)'=p'T^i+p*i*T'*T^(i-1) = p'T^i+i*p*u'*T^i since T'=u'T
                dPoly += Monomial(i, p.differentiate(), field=self.field,fieldTower=self.getFieldTower())
                dPoly += i*Monomial(i,p*fieldExtension.characteristicFunction.differentiate(),field=self.field,fieldTower=self.getFieldTower())
            elif fieldExtension.extensionType==FE.TRANS_LOG:
                dp = 0 if isNumber(p) else p.differentiate()
                dPoly += Monomial(i,dp,field=self.field)
                if i>0:
                    log_diff_u = Rat.RationalFunction(fieldExtension.characteristicFunction.differentiate(),fieldExtension.characteristicFunction,field=self.field-1)
                    dPoly += i*Monomial(i-1,p*log_diff_u, field=self.field,fieldTower=self.getFieldTower())
        return dPoly
    
    def differentiateWRTtoPolyVar(self):
        """
        differentiates function with respect to the field extension variable of the polynomial:
        p = 2*T**2+T+1
        -> differentiateWRTtoPolyVar(p) = (d/dT)p = 4T+1
        """
        dPoly = Polynomial(field=self.field,fieldTower=self.getFieldTower())
        for i in range(self.degree,0,-1):
            p = self.getCoefficient(i)
            dPoly.setCoefficient(p*i, i-1)
        return dPoly
    
    def asRational(self):
        return Rat.RationalFunction(self,1,field=self.field,fieldTower=self.getFieldTower())
    
    def makeMonic(self):
        if self.isZero():
            return self
        x = self.getLeadingCoefficient()
        o = self/Polynomial([x],field=self.field,fieldTower=self.getFieldTower())
        return o
    
    def isMonic(self):
        lcoeff = self.getLeadingCoefficient()
        if isNumber(lcoeff):
            return lcoeff==1
        return lcoeff.isConstant() and lcoeff.getConstant()==1
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __add__(self, other):
        if other ==0:
            return self
        if isNumber(other):
            return self.__add__(Polynomial([other],field=self.field,fieldTower=self.getFieldTower()))
        if type(other)==Rat.RationalFunction:
            return other.__add__(self)
        if (self.field!=other.field):
            if self.field<other.field:
                return other+self.increaseField(other.field)
            else:
                return self+other.increaseField(self.field)
            #raise Exception("Polynomials have to have coefficients in the same field in order to add them")
        tDeg = max(self.degree,other.degree)
        tPoly = Polynomial(field=self.field,fieldTower=self.getFieldTower())
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
        if isNumber(other):
            return self.__mul__(Polynomial(coefficients=[other],field=self.field,fieldTower=self.getFieldTower()))
        elif type(other) == Rat.RationalFunction:
            return other.__mul__(self)
        if other.isZero() or self.isZero():
            return 0
        
        if (self.field!=other.field):
            if other.isConstant():
                return self*other.getConstant()
            if self.isConstant():
                return other*self.getConstant()
            if self.field<other.field:
                return other*self.increaseField(other.field)
            else:
                return self*other.increaseField(self.field)
            #raise Exception("Polynomials have to have coefficients in the same field in order to multiply them")
        tDeg = self.degree + other.degree
        tPoly = Polynomial(field=self.field,fieldTower=self.getFieldTower())
        # calcs new coefficients using Cauchy-Product
        for i in range(tDeg,-1,-1):
            newCoeff = 0
            
            for j in range(0,i+1):
                a = self.getCoefficient(j)
                b = other.getCoefficient(i-j)
                newCoeff += a*b
            tPoly.setCoefficient(newCoeff, i)
            
        return tPoly
    def __pow__(self, other):
        if not isNumber(other) or int(other)!=other or other<0:
            raise Exception("power has to be an integer >= 0")
        if other==0:
            return 1
        return (self**(other-1))*self
    
    def __rtruediv__(self, other):
        return other.__truediv__(self)
    
    def __truediv__(self, other):
        if isNumber(other):
            return self.__mul__(1/other)
        (quot,rem)= PolyDiv(self, other)
        if rem==0:
            return quot
        
        return Rat.RationalFunction(self,other,field=self.field,fieldTower=self.getFieldTower())#Rat.RationalFunction(rem,other, field=self.field)+quot
    
    def __str__(self):
        return self.strCustomVar(self.getVariable())
    
    def __repr__(self):
        return self.__str__()

    def strCustomVar(self, variable):
        out = ""#FE.VARIABLES[self.field]
        for i in range(self.degree,-1,-1):
            if self.coeffIsZero(i):
                continue
            coeff = ""
            coeff_v = self.getCoefficient(i)
            if coeff_v != 1 or i==0:
                if isNumber(coeff_v) and coeff_v>0:
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant() and coeff_v.getConstant()>0:
                    if coeff_v.getConstant()!=1 or i==0:
                        coeff = str(coeff_v)
                else:
                    coeff = "({})".format(str(self.getCoefficient(i)))
                
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
    
    def printFull(self):
        """
        writes out field extension variables
        """
        out = ""
        fE = self.getFieldExtension()
        for i in range(self.degree,-1,-1):
            if self.coeffIsZero(i):
                continue
            coeff = ""
            coeff_v = self.getCoefficient(i)
            if coeff_v != 1 or i==0:
                if isNumber(coeff_v) and coeff_v>0:
                    coeff = str(coeff_v)
                elif not isNumber(coeff_v) and coeff_v.isConstant() and coeff_v.getConstant()>0:
                    if coeff_v.getConstant()!=1 or i==0:
                        coeff = str(coeff_v)
                else:
                    coeff = "({})".format(self.getCoefficient(i).printFull())

            power = ""
            if i>1:
                power = "**{}".format(i)
            elif i==1:
                power = ""
            var = ""
            if i>0:
                if fE==None:
                    var += "x"
                else:
                    var += fE.strFunc()
            
            out += coeff +"{}{}+".format(var,power)
        if len(out)==0:
            return "0"
        return out.strip("+")
  '''
    
def polyEqual(A,B):
    a = A
    if isNumber(A) and isNumber(B):
        return A==B
    return (A+(-1)*B).isZero()
def PolyDiv(polA, polB):
    if (polA.fieldTower!=polB.fieldTower):
        raise Exception("Polynomials have to have coefficients in the same field in order to apply PolyDiv")
    if polA.deg0() and polB.deg0():
        if isNumber(polB.getCoefficient(0)):
            return (Polynomial([polA.getCoefficient(0)/polB.getCoefficient(0)],fieldTower=polA.fieldTower),0)#Rat.RationalFunction(polA.getCoefficient(0),polB.getCoefficient(0),polA.field-1)],field=polA.field), 0)
        else:
            if type(polB.getCoefficient(0))==Polynomial or isPoly(polB.getCoefficient(0)):
                return (Polynomial([polA.getCoefficient(0)*polB.getCoefficient(0).asRational().Inverse()],fieldTower=polA.getFieldTower()),0) 
            else:
                return (Polynomial([polA.getCoefficient(0)*polB.getCoefficient(0).Inverse()],fieldTower=polA.getFieldTower()),0) 
    
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
        return (monomial+quot,remainder)
    
def PolyGCD(polA,polB):
    if polB==0:
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
        return 1
    
def __extendedEuclid(p,q):
    """
    finds x,y with
    xp+yq=gcd(p,q)
    """
    (P,Q) = (p,q)# if p.degree>=q.degree else (q,p)
    (s,r) = PolyDiv(P,Q)
    rm = 0 if r==0 else r.makeMonic()
    if polyEqual(PolyGCD(p,q),rm):
   # if r.deg0():
        if r.getFieldTower().towerHeight==0:
            r_inv_poly = Polynomial([1/r.getLeadingCoefficient()])
            return (r_inv_poly,r_inv_poly*(-1)*s)
        else:
            r_inv_poly = Polynomial([r.getLeadingCoefficient().Inverse()],fieldTower=r.getFieldTower())#Rat.RationalFunction(1,r,field=r.field)
            return (r_inv_poly,r_inv_poly*(-1)*s)
    (x,y) = extendedEuclid(Q, r)
    return (y,x+(-1)*s*y)
    
def extendedEuclid(p,q):
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
def PartialFraction(f,p,q):
    """
    f/(pq)=x/p+y/q
    gcd(p,q)=1
    """
    if f.getFieldTower()!=p.getFieldTower() or f.getFieldTower()!=q.getFieldTower():
        raise Exception()
    if PolyGCD(p,q)!=1:
        raise Exception()
    (sig,tau) = extendedEuclidGenF(p, q, f)
    return (tau,sig)


def PartialFractionWithPowerFactorization(numerator,factorization):
    if len(factorization)==1:
        return PartialFractionPower(numerator, factorization[0][0], factorization[0][1])
    
    otherFactors = 1
    for i in range(1,len(factorization)):
        otherFactors *= factorization[i][0]**factorization[i][1]
        
    if factorization[0][0]!=1:
        (r1,r2) = PartialFraction(numerator, factorization[0][0], otherFactors)
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
    
def Monomial(degree, coeff,fieldTower=None):
    ft = fieldTower
    if ft==None:
        ft = FE.FieldTower()
    return Polynomial([0]*degree+[coeff], fieldTower=ft)

def printFactorization(fact):
    out = ""
    for i in range(len(fact)):
        if fact[i][0]!=1:
            out += "({})".format(fact[i][0])
            if fact[i][1]!=1:
                out += "**{}".format(fact[i][1])
            out += "*"
    return out.strip("*")
            
def printPartialFraction(pfrac):
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
    from Parse import parseField0PolyFromStr,parseField0RatFromStr
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Polynomial([0,1]),"T1") # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Polynomial([0,1],fieldTower=FE.fieldTower.copy()),"T2") # field extension with log(log(x))
    FE.fieldTower.addFieldExtension(fieldExtension2)
    
    print("Field Tower: {}".format(FE.fieldTower))
    
    polA = Polynomial([1,1])
    print("polytest: should be x+1, got {}".format(polA))
    polB = Polynomial([1,1],fieldTower=FE.fieldTower.getStrippedTower(1))
    print("polytest: should be T1+1, got {}".format(polB))
    print("polytest: should be rationalfunction, got {}".format(type(polB.getCoefficient(0))))
    print("polytest: should be 0, got {}".format(polB.getCoefficient(0).getFieldTower().towerHeight))
    ratA = Rat.RationalFunction(1,polA)
    print("polytest: should be 1/[x+1], got {}".format(ratA))
    
    polD = Polynomial([ratA,3],fieldTower=FE.fieldTower.getStrippedTower(1))
    print("polytest: should be 3T1+(1/[x+1]), got {}".format(polD))
    print("polytest: should be 2x 0, got {}, {}".format(polD.getCoefficient(0).getFieldTower().towerHeight,polD.getCoefficient(1).getFieldTower().towerHeight))
    
    polC = Polynomial([ratA,1], fieldTower=FE.fieldTower)
    print("polytest: should be T2+(1/[x+1]), got {}".format(polC))
    print("polytest: should be 2x RationalFunction, got {},{}".format(type(polC.getCoefficient(0)),type(polC.getCoefficient(1))))
    print("polytest: should be 2, got {}".format(polC.getFieldTower().towerHeight))
    print("polytest: should be 2x 1, got {},{}".format(polC.getCoefficient(0).getFieldTower().towerHeight,polC.getCoefficient(1).getFieldTower().towerHeight))
    
    print("======================== Arithmetic stuff tests ==============================")
    polE = polB+polD
    print("polytest: should be 4T1+([x+2]/[x+1]), got {}".format(polE))
    polF = polE+polD
    print("polytest: should be 7T1+([x+3]/[x+1]), got {}".format(polF))


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
    
    polF_T = Polynomial(coefficients=[polA,polB],fieldTower=FE.fieldTower.getStrippedTower(1)) # polB*(T_1)**1+polA = 
    print(polytest("(2x**2+x)T1+(4x**3+2x**2+1)",str(polF_T)))#polF_T)
    polG_T = polF_T*polF_T
    print(polytest("(4x**4+4x**3+x**2)T1**2+(16x**5+16x**4+4x**3+4x**2+2x)T1+(16x**6+16x**5+4x**4+8x**3+4x**2+1)",polG_T))#print("[{}] * [{}] = {}".format(polF_T,polF_T,polG_T))
    
    polH = Polynomial(coefficients=[0,Polynomial([1])], fieldTower=FE.fieldTower.getStrippedTower(1))
    polG = Polynomial(coefficients=[polH,polH],fieldTower=FE.fieldTower.getStrippedTower(2))
    print(polytest("T1",polH))
    print(polytest("(T1)T2+(T1)",polG))
    print(polytest("(log(x))log(log(x))+(log(x))",polG.printFull()))
    
    ratA = parseField0RatFromStr("1/x")
    ratB = parseField0RatFromStr("x**2/x+1")
    polA = Polynomial(coefficients=[ratA,ratB],fieldTower=FE.fieldTower.getStrippedTower(1)) # x^2/(x+1)*T1+1/x
    print(polytest("T1+([x+1]/[x**3])",polA.makeMonic()))
    
    #polA = Polynomial([0,1]) # x in C(x)
    #polB = Polynomial([0,1],fieldTower=FE.fieldTower.getStrippedTower(1))# T in C(x,T)
    #print(polA,polB)
    #polC = polA.increaseField(1)
    #print(polC,polC.field)
    #print(polA*polB)
    print("---------- Diff poly tests") #################################################################################
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
    polA = Polynomial(coefficients=[parseField0PolyFromStr("(-1)*x**2"),0,1],fieldTower=FE.fieldTower.getStrippedTower(1))
    polB = Polynomial(coefficients=[parseField0PolyFromStr("x"),1],fieldTower=FE.fieldTower.getStrippedTower(1))
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
    polA = Polynomial([ratB,0,ratA],fieldTower=FE.fieldTower.getStrippedTower(1)) # = -x^2+T^2=(T-x)(T+x)
    ratC = parseField0RatFromStr("x**2+1/x")
    polB = Polynomial([ratA,ratC,ratA],fieldTower=FE.fieldTower.getStrippedTower(1))# = T^2+(x^2+1)/x*T+1=(T+x)(T+1/x)
    print(polytest("T1+(x)",PolyGCD(polA, polB)))
    
    (x,y) = extendedEuclid(polA,polB)
    print(polytest("([x]/[(-1.0)x**2+(-1.0)]), ([(-1.0)x]/[(-1.0)x**2+(-1.0)])",(x,y)))
    #print("extended euclid: {}, {}, {}".format(x,y,PolyGCD(polA, polB)))
    f = PolyGCD(polA, polB)*Polynomial([Polynomial([1,1]),1],fieldTower=FE.fieldTower.getStrippedTower(1)) #*(1+x+T+T^2)
    (x,y) = extendedEuclidGenF(polA, polB, f)
    print(polytest("([(-1.0)x**2+(-1.0)x+1.0]/[x**2+1.0]), ([2.0x**2+x]/[x**2+1.0])",(x,y)))

    print(polytest("True",polA.isSquareFree()))
    
    polC = Polynomial([Polynomial([0,1]).asRational(),ratA],fieldTower=FE.fieldTower.getStrippedTower(1))
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
    polA = Polynomial([ratA,1],fieldTower=FE.fieldTower.getStrippedTower(1))**2
    print(polytest("False",polA.isSquareFree()))
    #print("{} is squarefree: {}, should be false".format(polA,polA.isSquareFree()))
    fact = polA.factorSquareFree()
    print(polytest("(T1+([x+1.0]/[x+2.0]))**2",printFactorization(fact)))
    #print("{} = {}".format(polA,printFactorization(fact)))
    assert polA.isSquareFree()==False
    
    ratB = parseField0RatFromStr("(-1)*x/1")
    polB = Polynomial([ratB,1],fieldTower=FE.fieldTower.getStrippedTower(1))*polA
    print(polytest("False",polB.isSquareFree()))
    #print("{} is squarefree: {}, should be false".format(polB,polB.isSquareFree()))
    fact = polB.factorSquareFree()
    print(polytest("(T1+((-1.0)x))*(T1+([x+1.0]/[x+2.0]))**2",printFactorization(fact)))
    #print("{} = {}".format(polB,printFactorization(fact)))
    assert polB.isSquareFree()==False
    
    
    # http://www.wolframalpha.com/input/?i=expand+(T-x)(T%2B(x%2B1)%2F(x%2B2))%5E2+at+T%3D0
    ratA = parseField0RatFromStr("x+2*x**2+x**3/4+4*x+x**2")
    ratB = parseField0RatFromStr("1+(-2)*x+(-5)*x**2+(-2)*x**3/4+4*x+x**2")
    ratC = parseField0RatFromStr("2+(-1)*x**2/2+x")
    ratD = parseField0RatFromStr("1/1")
    
    polA = Polynomial([ratA*(-1),ratB,ratC,1],fieldTower=FE.fieldTower.getStrippedTower(1))
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
    polA = Polynomial([(-1)*x,one],fieldTower=FE.fieldTower.getStrippedTower(1)) # T-x
    polB = Polynomial([x,one],fieldTower=FE.fieldTower.getStrippedTower(1)) # T+x
    
    num = Polynomial([one,one],fieldTower=FE.fieldTower.getStrippedTower(1))#1+T
    denom = polA*(polB**2) # (T-x)(T+x)^2
    
    factorization = denom.factorSquareFree()
    #print("{} = {}".format(denom,printFactorization(factorization)))
    print(polytest("(T1+((-1.0)x))*(T1+(x))**2",printFactorization(factorization)))
    pfrac = PartialFractionWithPowerFactorization(num, factorization)
    print(polytest("[([0.25x+0.25]/[x**2])]/[T1+((-1.0)x)]+[([(-0.25)x+(-0.25)]/[x**2])]/[T1+(x)]+[([0.5x+(-0.5)]/[x])]/[T1+(x)]**2",printPartialFraction(pfrac)))
   # print("[{}]/[{}] = {}".format(num,denom,printPartialFraction(pfrac)))
        