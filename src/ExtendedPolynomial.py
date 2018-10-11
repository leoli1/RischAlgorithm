'''
Created on 07.10.2018

@author: Leonard
'''
import FieldExtension as FE
import Polynomial as Pol
import Number

class ExtendedPolynomial(object):
    """
    stores data of an extended polynomial: a_n * x^n + a_(n-1) * x^(n-1) + ... + a_0 + a_-1*x^(-1)+...a_(-m)*x^(-m) where x is exponential
    """


    def __init__(self, normCoeffs=None, princCoeffs=None, variable=None):
        self.variable = variable
        if variable==None:
            self.variable = FE.BASEVARIABLE
        elif type(variable)!=FE.Variable:
            raise TypeError()
        fieldTower = self.variable.fieldExtension.fieldTower
        if fieldTower==None or fieldTower.towerHeight==1:
            raise Exception("Extended Polynomials only for exponentials")
        fieldExtension = fieldTower.getLastExtension()
        if fieldExtension.extensionType!=FE.TRANS_EXP:
            raise Exception("Extended Polynomials only for exponentials")
        inverseExtension = FE.FieldExtension(FE.TRANS_EXP, -fieldExtension.argFunction,FE.Variable("{{1/{}}}".format(fieldExtension.variable)))
        self.princFieldTower = self.fieldTower.copy()
        self.princFieldTower = self.princFieldTower.addFieldExtension(inverseExtension)
        
        if princCoeffs!=None:
            princCoeffs = [0]+princCoeffs
        self.normPoly = Pol.Polynomial(normCoeffs, variable=self.variable)
        self.principalPart = Pol.Polynomial(princCoeffs, variable=self.variable)
        
        
    @property
    def degree(self):
        return self.normPoly.degree
    @property
    def princDegree(self):
        return self.principalPart.degree
    @property
    def fieldTower(self):
        return self.variable.fieldExtension.fieldTower
    
    def isZero(self):
        return self.normPoly.isZero() and self.principalPart.isZero()
    def asPolynomial(self):
        if self.principalPart.isZero():
            return self.normPoly
        return None
    
    def getCoefficient(self, power):
        if power>=0:
            return self.normPoly.getCoefficient(power)
        else:
            return self.principalPart.getCoefficient(-power)

    def setCoefficient(self, power, coeff,callUpdates=True):
        if type(power)!=int: raise TypeError()
        if power>=0:
            self.normPoly.setCoefficient(power,coeff,callUpdates=callUpdates)
        else:
            self.principalPart.setCoefficient(-power,coeff,callUpdates=callUpdates)
            
    def updateCoefficientsAll(self):
        self.normPoly.updateCoefficientsAll()
        self.principalPart.updateCoefficientsAll()
            
    def mulByMonicMonomial(self,power): # multiplies self by x^power
        newExtPoly = ExtendedPolynomial(variable=self.variable)
        for i in range(self.degree,-self.princDegree-1,-1):
            newExtPoly.setCoefficient(i+power, self.getCoefficient(i),callUpdates=False)

        return newExtPoly
    
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other==0:
            return self
        if type(other) == Pol.Polynomial:
            norm = self.normPoly+other
            pcoeffs = self.principalPart.getCoefficients()
            return ExtendedPolynomial(norm.getCoefficients(), pcoeffs[1:self.princDegree+1],variable=self.variable)
        if type(other) == ExtendedPolynomial:
            if other.fieldTower==self.fieldTower:
                return extPolyFromPolys(self.normPoly+other.normPoly, self.principalPart+other.principalPart, self.variable)
            else:
                raise Exception()
            
        raise NotImplementedError()
        
    def __str__(self):
        n = str(self.normPoly)
        p = str(self.principalPart)
        return "{}+{}".format(n,p)
    def __repr__(self):
        return self.__str__()
    def printFull(self):
        n = self.normPoly.printFull()
        p = self.principalPart.printFull(reverse=True)
        if self.normPoly.isZero():
            return p
        if self.principalPart.isZero():
            return n
        return "{}+{}".format(n,p)
    
def extPolyFromPolys(normPoly,princPoly,variable):
    return ExtendedPolynomial(normPoly.getCoefficients(),princPoly.getCoefficients()[1:princPoly.degree+1],variable=variable)
def extPolyFromPol_power(poly,power):  # poly/T^power
    extPoly = ExtendedPolynomial(variable=poly.variable)
    for i in range(poly.degree,-1,-1):
        extPoly.setCoefficient(i-power, poly.getCoefficient(i), callUpdates=False)
    return extPoly
if __name__ == '__main__':
    o = Number.Rational(1,1)
    fieldTower = FE.FieldTower(FE.FieldExtension(FE.TRANS_EXP, Pol.Polynomial([0,o]),"T_1"))
    pol = ExtendedPolynomial([o,2*o],[-o,42*o], fieldTower) # 2T_1 + 1 + -1/T_1 + 42*1/T_1**2
    print(pol.printFull())
    pol2 = pol.mulByMonicMonomial(1) #2T_1**2 + 1T_1 + -1 + 42*1/T_1
    print(pol2.printFull())
    pol3 = pol.mulByMonicMonomial(2)
    print(pol3.printFull())