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


    def __init__(self, normCoeffs=None, princCoeffs=None, fieldTower=None):
        self.fieldTower = fieldTower
        if fieldTower==None or fieldTower.towerHeight==0:
            raise Exception("Extended Polynomials only for exponentials")
        fieldExtension = fieldTower.getLastExtension()
        if fieldExtension.extensionType!=FE.TRANS_EXP:
            raise Exception("Extended Polynomials only for exponentials")
        inverseExtension = FE.FieldExtension(FE.TRANS_EXP, -fieldExtension.argFunction,"{{1/{}}}".format(fieldExtension.variable))
        self.princFieldTower = self.fieldTower.copy()
        self.princFieldTower.addFieldExtension(inverseExtension)
        
        if princCoeffs!=None:
            princCoeffs = [0]+princCoeffs
        self.normPoly = Pol.Polynomial(normCoeffs, fieldTower=self.fieldTower)
        self.principalPart = Pol.Polynomial(princCoeffs, fieldTower=self.princFieldTower)
        
        
    @property
    def degree(self):
        return self.normPoly.degree
    @property
    def princDegree(self):
        return self.principalPart.degree
    
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

    def setCoefficient(self, coeff, power,callUpdates=True):
        if power>=0:
            self.normPoly.setCoefficient(coeff, power,callUpdates=callUpdates)
        else:
            self.principalPart.setCoefficient(coeff, -power,callUpdates=callUpdates)
            
    def updateCoefficientsAll(self):
        self.normPoly.updateCoefficientsAll()
        self.principalPart.updateCoefficientsAll()
            
    def mulByMonicMonomial(self,power): # multiplies self by x^power
        newExtPoly = ExtendedPolynomial(fieldTower=self.fieldTower)
        for i in range(self.degree,-self.princDegree-1,-1):
            newExtPoly.setCoefficient(self.getCoefficient(i), i+power,callUpdates=False)

        return newExtPoly
    
    def __radd__(self, other):
        return self.__add__(other)
    def __add__(self, other):
        if other==0:
            return self
        if type(other) == Pol.Polynomial:
            norm = self.normPoly+other
            pcoeffs = self.principalPart.getCoefficients()
            return ExtendedPolynomial(norm.getCoefficients(), pcoeffs[1:self.princDegree+1],fieldTower=self.fieldTower)
        if type(other) == ExtendedPolynomial:
            if other.fieldTower==self.fieldTower:
                return extPolyFromPolys(self.normPoly+other.normPoly, self.principalPart+other.principalPart, self.fieldTower)
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
    
def extPolyFromPolys(normPoly,princPoly,fieldTower):
    return ExtendedPolynomial(normPoly.getCoefficients(),princPoly.getCoefficients()[1:princPoly.degree+1],fieldTower=fieldTower)
def extPolyFromPol_power(poly,power):  # poly/T^power
    extPoly = ExtendedPolynomial(fieldTower=poly.fieldTower)
    for i in range(poly.degree,-1,-1):
        extPoly.setCoefficient(poly.getCoefficient(i), i-power, callUpdates=False)
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