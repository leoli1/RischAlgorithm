'''
Created on 08.10.2018

@author: Leonard
'''
from __future__ import division,print_function

import FieldExtension as FE
import Number
from Utils import *
import Polynomial as Pol
import RationalFunction as Rat
from Tkconstants import EXTENDED



class MultiIndex(object):
    """
    stores data of the type (a,n) where a is a variable and n an integer
    used for storing exponents data in multivariate polynomials
    """
    def __init__(self,data=None):
        self.indices = []
        if data!=None:
            self.indices = data
            if type(data)!=list:
                raise TypeError("data attribute has to be a list with 2-tuple elements.")

    def __neg__(self):
        ninds = []
        for index in self.indices:
            ninds.append((index[0],-index[1]))
        return MultiIndex(data=ninds)
    
    def __add__(self, other):
        if type(other)!=MultiIndex:
            raise TypeError("Can only add two MultiIndex-instances")
        is1 = self.indices
        is2 = other.indices
        #new_is = []
        variables = []
        indices = []
        for index in is1+is2:
            try:
                i = variables.index(index[0])
                indices[i] += index[1]
            except ValueError: # new index, not in variables
                variables.append(index[0])
                indices.append(index[1])
                
        return MultiIndex(data=zip(variables,indices))
    
    def remove0Powers(self):
        newInds = []
        for ind in self.indices:
            if ind[1]!=0:
                newInds.append(ind)
                
        self.indices = newInds
    
    def sortTerms(self):
        self.indices.sort(key=lambda x: str(x[0]))
    def __eq__(self, other):
        """
        TODO: Decide whether '0' exponents should be deleted from the indices array
        """
        if type(other)!=MultiIndex:
            return False
        if len(self)!=len(other):
            return False
        self.sortTerms()
        other.sortTerms()
        for i in range(len(self)):
            if self[i]!=other[i]:
                return False
        return True
    def __ne__(self, other):
        return not self.__eq__(other)
        
    def __len__(self):
        return len(self.indices)
    def __getitem__(self, index):
        return self.indices[index]
        
    def powerNotation(self):
        out = ""
        for ind in self.indices:
            if ind[1]==0:
                continue
            exp = ""
            if ind[1]>1:
                exp = "**"+str(ind[1])
            elif ind[1]<0:
                exp = "**({})".format(str(ind[1]))
            out += "({})".format(str(ind[0])+exp)
        return out
    def __str__(self):
        return self.powerNotation()
    def __repr__(self):
        return self.__str__()
                

class ExtendedMultivariatePolynomial(object):
    """
    stores data for a multivariate polynomial that can have negative powers, e.g. 5x*y**3*z + (-3)*x**(-1)*y
    """
    def __init__(self, terms=None):
        """
        terms = list of ExtendedMultivariateMonomials
        """
        self.terms = [] if terms == None else terms
        
        self.updateTerms()
        
    def updateTerms(self):
        newTerms = []
        for term in self.terms:
            if term.coefficient==0:
                continue
            term.updateMultiIndex()
            newTerms.append(term)
        #newTerms.sort(key=lambda x: x.multiIndex.)
        if len(newTerms)==0:
            newTerms = [ExtendedMultivariateMonomial(coefficient=Number.ZERO,multiIndex=MultiIndex())]
        self.terms = newTerms
    
    @staticmethod
    def fromNormalPolynomial(polynomial):
        if isNumber(polynomial):
            return ExtendedMultivariateMonomial(coefficient=polynomial,multiIndex=MultiIndex()).asPolynomial()
        elif polynomial.isConstant():
            return ExtendedMultivariateMonomial(coefficient=polynomial.getConstant(),multiIndex=MultiIndex()).asPolynomial()
        if type(polynomial) == Rat.RationalFunction:
            num = ExtendedMultivariatePolynomial.fromNormalPolynomial(polynomial.numerator)
            denom = ExtendedMultivariatePolynomial.fromNormalPolynomial(polynomial.denominator)
            if len(denom.terms)!=1:
                return None
            if id(num)==None:
                return None
            return num*denom.asMonomial().Inverse()
        
        output = ExtendedMultivariatePolynomial()
        for i in range(polynomial.degree+1):
            coeff = polynomial.getCoefficient(i)
            polyVarIndex = MultiIndex(data=[(polynomial.variable,i)])
            polyVarMonomial = ExtendedMultivariateMonomial(coefficient=Number.ONE,multiIndex=polyVarIndex)
            
            if isNumber(coeff):
                output += polyVarMonomial*coeff
            else:
                if type(coeff)==Pol.Polynomial:
                    output += ExtendedMultivariatePolynomial.fromNormalPolynomial(coeff)*polyVarMonomial
                elif type(coeff)==Rat.RationalFunction:
                    o = ExtendedMultivariatePolynomial.fromNormalPolynomial(coeff)
                    if o!=None:
                        output += o*polyVarMonomial#num*denom.asMonomial().Inverse()*polyVarMonomial
                    else:
                        return None
                else:
                    raise Exception()
                
        return output
        
    def asMonomial(self):
        if len(self.terms)!=1:
            return None
        return self.terms[0]
    def __add__(self, other):
        if type(other)==ExtendedMultivariateMonomial:
            return self.__add__(other.asPolynomial())
        if type(other)!=ExtendedMultivariatePolynomial:
            raise TypeError()
        
        newMultiIndices = []
        newCoefficients = []
        for term in self.terms+other.terms:    
            try:
                ind = newMultiIndices.index(term.multiIndex)
                newCoefficients[ind] += term.coefficient
            except ValueError:
                newMultiIndices.append(term.multiIndex)
                newCoefficients.append(term.coefficient)
        
        newTerms = []
        for i in range(len(newMultiIndices)):
            newTerms.append(ExtendedMultivariateMonomial(coefficient=newCoefficients[i],multiIndex=newMultiIndices[i]))
        return ExtendedMultivariatePolynomial(terms=newTerms)
    
    def __mul__(self, other):
        if type(other)==ExtendedMultivariateMonomial:
            return self.__mul__(other.asPolynomial())
        if type(other)!=ExtendedMultivariatePolynomial:
            raise TypeError()
        
        result = ExtendedMultivariatePolynomial()
        for t1 in self.terms:
            for t2 in other.terms:
                result += t1*t2
                
        return result
    
    def __str__(self):
        out = ""
        for term in self.terms:
            out += str(term)+"+"
        return out.strip("+")
    def __repr__(self):
        return self.__str__()
    
class ExtendedMultivariateMonomial(object):
    """
    stores data of a multivariate monomial that can have negative powers
    """
    def __init__(self, coefficient=None,multiIndex=None):
        if id(coefficient)==id(None):
            coefficient = Number.ONE
            
        self.coefficient = coefficient
        
        if id(multiIndex)==id(None):
            self.multiIndex = MultiIndex()
        else:
            self.multiIndex = multiIndex
            
        self.updateMultiIndex()
            
    def updateMultiIndex(self):
        """
        removes 0-powers and sorts the multiindex
        """
        self.multiIndex.remove0Powers()
        self.multiIndex.sortTerms()
    
    def asPolynomial(self):
        return ExtendedMultivariatePolynomial(terms=[self])
    
    def Inverse(self):
        return ExtendedMultivariateMonomial(coefficient=1/self.coefficient,multiIndex=-self.multiIndex)
    def __add__(self, other):
        if type(other)==ExtendedMultivariatePolynomial:
            return other.__add__(self.asPolynomial())
        
        if type(other)!=ExtendedMultivariateMonomial:
            raise TypeError()
        if other.multiIndex!=self.multiIndex:
            raise Exception()
        return ExtendedMultivariateMonomial(coefficient=self.coefficient+other.coefficient,multiIndex=self.multiIndex)
    def __mul__(self, other):
        if type(other)==ExtendedMultivariatePolynomial:
            return other.__mul__(self.asPolynomial())
        
        if isNumber(other):
            return ExtendedMultivariateMonomial(coefficient=self.coefficient*other,multiIndex=self.multiIndex)
        
        if type(other)!=ExtendedMultivariateMonomial:
            raise TypeError()
        return ExtendedMultivariateMonomial(coefficient=self.coefficient*other.coefficient,multiIndex=(self.multiIndex+other.multiIndex))
    def __str__(self):
        return str(self.coefficient)+"*"+str(self.multiIndex)
    def __repr__(self):
        return self.__str__()
        
        
if __name__ == '__main__':
    from Parse import *
    
    [a,b,c] = FE.Variables(['a','b','c'])
    mi1 = MultiIndex(data = [(a,5),(b,4)])
    mi2 = MultiIndex(data=[(b,-1)])
    print(mi1,mi2,mi1+mi2,-mi2,mi1==MultiIndex(data = [(a,5),(b,4)]),sep="\n")
    mon1 = ExtendedMultivariateMonomial(coefficient=Number.Rational(42,1), multiIndex=mi1)
    mon2 = ExtendedMultivariateMonomial(coefficient=Number.Rational(-1,1), multiIndex=mi2)
    pol1 = ExtendedMultivariatePolynomial(terms=[mon1,mon2])
    pol2 = ExtendedMultivariatePolynomial(terms=[ExtendedMultivariateMonomial(coefficient=3,multiIndex=mi2)])
    print(mon1,mon2,pol1,pol2,pol1+pol2,pol1*pol2, sep="\n")
    

    fe1 = FE.FieldExtension(FE.TRANS_EXP, Pol.Polynomial([0,Number.ONE]), FE.Variable('T_1'))
    ft = FE.FieldTower(fe1)
    fe2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,Number.ONE],fieldTower=ft),FE.Variable('T_2'))
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fe1,fe2])
    
    pol = parseExpressionFromStr("x*T_2**3+(x*T_1/T_2)+5*T_2/x", FE.fieldTower)
    emp = ExtendedMultivariatePolynomial.fromNormalPolynomial(pol)
    print(emp)
    