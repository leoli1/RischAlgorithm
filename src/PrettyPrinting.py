'''
Created on 23.10.2018

@author: Leonard
'''
from __future__ import division
from Number import Rational
from PrettyPrintingHelper import *
from Utils import *
import FieldExtension as FE
import Polynomial as Pol

class StringMatrix(object):
    
    def __init__(self, width=1,height=1):
        self._data = [[' ' for _ in range(width)] for _ in range(height)]
        self.width = width
        self.height = height
        self.mainRow = 0
        self.hasBrackets = False
        
        """
            2
        541
        ----- <- mainRow
         2 
        """
        
    @staticmethod
    def fromString(s):
        l = len(s)
        sm = StringMatrix(width=l,height=1)
        for i in range(l):
            sm._data[0][i] = s[i]
        return sm
    
    def setRowAtCenter(self, row, data):
        data = list(data)
        d = len(data)
        total = self.width
        start = total//2-d//2
        self._data[row] = [self._data[row][i] if (i<start or i>=start+d) else data[i-start] for i in range(total)]
        
    def insertColumn(self, index, col):
        self.width += 1
        for i in range(self.height):
            self._data[i].insert(index,col[i])
    def insertEmptyRow(self, index):
        self.height+=1
        self._data.insert(index, [" "]*self.width)
        if index<=self.mainRow:
            self.mainRow += 1
    def setSubMatrix(self, matr, i,j):
        if matr.height>(self.height-i) or matr.width>(self.width-j):
            raise Exception()
        for k in range(i,i+matr.height):
            for l in range(j, j+matr.width):
                i1 = k-i
                j1 = l-j
                self._data[k][l] = matr._data[i1][j1]
              
        
    def __add__(self, other):
        #if self.height!=other.height:
         #   raise Exception()
        if self.width==1 and self.height==1 and self._data[0][0]==" ":
            return other
        
        mainA = self.mainRow
        mainB = other.mainRow
        if mainA<mainB:
            while self.mainRow<other.mainRow:
                self.insertEmptyRow(0)
        elif mainB<mainA:
            while other.mainRow< self.mainRow:
                other.insertEmptyRow(0)
        hA = self.height
        hB = other.height
        
        if hA<hB:
            while self.height<other.height:
                self.insertEmptyRow(self.height)
        elif hB<hA:
            while other.height<self.height:
                other.insertEmptyRow(other.height)
        
        
        dat = []
        for i in range(self.height):
            dat.append(self._data[i]+other._data[i])
        sm = StringMatrix(width=self.width+other.width,height=self.height)
        sm._data = dat
        sm.mainRow = self.mainRow
        return sm
    
    def __str__(self):
        s = u""
        for row in self._data:
            for c in row:
                s += c
            s += u"\n"
        return s.strip("\n")
    
    def encapsulateWithBrackets(self, bracketType="()"):
        if bracketType=="()":
            self.insertColumn(0, LeftBracket1(self.height))
            self.insertColumn(self.width, RightBracket1(self.height))
            self.hasBrackets = True
            
    def pprint(self):
        print(self.__str__())
        
    
    
    
def getStringMatrix(obj):
    if obj==None:
        return StringMatrix()
    if isNumber(obj) and type(obj)!=Rational:
        obj = Rational.fromFloat(obj)
    if type(obj)==Rational:
        return getStringMatrixFromRational(obj)
    
    if type(obj)==FE.Variable:
        fe = obj.fieldExtension
        if fe.argFunction==None:
            return StringMatrix.fromString(obj.stringRepr) 
        func = fe.getFVar()
        argSM = getStringMatrix(fe.argFunction)
        argSM.encapsulateWithBrackets()
        h = argSM.height//2
        funcSM = StringMatrix(width=len(func),height=argSM.height)
        funcSM.setRowAtCenter(h, func)
        return funcSM+argSM
    
    if type(obj)==Pol.Polynomial:
        return getStringMatrixFromPolynomial(obj)

def getStringMatrixFromRational(rational):
    p = rational._p
    q = rational._q
    if q==1:
        s = str(p)
        sm = StringMatrix(width=len(s),height=1)
        sm.setRowAtCenter(0, s)
        return sm
    ps = str(p)
    qs = str(q)
    w = max(len(ps),len(qs))
    sm = StringMatrix(width=w,height=3)
    sm.setRowAtCenter(0, ps)
    sm.setRowAtCenter(1,FRACTION_MID*w)
    sm.setRowAtCenter(2,qs)
    sm.mainRow = 1
    return sm
    
def getStringMatrixFromPolynomial(polynomial):
    sm = StringMatrix()
    
    for i in range(polynomial.degree+1):
        if polynomial.coeffIsZero(i):
            continue
        
        mSM = getStringMatrixFromMonomial(polynomial.getCoefficient(i), polynomial.variable, i)
        sm += mSM
        if i!=polynomial.degree:
            p = StringMatrix(width=1,height=1)
            p._data[0][0] = "+"
            sm += p
    return sm
        
def getStringMatrixFromMonomial(coeff, var, power):
    coeffSM = getStringMatrix(coeff)
    if power==0:
        return coeffSM
    if coeffSM.height>1 and not coeffSM.hasBrackets:
        coeffSM.encapsulateWithBrackets()
    var_pow_SM = getStringMatrixFromPower(var, power)
    if coeff==1:
        return var_pow_SM
    return coeffSM + var_pow_SM
    
def getStringMatrixFromPower(base, power):
    baseSM = getStringMatrix(base)
    if power==1:
        return baseSM
    if baseSM.height>1 and not baseSM.hasBrackets:
        baseSM.encapsulateWithBrackets()
    powerSM = getStringMatrix(power)
    if powerSM.height>1 and not powerSM.hasBrackets:
        powerSM.encapsulateWithBrackets()
    rowsNeeded = powerSM.height
    rowsToAdd = rowsNeeded-baseSM.mainRow
    sm = StringMatrix(height=baseSM.height+rowsToAdd,width=baseSM.width+powerSM.width)
    sm.setSubMatrix(baseSM, rowsToAdd, 0)
    sm.setSubMatrix(powerSM, 0, baseSM.width)
    sm.mainRow = rowsToAdd+baseSM.mainRow
    return sm
    

if __name__ == "__main__":
    [T1,T2] = FE.Variables(['T_1','T_2'])
    x = FE.getVariable("x")
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,Number.ONE]),T1) # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,Number.ONE],variable=T1),T2) # field extension with log(log(x))
    FE.fieldTower = FE.fieldTower.addFieldExtension(fieldExtension2)
    
    
    """sm = getStringMatrixFromRational(Rational(5,3))
    sm.pprint()
    sm = getStringMatrixFromRational(Rational(123,7))
    sm.encapsulateWithBrackets()
    sm.pprint()
    sm = getStringMatrixFromPower(Rational(5,3),45)
    sm.pprint()
    
    sm = getStringMatrixFromMonomial(Rational(100,3), x, Rational(42,5))
    sm.pprint()"""
    
    poly = Pol.Polynomial([0,Pol.Polynomial([0,1],x),Rational(5,2)],variable=T1)
    sm = getStringMatrixFromPolynomial(poly)
    sm.pprint()