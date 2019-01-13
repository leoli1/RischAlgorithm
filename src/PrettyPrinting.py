'''
Created on 23.10.2018

@author: Leonard
'''
from __future__ import division
from Number import Rational,SqrRootPolynomial
from PrettyPrintingHelper import *
from Utils import *
import FieldExtension as FE
import Polynomial as Pol
import RationalFunction as Rat
import Integral as Int



class StringMatrix(object): 
    def __init__(self, width=1,height=1,data=None):
        if data!=None:
            self._data = data
            self.width = len(self._data[0])
            self.height = len(self._data)
        else:
            self._data = [[' ' for _ in range(width)] for _ in range(height)]
            self.width = width
            self.height = height
        self.mainRow = 0
        
        self.hasBrackets = False

    @staticmethod
    def fromString(s):
        l = len(s)
        sm = StringMatrix(width=l,height=1)
        for i in range(l):
            sm._data[0][i] = s[i]
        return sm
    
    def setChar(self, row,column, char):
        self._data[row][column] = char
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
       
    def replace(self, old, new):
        for x in range(self.width):
            for y in range(self.height):
                if self._data[y][x] == old:
                    self._data[y][x] = new     
    def __add__(self, other):
        if other==None:
            return self
        if self.width==1 and self.height==1 and self._data[0][0]==" ":
            return other
        
        while self.mainRow<other.mainRow:
            self.insertEmptyRow(0)
        while other.mainRow< self.mainRow:
            other.insertEmptyRow(0)

        while self.height<other.height:
            self.insertEmptyRow(self.height)
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
    
    def surroundWithBrackets(self, bracketType="()"):
        if bracketType=="()":
            self.insertColumn(0, LeftBracket1(self.height))
            self.insertColumn(self.width, RightBracket1(self.height))
            self.hasBrackets = True
            
    def pprint(self):
        print(self.__str__())
        
def pp(obj):
    if obj==None:
        return StringMatrix()
    if isNumber(obj) and type(obj)!=Rational and type(obj)!=SqrRootPolynomial:
        obj = Rational.fromFloat(obj)
    if type(obj)==SqrRootPolynomial:
        return ppSqrRootPolynomial(obj)
    if type(obj)==Rational:
        return ppRational(obj)
    
    if type(obj)==FE.Variable:
        fe = obj.fieldExtension
        if fe.argFunction==None:
            return StringMatrix.fromString(obj.stringRepr) 
        func = fe.getFVar()
        argSM = pp(fe.argFunction)
        argSM.surroundWithBrackets()
        h = argSM.height//2
        funcSM = StringMatrix(width=len(func),height=argSM.height)
        funcSM.setRowAtCenter(h, func)
        funcSM.mainRow = h
        return funcSM+argSM
    
    if type(obj)==Pol.Polynomial:
        return ppPolynomial(obj)
    if type(obj)==Rat.RationalFunction:
        return ppRationalFunction(obj)
    if type(obj)==str:
        return StringMatrix.fromString(obj)
    
def ppRootSums(rootSums):
    if len(rootSums)==0:
        return None
    sm = StringMatrix()
    for rs in rootSums:
        sm += ppRootSum(rs)
    return sm

def ppRootSum(rootSum):
    #rootSum.poly.variable = FE.BASEVARIABLE
    polySM = ppPolynomial(rootSum.poly) + StringMatrix.fromString("=0")
    polySM.replace(rootSum.poly.variable.stringRepr,"w")
    w = polySM.width
    sigma = bigSigma()
    sm = StringMatrix(height=sigma.height+polySM.height,width=w+1)
    sm.setSubMatrix(sigma,0,w//2-1)
    sm.setSubMatrix(polySM,2,0)
    sm.mainRow = sigma.height/2
    if type(rootSum.expr)==Int.LogFunction:
        argSM = pp(rootSum.expr.argFunction)
        argSM.surroundWithBrackets()
        exprSM = pp(rootSum.expr.factor) + StringMatrix.fromString("log")+argSM
        exprSM.replace(rootSum.exprVar,"w")
    else:
        exprSM = pp(rootSum.expr)
        exprSM.replace(rootSum.exprVar,"w")
    return sm + exprSM
    
    

def ppRational(rational):
    p = rational._p
    q = rational._q
    if q==1:
        sm = StringMatrix.fromString(str(p))
        if p<0:
            sm.surroundWithBrackets()
        return sm
    ps = str(p)
    qs = str(q)
    width = max(len(ps),len(qs))
    sm = StringMatrix(width=width,height=3)
    sm.setRowAtCenter(0, ps)
    sm.setRowAtCenter(1,FRACTION_MID*width)
    sm.setRowAtCenter(2,qs)
    sm.mainRow = 1
    if rational<0:
        sm.surroundWithBrackets()
    return sm

def ppRationalFunction(rational):
    numSM = pp(rational.numerator)
    if rational.denominator==1:
        return numSM
    denomSM = pp(rational.denominator)
    h1 = numSM.height
    h2 = denomSM.height
    width = max(numSM.width,denomSM.width)
    sm = StringMatrix(width=width,height=1+h1+h2)
    sm.setSubMatrix(numSM,0,(width-numSM.width)//2)
    sm.setRowAtCenter(h1,FRACTION_MID*width)
    sm.setSubMatrix(denomSM,h1+1,(width-denomSM.width)//2)
    sm.mainRow = h1
    return sm
def ppPolynomial(polynomial,customVar=None):
    sm = StringMatrix()
    
    p = StringMatrix(width=1,height=1)
    p._data[0][0] = "+"
    m = StringMatrix(width=1,height=1)
    m._data[0][0] = "-"
            
    for i in range(polynomial.degree+1):
        if polynomial.coeffIsZero(i):
            continue
        mSM = ppMonomial(polynomial.getCoefficient(i), polynomial.variable, i,customVar=customVar)
        sm += mSM
        if i!=polynomial.degree:
            sm += p
    return sm

def ppSqrRootPolynomial(rootPol):
    a = rootPol.a
    b = rootPol.b
    aSM = pp(a)
    bSM = pp(b)
    rootSM = ppSquareRoot(rootPol.radicand)
    p = StringMatrix(width=1,height=1)
    p._data[0][0] = "+"
    m = StringMatrix(width=1,height=1)
    m._data[0][0] = "-"
    if a==0:
        if b==1:
            return rootSM
        if b==-1:
            return m+rootSM
        return bSM+rootSM
    if b==1:
        return aSM+p+rootSM
    if b==-1:
        return aSM+m+rootSM
    return aSM+p+bSM+rootSM
def ppSquareRoot(rad):
    radSM = pp(rad)
    linedSM = StringMatrix(width=radSM.width,height=radSM.height+1)
    linedSM.setSubMatrix(radSM, 1,0)
    for x in range(radSM.width):
        linedSM.setChar(0, x, LOW_LINE)
    linedSM.mainRow = radSM.mainRow+1
    diag = UpDiag(radSM.height)
    diag.mainRow=linedSM.mainRow-1
    sm = diag+linedSM
    col = StringMatrix(width=1,height=sm.height)
    col.mainRow = sm.mainRow
    col.setChar(sm.height-1, 0, DIAG_DOWN)
    
    return col+sm


def ppMonomial(coeff, base, power,customVar=None):
    coeffSM = pp(coeff)
    if power==0:
        return coeffSM
    if coeffSM.height>1 and not coeffSM.hasBrackets:
        coeffSM.surroundWithBrackets()
    
    powerSM = ppPower(base,power)
    if coeff==1:
        return powerSM
    return coeffSM + powerSM

def ppPower(base, power,customVar=None):
    if customVar!=None:
        baseSM = StringMatrix.fromString(customVar)
    else:
        baseSM = pp(base)
    if power==1:
        return baseSM
    if baseSM.height>1 and not baseSM.hasBrackets:
        baseSM.surroundWithBrackets()
    powerSM = pp(power)
    if powerSM.height>1 and not powerSM.hasBrackets:
        powerSM.surroundWithBrackets()
        
    rowsNeeded = powerSM.height
    rowsToAdd = rowsNeeded-baseSM.mainRow
    sm = StringMatrix(height=baseSM.height+rowsToAdd,width=baseSM.width+powerSM.width)
    sm.setSubMatrix(baseSM,rowsToAdd,0)
    sm.setSubMatrix(powerSM,0,baseSM.width)
    sm.mainRow = rowsToAdd+baseSM.mainRow
    return sm

def ppIntegralExpression(func, var):
    funcSM = pp(func)
    funcSM.insertColumn(0,IntegralSymbol(funcSM.height))
    varSM = StringMatrix.fromString(" d"+str(var))
    return funcSM+varSM

def ppEquation(left, right):
    eqSM = StringMatrix.fromString(" = ")
    return left+eqSM+right

    

if __name__ == "__main__":
    from Parse import parseExpressionFromStr
    [T1,T2] = FE.Variables(['T_1','T_2'])
    x = FE.getVariable("x")
    fieldExtension1 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,Number.ONE]),T1) # field extension with log(x)
    FE.fieldTower = FE.FieldTower(fieldExtensions=[fieldExtension1])
    fieldExtension2 = FE.FieldExtension(FE.TRANS_LOG,Pol.Polynomial([0,Number.ONE],variable=T1),T2) # field extension with log(log(x))
    FE.fieldTower = FE.fieldTower.addFieldExtension(fieldExtension2)
    
    
    sm = ppRational(Rational(5,3))
    sm.pprint()
    sm = ppRational(Rational(123,7))
    sm.surroundWithBrackets()
    sm.pprint()
    sm = ppPower(Rational(5,3),45)
    sm.pprint()
    
    sm = ppMonomial(Rational(100,3), x, Rational(42,5))
    sm.pprint()
    
    poly = Pol.Polynomial([0,Pol.Polynomial([0,1],x),Rational(5,2)],variable=T1)
    sm = ppPolynomial(poly)
    sm.pprint()
    
    rat = parseExpressionFromStr("(x+1)/(x**2+1)*T_1/T_2+T_1*x**2", fieldTower=FE.fieldTower)
    sm = pp(rat)
    sm.pprint()
    
    ppSquareRoot(-5).pprint()
    ppSquareRoot(rat).pprint()