'''
Created on 28.09.2018

@author: Leonard
'''
from __future__ import division
from Utils import *
import RationalFunction as Rat

class QuadMatrix(object):
    """
    Matrix class
    """

    def __init__(self, size, data=None):
        self.size = size
        if data==None:
            self.__data = [[0 for i in range(size)] for j in range(size)]
        else:
            self.__data = data
        
    def getSubMatrix(self, i, j):
        M = QuadMatrix(self.size-1)
        for n in range(self.size-1):
            for m in range(self.size-1):
                a = n
                if n>=i:
                    a = n+1
                b = m
                if m>=j:
                    b = m+1
                
                M.setElement(n,m,self.getElement(a,b))
        return M
    #def setSubMatrix(self, i,j, matrix):
        
    
    def determinant(self):
        return self.determinantGauss()#self.determinantLaplace()
    def determinantGauss(self):
        if self.size==1:
            return self.getElement(0, 0)
        factor = 1
        #i_max  = argmax(list(map(abs,[self.getElement(i, 0) for i in range(self.size)])))
        i_max = 0
        for i in range(self.size):
            if self.getElement(i,0)!=0:
                i_max = i
                break
        M = self.copy()
        if self.getElement(i_max, 0)==0:
            return 0
        else:
            M.swapRows(0,i_max)
            if i_max!=0:
                factor *= (-1)
            subM = M.getSubMatrix(0, 0)
            mainRow = M.getRow(0)[1:self.size]
            lcoeff = M.getElement(0, 0)
            for i in range(1,self.size):
                f = -M.getElement(i, 0)/lcoeff
                if type(f) is Rat.RationalFunction:
                    f.makeDenominatorMonicInPlace()
                s = mulObjectToListElements(f, mainRow)
                newRow = addListToList(subM.getRow(i-1), s)
                subM.setRow(i-1,newRow)
            r = M.getElement(0, 0)*subM.determinantGauss()*factor
            if type(r)is Rat.RationalFunction:
                p = r.asPolynomial()
                if p!=None:
                    return p
            return r
            #(subTri, vecn) = subM.getTriangularMatrix(vec[1:self.size])
            
            
    def determinantLaplace(self): # determinant calculation using laplace theorem, complexity = O(n!) (n=size)
        if self.size==1:
            return self.getElement(0,0)
        s = 0
        for j in range(self.size):
            s += (-1)**j*self.getElement(0, j)*self.getSubMatrix(0,j).determinantLaplace()
            
        return s
    
    def getTriangularMatrix(self,vec):
        """
        converts matrix A=
        ( a11 a12 ... a1n )
        ( a21 a22 ... a2n )
              .
              .
              .
        ( an1 an2 ... ann )
        into a triangular matrix:
        ( b11 b12 ... b1n )
        ( 0   b22 ... b2n )
        ( 0   0 b33...b3n )
              .
              .
              .
        ( 0   0 ... 0 bnn )
        
        can be used for solving linear system of equations
        """
        if self.size==1:
            return (self, vec)
        M = self.copy()
        subM = self.getSubMatrix(0, 0)
        i_max  = argmax(list(map(abs,[self.getElement(i, 0) for i in range(self.size)])))
        if self.getElement(i_max, 0)==0:
            (subTri, vecn) = subM.getTriangularMatrix(vec[1:self.size])
            vec = [vec[0]]+vecn
        else:
            M.swapRows(0,i_max)
            t = vec[i_max]
            vec[i_max] = vec[0]
            vec[0] = t
            subM = M.getSubMatrix(0, 0)
            mainRow = M.getRow(0)[1:self.size]
            lcoeff = M.getElement(0, 0)
            for i in range(1,self.size):
                f = -M.getElement(i, 0)/lcoeff
                s = mulObjectToListElements(f, mainRow)
                newRow = addListToList(subM.getRow(i-1), s)
                subM.setRow(i-1,newRow)
                vec[i] += vec[0]*f
            (subTri, vecn) = subM.getTriangularMatrix(vec[1:self.size])
            vec = [vec[0]]+vecn
            
        for i in range(self.size-1):
            for j in range(self.size-1):
                M.setElement(i+1, j+1, subTri.getElement(i,j))
            M.setElement(i+1,0,0)
        return (M, vec)                
        
    def getColumn(self, j):
        col = []
        for i in range(self.size):
            col.append(self.getElement(i, j))
        return col
    def setColumn(self, j, col):
        if len(col)!=self.size:
            raise TypeError("input column has to have the same length as the matrix size.")
        for i in range(self.size):
            self.setElement(i, j, col[i])
            
    def getRow(self,i):
        row = []
        for j in range(self.size):
            row.append(self.getElement(i, j)) 
        return row
    def setRow(self, i, row):
        if len(row)!=self.size:
            raise TypeError("input row has to have the same length as the matrix size.")
        for j in range(self.size):
            self.setElement(i, j, row[j])
    def swapRows(self, a,b):
        tr = self.getRow(a)
        self.setRow(a, self.getRow(b))
        self.setRow(b, tr)
        
    def copy(self):
        data = [[x for x in self.__data[i]] for i in range(self.size)]
        return QuadMatrix(self.size,data=data)
        
    def getElement(self, i,j): # i=row, j = column
        return self.__data[i][j]
    def setElement(self, i,j,value):
        self.__data[i][j] = value
    def __str__(self):
        out = ""
        for i in range(self.size):
            r = "["
            for j in range(self.size):
                r += str(self.getElement(i,j))+","
            r = r.strip(",")+"]"
            out += r+"\n"
        return out.strip("\n")
    def __repr__(self):
        return self.__str__()

class SylvesterMatrix(QuadMatrix):
    def __init__(self, coeffsA,coeffsB):
        self.degA = len(coeffsA)-1
        self.degB = len(coeffsB)-1
        super(SylvesterMatrix,self).__init__(self.degA+self.degB)
        self.coeffsA = coeffsA
        self.coeffsB = coeffsB
        cA = _reverseList(coeffsA)
        cB = _reverseList(coeffsB)
        for j in range(self.degA+1):
            for i in range(self.degB):
                self.setElement(i, i+j, cA[j])
        for j in range(self.degB+1):
            for i in range(self.degA):
                self.setElement(i+self.degB,j+i,cB[j])
        
        
def Resultant(coeffsA,coeffsB):
    if len(coeffsA)==1:
        if isNumber(coeffsA[0]):
            if coeffsA[0]==0:
                return 0
            else:
                if len(coeffsB)==0:
                    if isNumber(coeffsB[0]):
                        if coeffsB[0]!=0:
                            return 1
                    elif not coeffsB[0].isZero():
                        return 1
                    
        else:
            if coeffsA[0].isZero():
                return 0
            else:
                if len(coeffsB)==0:
                    if isNumber(coeffsB[0]):
                        if coeffsB[0]!=0:
                            return 1
                    elif not coeffsB[0].isZero():
                        return 1
                    
    if len(coeffsB)==1:
        if isNumber(coeffsB[0]):
            if isNumber==0:
                return 0
        elif coeffsB[0].isZero():
            return 0
        
    if len(coeffsA)==0 or len(coeffsB)==0:
        raise Exception()
    sylv = SylvesterMatrix(coeffsA,coeffsB)
    res = timeMethod(sylv.determinant)
    if type(res) is Rat.RationalFunction:
        r = res.asPolynomial()
        if r!=None:
            return r
    return res
            
def solveLinearSystem(coeffs, vec):
    (tri, vec) = coeffs.getTriangularMatrix(vec)
    n = coeffs.size
    xs = []
    if tri.getElement(n-1,n-1)==0 and vec[n-1]!=0:
        return None # last row : 0 0 0 ... 0 | x, x!=0 => no solution
    unique = True
    for i in range(n):
        if tri.getElement(n-1-i,n-1-i)==0:
            xs.append(1)
            unique = False
            continue
        s = 0
        for j in range(i):
            s += xs[j]*tri.getElement(n-1-i,n-1-j)
        d = (vec[n-1-i]+(-s))
        xs.append(d/tri.getElement(n-1-i,n-1-i))
    return (_reverseList(xs), unique)
    

def _reverseList(a):
    b = []
    for i in range(len(a)-1,-1,-1):
        b.append(a[i])
    return b

def matrixtest(expected, got):
    print("matrixtest: should be {}, got {}".got(expected,got))
if __name__ == '__main__':
    from Number import Rational
    """ M = QuadMatrix(4)
    M.setElement(2,3,42)
    M.setElement(3,1,-1)
    M.setElement(3,3,2)
    M.setRow(0, [1,2,3,4])
    M.setColumn(1, [101,102,1000,-1])
    M.swapRows(1,2)
    print(M)
    print(M.copy())
    print(M.getSubMatrix(2, 3))
    print("-----")"""
    A = QuadMatrix(2,data=[[3,2],[43,101]])
    #print(A)
    print(A.determinant()) # 217
    B = QuadMatrix(3,data=[[-2,2,-3],[-1,1,3],[2,0,1]])
    print(B.determinant()) # 18
    C = SylvesterMatrix([-2,-1,1,3,3],[5,1,-3,1])
    #print(C)
    print(Resultant([-2,-1,1,3,3],[5,1,-3,1])) # 0
    
    M = QuadMatrix(3,data=[[2,1,-1],[-3,-1,2],[-2,1,2]])
    #print(M.getTriangularMatrix([8,-11,-3]))
    print(solveLinearSystem(M, [8,-11,-3]))
    
    M = QuadMatrix(2,data=[[Rational(1,1),Rational(1,2)],[Rational(1,2),Rational(1,3)]])
    print(solveLinearSystem(M, [Rational(1,1),Rational(5,7)]))
    #print(M.getTriangularMatrix([Rational(1,1),Rational(5,7)]))
    
    M = QuadMatrix(2,data=[[1,1],[2,2]])
    print(solveLinearSystem(M, [1,2]))
    
    M = QuadMatrix(4,data=[[3,5,-1,0],[0,2,0,3],[-1,-2,0,1],[3,3,1,2]])
    print(M.determinant())