'''
Created on 28.09.2018

@author: Leonard
'''
from __future__ import division
from Utils import *
import RationalFunction as Rat
from Number import Rational,ONE,ZERO


class Matrix(object):
    def __init__(self, rows=1, columns=1,data=None):
        if data==None:
            self._data = [[ZERO for _ in range(columns)] for _ in range(rows)]
            self.rows = rows
            self.columns = columns
        else:
            self._data = data
            self.rows = len(data)
            self.columns = len(data[0])
        
        
    @property
    def dimension(self):
        return (self.rows,self.columns)
    
    def replaceNumbersWithRationals(self):
        for i in range(self.rows):
            for j in range(self.columns):
                if type(self[i,j])==float or type(self[i,j])==int:
                    self[i,j] = Rational.fromFloat(self[i,j])
    
    def getSubMatrix(self,i,j):
        M = Matrix(self.rows-1,self.columns-1)
        for n in range(self.rows-1):
            for m in range(self.columns-1):
                a = n
                if n>=i:
                    a = n+1
                b = m
                if m>=j:
                    b = m+1
                
                M[n,m] = self[a,b]
        return M
    
    # ========================================== set/get elements =========================================
    def getElement(self, i,j):
        return self._data[i][j]
    def setElement(self, i,j,item):
        self._data[i][j] = item
    def __getitem__(self, index):
        if type(index)!=tuple:
            return self.getRow(index)
        return self.getElement(index[0], index[1])
    def __setitem__(self, index, item):
        self._data[index[0]][index[1]] = item
        
    # ========================================== set/get/adjust rows/columns =========================================
    def getColumn(self, j):
        return [self[i,j] for i in range(self.rows)]
    def setColumn(self, j, col):
        if len(col)!=self.rows:
            raise TypeError()
        for i in range(self.rows):
            self[i,j]=col[i]
            
    def getRow(self,i):
        return [self[i,j] for j in range(self.columns)]
    def setRow(self, i, row):
        if len(row)!=self.columns:
            raise TypeError()
        for j in range(self.columns):
            self[i, j] = row[j]
    
    def rowIsZero(self,i):
        for j in range(self.columns):
            if self[i,j]!=0:
                return False
        return True
    def swapRows(self, a,b):
        tr = self.getRow(a)
        self.setRow(a, self.getRow(b))
        self.setRow(b, tr)
  
    # ========================================== some matrix stuff =========================================
    def determinant(self):
        if self.rows!=self.columns:
            raise ValueError("Can't calculate the determinant of a non-square matrix")
        sqm = QuadMatrix(self.rows,data=self._data)
        return sqm.determinant()
    def getRowEchelonForm(self):
        """
        converts matrix A=
        ( a11 a12 ... a1m )
        ( a21 a22 ... a2m )
              .
              .
              .
        ( an1 an2 ... anm )
        into a triangular matrix:
        ( b11 b12 ... b1m )
        ( 0   b22 ... b2m )
        ( 0   0 b33...b3m )
              .
              .
              .
        ( 0   0 ... 0 bnm )
        can be used for solving linear system of equations
        """
        #if self.rows>=self.columns:
        #    raise NotImplementedError()
        if self.columns==1 or self.rows==1:
            return self
        M = self.copy()
        subM = self.getSubMatrix(0, 0)
        i_max  = argmax(list(map(abs,[self[i,0] for i in range(self.rows)])))
        if self.getElement(i_max, 0)==0:
            subTri = subM.getRowEchelonForm()
        else:
            M.swapRows(0,i_max)
            subM = M.getSubMatrix(0, 0)
            mainRow = M.getRow(0)[1:self.columns]
            lcoeff = M.getElement(0, 0)
            for i in range(1,self.rows):
                f = -M[i,0]/lcoeff
                s = mulObjectToListElements(f, mainRow)
                newRow = addListToList(subM.getRow(i-1), s)
                subM.setRow(i-1,newRow)
            subTri = subM.getRowEchelonForm()
            
        for i in range(self.rows-1):
            for j in range(self.columns-1):
                M.setElement(i+1, j+1, subTri.getElement(i,j))
            M.setElement(i+1,0,0)
        return M
    
    def getReducedRowEchelonForm(self):
        tri = self.getRowEchelonForm()
        for i in range(self.rows-1,-1,-1):
            if tri[i,i]!=0:
                if i>0:
                    for ii in range(i):
                        f = -tri[ii,i]/tri[i,i]
                        s = mulObjectToListElements(f, tri.getRow(i))
                        tri.setRow(ii,addListToList(s, tri.getRow(ii)))
                tri.setRow(i, mulObjectToListElements(1/tri[i,i], tri.getRow(i)))
        return tri
        
    def copy(self):
        data = [[x for x in self._data[i]] for i in range(self.rows)]
        return Matrix(self.rows,self.columns,data=data)
    
    def __add__(self, other):
        if not issubclass(type(other), Matrix):
            raise TypeError()
        if self.dimension!=other.dimension:
            raise TypeError("The dimensions of the matrices have to be identical")
        
    def __str__(self):
        out = ""
        for i in range(self.rows):
            r = "["
            for j in range(self.columns):
                r += str(self.getElement(i,j))+","
            r = r.strip(",")+"]"
            out += r+"\n"
        return out.strip("\n")
    def __repr__(self):
        return self.__str__()
        
        
class QuadMatrix(Matrix):
    """
    Matrix class
    """

    def __init__(self, size, data=None):
        self.size = size
        self.rows = size
        self.columns = size
        if data==None:
            self._data = [[0 for _ in range(size)] for _ in range(size)]
        else:
            self._data = data
           
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
            r = M.getElement(0, 0)*subM.determinant()*factor
            if type(r)is Rat.RationalFunction:
                p = r.asPolynomial()
                if p!=None:
                    return p
            return r
            
            
    def determinantLaplace(self): # determinant calculation using laplace theorem, complexity = O(n!) (n=size)
        if self.size==1:
            return self.getElement(0,0)
        s = 0
        for j in range(self.size):
            s += (-1)**j*self.getElement(0, j)*self.getSubMatrix(0,j).determinantLaplace()
            
        return s
    
        
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
        cA = list(reversed(coeffsA))
        cB = list(reversed(coeffsB))
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
    """
    Equation: coeffs * (x0,x1,...,xn)=(v0,v1,...,vn)
    returns None if the equation has no solution and a 4-tuple if it has. The tuple contains a particular solution, uniqueness of the solution, general solution, dimension of the solution
    """
    
    augmentedMatrix = Matrix(rows=coeffs.rows,columns=coeffs.columns+1)
    for i in range(coeffs.rows):
        for j in range(coeffs.columns):
            augmentedMatrix[i,j] = coeffs[i,j]
    augmentedMatrix.setColumn(coeffs.columns, vec)
    
    rref = augmentedMatrix.getReducedRowEchelonForm()
    vec = rref.getColumn(rref.columns-1)
    
    for i in range(rref.rows-1,-1,-1):
        if any(rref[i,j]!=0 for j in range(rref.columns-1)):
            break
        if vec[i]!=0:
            return None
    
    xs = [None]*coeffs.columns
    
    # solve for all variables, that are uniquely determined, i.e. their corresponding row is in rref is (0,...,0,1,0,...|a)
    for i in range(coeffs.rows):
        oneVar = False
        varIndex = 0
        for j in range(coeffs.columns):
            if rref[i,j]!=0:
                if oneVar:
                    oneVar = False
                    break
                else:
                    varIndex = j
                    oneVar = True
        if oneVar:
            xs[varIndex] = lambda _,varIndex=varIndex,i=i: vec[i]
            
            
    # determine indices of variables that can be set arbitrarily (independent from each other, note: such a set is not uniquely determined but that doesn't matter)
    # the values of the other variables can then be computed with those. (such
    #arb_vars = []
    non_arb_vars = []
    for i in range(rref.rows):
        row = rref.getRow(i)
        non_arb_vars.append(min([rref.columns]+[k for k in range(coeffs.columns) if row[k]!=0]))
    arb_vars = [k for k in range(coeffs.columns) if k not in non_arb_vars]
    
    # the solution space of a system of linear equation is a 'translated' sub-vectorspace of R^n, ( it's translated, because the system might be inhomogeneous -> the zero vector is no solution)
    # it's dimension equals to the maxmimal amount of variables that can be set arbitrarily
    dimension = len(arb_vars)
    
    for k in range(dimension):
        xs[arb_vars[k]] = lambda v,k=k: v[k]
        
    for i in range(rref.rows):
        row = rref.getRow(i)
        non_zero_args = [k for k in range(coeffs.columns) if row[k]!=0]
        if len(non_zero_args)>0:
            j = min(non_zero_args)
            if xs[j]==None:
                xs[j] = lambda v,i=i: (vec[i]-sum(rref[i,arb_vars[k]]*v[k] for k in range(dimension)))
    
    general_solution = lambda v: [xs[i](v) for i in range(coeffs.columns)]
    zeros = [0]*dimension
    
    if dimension==0:
        unique = True
    else:
        unique = False
    
    return (general_solution(zeros), unique, general_solution, dimension) 
  
    
def matrixtest(expected, got):
    print("matrixtest: should be {}, got {}".got(expected,got))
if __name__ == '__main__':
    M = Matrix(2,3)
    M[1,2] = 5
    M[0,1] = 42
    print(M.getRow(1))
    print(M.getColumn(1))
    print(M)
    print("----")
    M = Matrix(3,4)
    M.setColumn(2,[1,2,3])
    M.setRow(1,[42,-1,2,3])
    print(M)
    print("--")
    print(M.getSubMatrix(2,2))
    print("-----")
    M = Matrix(data=[[2,1,-1,8],[-3,-1,2,-11],[-2,1,2,-3]])
    M.replaceNumbersWithRationals()
    print(M.getReducedRowEchelonForm())
    print("----")
    M = Matrix(data=[[2,5,2,-38],[3,-2,4,17],[-6,1,-7,-12]])
    print(M.getReducedRowEchelonForm())
    print("---")
    M = Matrix(data=[[2,2,0],[0,1,0]])
    print(M.getReducedRowEchelonForm())
    M = Matrix(data=[[2,2,2,0],[0,0,1,0],[1,1,0,0]])
    M.replaceNumbersWithRationals()
    print(M.getReducedRowEchelonForm())
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
    #print(A.determinant()) # 217
    B = QuadMatrix(3,data=[[-2,2,-3],[-1,1,3],[2,0,1]])
    #print(B.determinant()) # 18
    C = SylvesterMatrix([-2,-1,1,3,3],[5,1,-3,1])
    #print(C)
    #print(Resultant([-2,-1,1,3,3],[5,1,-3,1])) # 0
    
    M = QuadMatrix(3,data=[[2,0,2],[0,1,-2],[0,0,0]])
    
    M = QuadMatrix(3,data=[[2*ONE,ONE,-ONE],[-3*ONE,-1*ONE,2*ONE],[-2*ONE,1*ONE,2*ONE]])
    #print(M.getTriangularMatrix([8,-11,-3]))
    print(solveLinearSystem(M, [8*ONE,-11*ONE,-3*ONE]))
    
    M = QuadMatrix(2,data=[[Rational(1,1),Rational(1,2)],[Rational(1,2),Rational(1,3)]])
    print(solveLinearSystem(M, [Rational(1,1),Rational(5,7)]))
    #print(M.getTriangularMatrix([Rational(1,1),Rational(5,7)]))
    
    M = QuadMatrix(2,data=[[1,1],[2,2]])
    print(solveLinearSystem(M, [1,2]))
    
    M = QuadMatrix(4,data=[[3,5,-1,0],[0,2,0,3],[-1,-2,0,1],[3,3,1,2]])
    print(M.determinant())