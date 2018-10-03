'''
Created on 28.09.2018

@author: Leonard
'''
from Utils import isNumber

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
        
    def subMatrix(self, i, j):
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
    
    def determinant(self):
        return self.determinantLaplace()
    def determinantGauss(self):
        pass
    def determinantLaplace(self): # determinant calculation using laplace theorem, complexity = O(n!) (n=size)
        if self.size==1:
            return self.getElement(0,0)
        s = 0
        for j in range(self.size):
            s += (-1)**j*self.getElement(0, j)*self.subMatrix(0,j).determinantLaplace()
            
        return s
    
    def getColumn(self, j):
        col = []
        for i in range(self.size):
            col.append(self.getElement(i, j))
        return col
    def getRow(self,i):
        row = []
        for j in range(self.size):
            row.append(self.getElement(i, j)) 
        return row
    
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

class SylvesterMatrix(QuadMatrix):
    def __init__(self, coeffsA,coeffsB):
        self.degA = len(coeffsA)-1
        self.degB = len(coeffsB)-1
        super(SylvesterMatrix,self).__init__(self.degA+self.degB)
        self.coeffsA = coeffsA
        self.coeffsB = coeffsB
        cA = reverseList(coeffsA)
        cB = reverseList(coeffsB)
        for j in range(self.degA+1):
            for i in range(self.degB):
                self.setElement(i, i+j, cA[j])
        for j in range(self.degB+1):
            for i in range(self.degA):
                self.setElement(i+self.degA-1,j+i,cB[j])
        
        
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
    return SylvesterMatrix(coeffsA,coeffsB).determinant()
            
        
def reverseList(a):
    b = []
    for i in range(len(a)-1,-1,-1):
        b.append(a[i])
    return b
if __name__ == '__main__':
    M = QuadMatrix(4)
    M.setElement(2,3,42)
    M.setElement(3,1,-1)
    M.setElement(3,3,2)
    print(M)
    print(M.subMatrix(2, 3))
    print("-----")
    A = QuadMatrix(2,data=[[3,2],[43,101]])
    print(A.getElement(1,0))
    print(A)
    print(A.determinant())
    B = QuadMatrix(3,data=[[-2,2,-3],[-1,1,3],[2,0,1]])
    print(B.determinant())
    C = SylvesterMatrix([-2,-1,1,3,3],[5,1,-3,1])
    print(C)
    print(Resultant([-2,-1,1,3,3],[5,1,-3,1]))